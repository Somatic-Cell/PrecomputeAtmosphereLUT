#include "bake_common.cuh"
#include "device_common.cuh"

#include <cuda_runtime.h>

#include <iostream>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace bake {
namespace {

void checkCuda(cudaError_t err, const char* what) {
    if (err != cudaSuccess) {
        throw std::runtime_error(std::string(what) + ": " + cudaGetErrorString(err));
    }
}

struct DeviceSpectralData {
    float* solarIrradiance      = nullptr;
    float* rayleighScattering   = nullptr;
    float* mieScattering        = nullptr;
    float* mieExtinction        = nullptr;
    float* absorptionExtinction = nullptr;
};

void freeDeviceSpectralData(DeviceSpectralData& d) {
    cudaFree(d.solarIrradiance);
    cudaFree(d.rayleighScattering);
    cudaFree(d.mieScattering);
    cudaFree(d.mieExtinction);
    cudaFree(d.absorptionExtinction);
    d = {};
}

DeviceSpectralData uploadDeviceSpectralData(const atmo::AtmosphereMetadata& meta) {
    DeviceSpectralData d{};
    const std::size_t bytes = meta.wavelengths_nm.size() * sizeof(float);

    auto upload = [bytes](float*& dst,
                          const std::vector<float>& src,
                          const char* name) {
        checkCuda(cudaMalloc(&dst, bytes), name);
        checkCuda(cudaMemcpy(dst, src.data(), bytes, cudaMemcpyHostToDevice), name);
    };

    try {
        upload(d.solarIrradiance,      meta.solarIrradiance,      "cudaMemcpy solarIrradiance");
        upload(d.rayleighScattering,   meta.rayleighScattering,   "cudaMemcpy rayleighScattering");
        upload(d.mieScattering,        meta.mieScattering,        "cudaMemcpy mieScattering");
        upload(d.mieExtinction,        meta.mieExtinction,        "cudaMemcpy mieExtinction");
        upload(d.absorptionExtinction, meta.absorptionExtinction, "cudaMemcpy absorptionExtinction");
    } catch (...) {
        freeDeviceSpectralData(d);
        throw;
    }

    return d;
}

ATMO_HD inline float decodeScatteringMuFromUnit(float u) {
    // internal scattering は全方向を持つ
    return -1.0f + 2.0f * atmo::clampf(u, 0.0f, 1.0f);
}

ATMO_HD inline float decodeMuSFromUnit(float u, float muSMin) {
    return muSMin + atmo::clampf(u, 0.0f, 1.0f) * (1.0f - muSMin);
}

ATMO_HD inline float decodeNuFromUnit(float u) {
    return -1.0f + 2.0f * atmo::clampf(u, 0.0f, 1.0f);
}

ATMO_HD inline float radiusFromScatteringV(
    const DeviceAtmosphereParameters& atm,
    float v) {

    const float bottom = atm.bottomRadius_m;
    const float top    = atm.topRadius_m;

    const float H   = atmo::safeSqrt(top * top - bottom * bottom);
    const float rho = H * atmo::clampf(v, 0.0f, 1.0f);

    return atmo::safeSqrt(rho * rho + bottom * bottom);
}

ATMO_HD inline float lerpfLocal(float a, float b, float t) {
    return a + (b - a) * t;
}

ATMO_HD inline float encodeScatteringMuToUnit(float mu) {
    return 0.5f * (atmo::clampCosine(mu) + 1.0f);
}


ATMO_HD inline float encodeMuSToUnit(float muS, float muSMin) {
    const float denom = 1.0f - muSMin;
    if (denom <= 0.0f) {
        return 0.0f;
    }
    return atmo::clampf((muS - muSMin) / denom, 0.0f, 1.0f);
}

ATMO_HD inline float encodeNuToUnit(float nu) {
    return 0.5f * (atmo::clampCosine(nu) + 1.0f);
}

ATMO_HD inline float encodeRadiusToUnit(
    const DeviceAtmosphereParameters& atm,
    float r_m) {

    const float bottom = atm.bottomRadius_m;
    const float top    = atm.topRadius_m;

    const float H = atmo::safeSqrt(top * top - bottom * bottom);
    if (H <= 0.0f) {
        return 0.0f;
    }

    const float rr  = atmo::clampf(r_m, bottom, top);
    const float rho = atmo::safeSqrt(rr * rr - bottom * bottom);
    return atmo::clampf(rho / H, 0.0f, 1.0f);
}

// 同じ (r, muS, lambda) スライス上で、(mu, nu) を bilinear sample
ATMO_HD inline float sampleScatteringSameSlice(
    const float* table,
    const DeviceLutDims& dims,
    int iR,
    int iMuS,
    int iLambda,
    float mu,
    float nu) {

    const float uMu = encodeScatteringMuToUnit(mu);
    const float uNu = encodeNuToUnit(nu);

    const float xMu = uMu * float(dims.skyMu - 1);
    const float xNu = uNu * float(dims.skyNu - 1);

    const int mu0 = clampIndex(static_cast<int>(xMu), dims.skyMu);
    const int nu0 = clampIndex(static_cast<int>(xNu), dims.skyNu);
    const int mu1 = clampIndex(mu0 + 1, dims.skyMu);
    const int nu1 = clampIndex(nu0 + 1, dims.skyNu);

    const float tMu = xMu - float(mu0);
    const float tNu = xNu - float(nu0);

    const float f00 = table[scatteringIndex(iR, mu0, iMuS, nu0, iLambda, dims)];
    const float f10 = table[scatteringIndex(iR, mu1, iMuS, nu0, iLambda, dims)];
    const float f01 = table[scatteringIndex(iR, mu0, iMuS, nu1, iLambda, dims)];
    const float f11 = table[scatteringIndex(iR, mu1, iMuS, nu1, iLambda, dims)];

    const float fx0 = lerpfLocal(f00, f10, tMu);
    const float fx1 = lerpfLocal(f01, f11, tMu);
    return lerpfLocal(fx0, fx1, tNu);
}

// 4 次元 (r, mu, muS, nu) + lambda 固定の 5D LUT を quadlinear sample
ATMO_HD inline float sampleScatteringLUT(
    const float* table,
    const DeviceLutDims& dims,
    const DeviceAtmosphereParameters& atm,
    int iLambda,
    float r_m,
    float mu,
    float muS,
    float nu) {

    const float uR   = encodeRadiusToUnit(atm, r_m);
    const float uMu  = encodeScatteringMuToUnit(mu);
    const float uMuS = encodeMuSToUnit(muS, atm.muSMin);
    const float uNu  = encodeNuToUnit(nu);

    const float xR   = uR   * float(dims.scatteringR - 1);
    const float xMu  = uMu  * float(dims.skyMu - 1);
    const float xMuS = uMuS * float(dims.skyMuS - 1);
    const float xNu  = uNu  * float(dims.skyNu - 1);

    const int r0   = clampIndex(static_cast<int>(xR),   dims.scatteringR);
    const int mu0  = clampIndex(static_cast<int>(xMu),  dims.skyMu);
    const int muS0 = clampIndex(static_cast<int>(xMuS), dims.skyMuS);
    const int nu0  = clampIndex(static_cast<int>(xNu),  dims.skyNu);

    const int r1   = clampIndex(r0 + 1,   dims.scatteringR);
    const int mu1  = clampIndex(mu0 + 1,  dims.skyMu);
    const int muS1 = clampIndex(muS0 + 1, dims.skyMuS);
    const int nu1  = clampIndex(nu0 + 1,  dims.skyNu);

    const float tR   = xR   - float(r0);
    const float tMu  = xMu  - float(mu0);
    const float tMuS = xMuS - float(muS0);
    const float tNu  = xNu  - float(nu0);

    float accum = 0.0f;

    for (int br = 0; br < 2; ++br) {
        const int ir = br ? r1 : r0;
        const float wr = br ? tR : (1.0f - tR);

        for (int bm = 0; bm < 2; ++bm) {
            const int imu = bm ? mu1 : mu0;
            const float wmu = bm ? tMu : (1.0f - tMu);

            for (int bs = 0; bs < 2; ++bs) {
                const int imuS = bs ? muS1 : muS0;
                const float wmuS = bs ? tMuS : (1.0f - tMuS);

                for (int bn = 0; bn < 2; ++bn) {
                    const int inu = bn ? nu1 : nu0;
                    const float wnu = bn ? tNu : (1.0f - tNu);

                    const float w = wr * wmu * wmuS * wnu;
                    accum += w * table[scatteringIndex(ir, imu, imuS, inu, iLambda, dims)];
                }
            }
        }
    }

    return accum;
}

struct ScatteringFlatIndex {
    int iR;
    int iMu;
    int iMuS;
    int iNu;
    int iLambda;
};

ATMO_HD inline std::size_t scatteringElementCount(const DeviceLutDims& dims) {
    return static_cast<std::size_t>(dims.scatteringR) *
           static_cast<std::size_t>(dims.skyMu) *
           static_cast<std::size_t>(dims.skyMuS) *
           static_cast<std::size_t>(dims.skyNu) *
           static_cast<std::size_t>(dims.wavelengthCount);
}

ATMO_HD inline ScatteringFlatIndex decodeScatteringFlatIndex(
    std::size_t flat,
    const DeviceLutDims& dims) {

    ScatteringFlatIndex out{};

    out.iLambda = static_cast<int>(flat % dims.wavelengthCount);
    flat /= dims.wavelengthCount;

    out.iNu = static_cast<int>(flat % dims.skyNu);
    flat /= dims.skyNu;

    out.iMuS = static_cast<int>(flat % dims.skyMuS);
    flat /= dims.skyMuS;

    out.iMu = static_cast<int>(flat % dims.skyMu);
    flat /= dims.skyMu;

    out.iR = static_cast<int>(flat);
    return out;
}

inline dim3 make1DGrid(std::size_t total, int blockSize) {
    const std::size_t gridX = (total + static_cast<std::size_t>(blockSize) - 1) /
                              static_cast<std::size_t>(blockSize);
    return dim3(static_cast<unsigned int>(gridX), 1u, 1u);
}

struct SkyFlatIndex {
    int iMu;
    int iMuS;
    int iNu;
    int iLambda;
};

ATMO_HD inline std::size_t skyElementCount(const DeviceLutDims& dims) {
    return static_cast<std::size_t>(dims.skyMu) *
           static_cast<std::size_t>(dims.skyMuS) *
           static_cast<std::size_t>(dims.skyNu) *
           static_cast<std::size_t>(dims.wavelengthCount);
}

ATMO_HD inline SkyFlatIndex decodeSkyFlatIndex(
    std::size_t flat,
    const DeviceLutDims& dims) {

    SkyFlatIndex out{};

    out.iLambda = static_cast<int>(flat % dims.wavelengthCount);
    flat /= dims.wavelengthCount;

    out.iNu = static_cast<int>(flat % dims.skyNu);
    flat /= dims.skyNu;

    out.iMuS = static_cast<int>(flat % dims.skyMuS);
    flat /= dims.skyMuS;

    out.iMu = static_cast<int>(flat);
    return out;
}

__global__ void kernelBuildDeltaSingleScattering(
    DeviceAtmosphereParameters atm,
    DeviceLutDims dims,
    const float* transmittanceLUT,
    std::uint32_t integrationStepsView,
    float* outDeltaRayleigh,
    float* outDeltaMie) {

    const std::size_t flat =
        static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;

    const std::size_t total =
        static_cast<std::size_t>(dims.scatteringR) *
        static_cast<std::size_t>(dims.skyMu) *
        static_cast<std::size_t>(dims.skyMuS) *
        static_cast<std::size_t>(dims.skyNu) *
        static_cast<std::size_t>(dims.wavelengthCount);

    if (flat >= total) {
        return;
    }

    const ScatteringFlatIndex idx = decodeScatteringFlatIndex(flat, dims);
    
    const int iR      = idx.iR;
    const int iMu     = idx.iMu;
    const int iMuS    = idx.iMuS;
    const int iNu     = idx.iNu;
    const int iLambda = idx.iLambda;

    const float uMu = (static_cast<float>(iMu) + 0.5f) / static_cast<float>(dims.skyMu);
    const float uMuS = (static_cast<float>(iMuS) + 0.5f) / static_cast<float>(dims.skyMuS);
    const float uNu = (static_cast<float>(iNu) + 0.5f) / static_cast<float>(dims.skyNu);
    const float vR = (static_cast<float>(iR) + 0.5f) / static_cast<float>(dims.scatteringR);

    const float mu  = decodeScatteringMuFromUnit(uMu);
    const float muS = decodeMuSFromUnit(uMuS, atm.muSMin);
    const float nu  = decodeNuFromUnit(uNu);
    const float r_m = radiusFromScatteringV(atm, vR);

    const std::size_t outIdx = scatteringIndex(iR, iMu, iMuS, iNu, iLambda, dims);

    // 視線が地面と交差するかどうかの判定
    const bool viewHitsGround = atmo::rayIntersectsGround(atm.bottomRadius_m, r_m, mu); 
  
    // 視線方向の積分区間の決定
    // 上向きなら top atmosphere まで，下向きなら ground hit まで積分
    const float sMax = atmo::distanceToNearestAtmosphereBoundary(
        atm.bottomRadius_m,
        atm.topRadius_m,
        r_m,
        mu,
        viewHitsGround);

    if (sMax <= 0.0f) {
        outDeltaRayleigh[outIdx] = 0.0f;
        outDeltaMie[outIdx] = 0.0f;
        return;
    }

    const float ds = sMax / static_cast<float>(integrationStepsView);

    float opticalDepthView = 0.0f;
    float accumRayleigh = 0.0f;
    float accumMie = 0.0f;

    for (std::uint32_t i = 0; i < integrationStepsView; ++i) {
        const float sMid = (static_cast<float>(i) + 0.5f) * ds;
        const float rSample = atmo::pointRadiusAlongRay(r_m, mu, sMid);
        const float altitude_m = rSample - atm.bottomRadius_m;

        const float sigmaExtMid = extinctionAt(atm, altitude_m, iLambda);
        const float tauToMid = opticalDepthView + 0.5f * sigmaExtMid * ds;
        const float tView = ::expf(-tauToMid);

        // 散乱点から見た太陽方向余弦
        const float muSSample =
            atmo::clampCosine((r_m * muS + sMid * nu) / rSample);

        const float tSun = sampleTransmittanceLUT(
            transmittanceLUT,
            dims,
            atm,
            iLambda,
            rSample,
            muSSample);

        const float rhoR = numberDensityRayleigh(atm, altitude_m);
        const float rhoM = numberDensityMie(atm, altitude_m);

        accumRayleigh += tView * tSun * rhoR * ds;
        accumMie      += tView * tSun * rhoM * ds;

        opticalDepthView += sigmaExtMid * ds;
    }

    // phase-free internal single scattering
    outDeltaRayleigh[outIdx] =
        accumRayleigh *
        atm.solarIrradiance[iLambda] *
        atm.rayleighScattering[iLambda];

    outDeltaMie[outIdx] =
        accumMie *
        atm.solarIrradiance[iLambda] *
        atm.mieScattering[iLambda];
}

__global__ void kernelAddSingleScattering(
    DeviceAtmosphereParameters atm,
    DeviceLutDims dims,
    const float* deltaRayleighSingle,
    const float* deltaMieSingle,
    float* outPrevScattering) {

    const std::size_t flat =
        static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;

    const std::size_t total = scatteringElementCount(dims);
    if (flat >= total) {
        return;
    }

    const ScatteringFlatIndex idx = decodeScatteringFlatIndex(flat, dims);

    const int iR      = idx.iR;
    const int iMu     = idx.iMu;
    const int iMuS    = idx.iMuS;
    const int iNu     = idx.iNu;
    const int iLambda = idx.iLambda;

    const float uNu = (static_cast<float>(iNu) + 0.5f) / static_cast<float>(dims.skyNu);
    const float nu  = decodeNuFromUnit(uNu);

    const std::size_t index = scatteringIndex(iR, iMu, iMuS, iNu, iLambda, dims);

    const float phaseR = atmo::rayleighPhaseFunction(nu);
    const float phaseM = atmo::miePhaseFunction(atm.miePhaseFunctionG, nu);

    outPrevScattering[index] =
        deltaRayleighSingle[index] * phaseR +
        deltaMieSingle[index] * phaseM;
}

__global__ void kernelBuildScatteringDensity(
    DeviceAtmosphereParameters atm,
    DeviceLutDims dims,
    const float* prevScattering,
    float* outScatteringDensity) {

    const std::size_t flat =
        static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;

    const std::size_t total = scatteringElementCount(dims);
    if (flat >= total) {
        return;
    }

    const ScatteringFlatIndex idx = decodeScatteringFlatIndex(flat, dims);

    const int iR      = idx.iR;
    const int iMu     = idx.iMu;
    const int iMuS    = idx.iMuS;
    const int iNu     = idx.iNu;
    const int iLambda = idx.iLambda;

    const float uMu  = (static_cast<float>(iMu) + 0.5f)     / static_cast<float>(dims.skyMu);
    const float uMuS = (static_cast<float>(iMuS) + 0.5f)    / static_cast<float>(dims.skyMuS);
    const float uNu  = (static_cast<float>(iNu) + 0.5f)     / static_cast<float>(dims.skyNu);
    const float vR   = (static_cast<float>(iR) + 0.5f)      / static_cast<float>(dims.scatteringR);

    const float mu   = decodeScatteringMuFromUnit(uMu);
    const float muS  = decodeMuSFromUnit(uMuS, atm.muSMin);
    const float nu   = decodeNuFromUnit(uNu);
    const float r_m  = radiusFromScatteringV(atm, vR);
    const float altitude_m = r_m - atm.bottomRadius_m;

    const float sigmaR = numberDensityRayleigh(atm, altitude_m) * atm.rayleighScattering[iLambda];
    const float sigmaM = numberDensityMie(atm, altitude_m)      * atm.mieScattering[iLambda];

    const float sinTheta  = atmo::safeSqrt(fmaxf(0.0f, 1.0f - mu  * mu));
    const float sinThetaS = atmo::safeSqrt(fmaxf(0.0f, 1.0f - muS * muS));

    float cosPhiOut = 1.0f;
    if (sinTheta > 1.0e-6f && sinThetaS > 1.0e-6f) {
        cosPhiOut = (nu - mu * muS) / (sinTheta * sinThetaS);
        cosPhiOut = atmo::clampf(cosPhiOut, -1.0f, 1.0f);
    }
    const float sinPhiOut = atmo::safeSqrt(fmaxf(0.0f, 1.0f - cosPhiOut * cosPhiOut));

    constexpr int kMuSamples  = 16;
    constexpr int kPhiSamples = 32;
    constexpr float kTwoPi = 6.28318530717958647692f;

    const float dMu  = 2.0f / float(kMuSamples);
    const float dPhi = kTwoPi / float(kPhiSamples);

    float accum = 0.0f;

    for (int a = 0; a < kMuSamples; ++a) {
        const float muIn = -1.0f + (float(a) + 0.5f) * dMu;
        const float sinThetaIn = atmo::safeSqrt(fmaxf(0.0f, 1.0f - muIn * muIn));

        for (int b = 0; b < kPhiSamples; ++b) {
            const float phiIn = (float(b) + 0.5f) * dPhi;
            const float cosPhiIn = cosf(phiIn);
            const float sinPhiIn = sinf(phiIn);

            const float nuIn =
                atmo::clampCosine(muIn * muS + sinThetaIn * sinThetaS * cosPhiIn);

            const float cosTheta =
                atmo::clampCosine(
                    mu * muIn +
                    sinTheta * sinThetaIn * (cosPhiOut * cosPhiIn + sinPhiOut * sinPhiIn));

            const float Lin = sampleScatteringSameSlice(
                prevScattering, dims, iR, iMuS, iLambda, muIn, nuIn);

            const float phaseR = atmo::rayleighPhaseFunction(cosTheta);
            const float phaseM = atmo::miePhaseFunction(atm.miePhaseFunctionG, cosTheta);

            const float sigmaScatter = sigmaR * phaseR + sigmaM * phaseM;

            accum += Lin * sigmaScatter * dMu * dPhi;
        }
    }

    outScatteringDensity[scatteringIndex(iR, iMu, iMuS, iNu, iLambda, dims)] = accum;
}

__global__ void kernelBuildDeltaMultipleScattering(
    DeviceAtmosphereParameters atm,
    DeviceLutDims dims,
    const float* transmittanceLUT,
    const float* scatteringDensity,
    std::uint32_t integrationStepsView,
    float* outDeltaMultiple) {

    const std::size_t flat =
        static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;

    const std::size_t total = scatteringElementCount(dims);
    if (flat >= total) {
        return;
    }

    const ScatteringFlatIndex idx = decodeScatteringFlatIndex(flat, dims);

    const int iR      = idx.iR;
    const int iMu     = idx.iMu;
    const int iMuS    = idx.iMuS;
    const int iNu     = idx.iNu;
    const int iLambda = idx.iLambda;

    const float uMu  = (static_cast<float>(iMu) + 0.5f) / static_cast<float>(dims.skyMu);
    const float uMuS = (static_cast<float>(iMuS) + 0.5f) / static_cast<float>(dims.skyMuS);
    const float uNu  = (static_cast<float>(iNu) + 0.5f) / static_cast<float>(dims.skyNu);
    const float vR   = (static_cast<float>(iR) + 0.5f) / static_cast<float>(dims.scatteringR);

    const float mu  = decodeScatteringMuFromUnit(uMu);
    const float muS = decodeMuSFromUnit(uMuS, atm.muSMin);
    const float nu  = decodeNuFromUnit(uNu);
    const float r_m = radiusFromScatteringV(atm, vR);

    const std::size_t outIdx = scatteringIndex(iR, iMu, iMuS, iNu, iLambda, dims);

    // 視線が地面と交差するかどうかの判定
    const bool viewHitsGround = atmo::rayIntersectsGround(atm.bottomRadius_m, r_m, mu); 
  
    // 視線方向の積分区間の決定
    // 上向きなら top atmosphere まで，下向きなら ground hit まで積分
    const float sMax = atmo::distanceToNearestAtmosphereBoundary(
        atm.bottomRadius_m,
        atm.topRadius_m,
        r_m,
        mu,
        viewHitsGround);

    if (sMax <= 0.0f) {
        outDeltaMultiple[outIdx] = 0.0f;
        return;
    }

    const float ds = sMax / static_cast<float>(integrationStepsView);

    float opticalDepthView = 0.0f;
    float accum = 0.0f;

    for (std::uint32_t i = 0; i < integrationStepsView; ++i) {
        const float sMid = (static_cast<float>(i) + 0.5f) * ds;
        const float rSample = atmo::pointRadiusAlongRay(r_m, mu, sMid);
        const float altitude_m = rSample - atm.bottomRadius_m;

        const float sigmaExtMid = extinctionAt(atm, altitude_m, iLambda);
        const float tauToMid = opticalDepthView + 0.5f * sigmaExtMid * ds;
        const float tView = ::expf(-tauToMid);

        const float muSample =
            atmo::clampCosine((r_m * mu + sMid) / rSample);

        const float muSSample =
            atmo::clampCosine((r_m * muS + sMid * nu) / rSample);

        const float density = sampleScatteringLUT(
            scatteringDensity,
            dims,
            atm,
            iLambda,
            rSample,
            muSample,
            muSSample,
            nu);

        accum += tView * density * ds;
        opticalDepthView += sigmaExtMid * ds;
    }

    outDeltaMultiple[outIdx] = accum;
}

__global__ void kernelProjectScatteringToObserverSky(
    DeviceAtmosphereParameters atm,
    DeviceLutDims dims,
    const float* deltaMultiple,
    float* outSkyAccum) {

    const int iNu    = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    const int iMu    = static_cast<int>(blockIdx.y * blockDim.y + threadIdx.y);
    const int packed = static_cast<int>(blockIdx.z * blockDim.z + threadIdx.z);

    const int iLambda = packed % dims.wavelengthCount;
    const int iMuS    = packed / dims.wavelengthCount;

    if (iNu >= dims.skyNu ||
        iMu >= dims.skyMu ||
        iMuS >= dims.skyMuS ||
        iLambda >= dims.wavelengthCount) {
        return;
    }

    const float uMu  = (static_cast<float>(iMu) + 0.5f) / static_cast<float>(dims.skyMu);
    const float uMuS = (static_cast<float>(iMuS) + 0.5f) / static_cast<float>(dims.skyMuS);
    const float uNu  = (static_cast<float>(iNu) + 0.5f) / static_cast<float>(dims.skyNu);

    const float mu  = -1.0f + 2.0f * atmo::clampf(uMu, 0.0f, 1.0f);
    const float muS = decodeMuSFromUnit(uMuS, atm.muSMin);
    const float nu  = decodeNuFromUnit(uNu);

    const float rObserver = atm.bottomRadius_m + atm.observerAltitude_m;

    const float value = sampleScatteringLUT(
        deltaMultiple,
        dims,
        atm,
        iLambda,
        rObserver,
        mu,
        muS,
        nu);

    const std::size_t outIdx = skyIndex(iMu, iMuS, iNu, iLambda, dims);
    outSkyAccum[outIdx] += value;
}

} // namespace

void buildDeltaSingleScattering(
    const BakeConfig& cfg,
    const std::vector<float>& transmittance,
    std::vector<float>& outDeltaRayleigh,
    std::vector<float>& outDeltaMie) {

    DeviceLutDims dims{};
    dims.transmittanceR  = cfg.lut.transmittanceR;
    dims.transmittanceMu = cfg.lut.transmittanceMu;
    dims.scatteringR     = cfg.lut.irradianceR; // 最初は irradianceR を流用
    dims.skyMu           = cfg.lut.nMu;
    dims.skyMuS          = cfg.lut.nMuS;
    dims.skyNu           = cfg.lut.nNu;
    dims.irradianceR     = cfg.lut.irradianceR;
    dims.irradianceMuS   = cfg.lut.irradianceMuS;
    dims.wavelengthCount = static_cast<int>(cfg.meta.wavelengths_nm.size());

    const std::size_t expectedTransmittanceSize =
        static_cast<std::size_t>(dims.transmittanceR) *
        static_cast<std::size_t>(dims.transmittanceMu) *
        static_cast<std::size_t>(dims.wavelengthCount);

    if (transmittance.size() != expectedTransmittanceSize) {
        throw std::runtime_error(
            "Transmittance LUT size mismatch: expected " +
            std::to_string(expectedTransmittanceSize) +
            ", got " + std::to_string(transmittance.size()));
    }

    const std::size_t total = scatteringElementCount(dims);
    const int blockSize = 256;
    const dim3 block(blockSize, 1, 1);
    const dim3 grid = make1DGrid(total, blockSize);

    outDeltaRayleigh.resize(total);
    outDeltaMie.resize(total);

    DeviceSpectralData dSpectra = uploadDeviceSpectralData(cfg.meta);

    DeviceAtmosphereParameters atm = makeDeviceAtmosphericParameters(cfg.meta);
    atm.wavelengths_nm        = nullptr;
    atm.solarIrradiance       = dSpectra.solarIrradiance;
    atm.rayleighScattering    = dSpectra.rayleighScattering;
    atm.mieScattering         = dSpectra.mieScattering;
    atm.mieExtinction         = dSpectra.mieExtinction;
    atm.absorptionExtinction  = dSpectra.absorptionExtinction;

    float* dTrans = nullptr;
    float* dDeltaRayleigh = nullptr;
    float* dDeltaMie = nullptr;

    checkCuda(cudaMalloc(&dTrans, transmittance.size() * sizeof(float)),
              "cudaMalloc transmittance LUT");
    checkCuda(cudaMalloc(&dDeltaRayleigh, total * sizeof(float)),
              "cudaMalloc delta rayleigh single");
    checkCuda(cudaMalloc(&dDeltaMie, total * sizeof(float)),
              "cudaMalloc delta mie single");

    try {
        checkCuda(cudaMemcpy(dTrans,
                             transmittance.data(),
                             transmittance.size() * sizeof(float),
                             cudaMemcpyHostToDevice),
                  "cudaMemcpy transmittance LUT");

        std::cout
            << "HOST dims: "
            << "transR="    << dims.transmittanceR
            << " transMu="  << dims.transmittanceMu
            << " scatR="    << dims.scatteringR
            << " skyMu="    << dims.skyMu
            << " skyMuS="   << dims.skyMuS
            << " skyNu="    << dims.skyNu
            << " irrR="     << dims.irradianceR
            << " irrMuS="   << dims.irradianceMuS
            << " lambda="   << dims.wavelengthCount
            << std::endl;

        kernelBuildDeltaSingleScattering<<<grid, block>>>(
            atm,
            dims,
            dTrans,
            static_cast<std::uint32_t>(cfg.lut.integrationStepsView),
            dDeltaRayleigh,
            dDeltaMie);

        checkCuda(cudaGetLastError(), "kernelBuildDeltaSingleScattering launch");
        checkCuda(cudaDeviceSynchronize(), "kernelBuildDeltaSingleScattering sync");

        checkCuda(cudaMemcpy(outDeltaRayleigh.data(),
                             dDeltaRayleigh,
                             total * sizeof(float),
                             cudaMemcpyDeviceToHost),
                  "cudaMemcpy delta rayleigh single");
        checkCuda(cudaMemcpy(outDeltaMie.data(),
                             dDeltaMie,
                             total * sizeof(float),
                             cudaMemcpyDeviceToHost),
                  "cudaMemcpy delta mie single");
    } catch (...) {
        cudaFree(dTrans);
        cudaFree(dDeltaRayleigh);
        cudaFree(dDeltaMie);
        freeDeviceSpectralData(dSpectra);
        throw;
    }

    cudaFree(dTrans);
    cudaFree(dDeltaRayleigh);
    cudaFree(dDeltaMie);
    freeDeviceSpectralData(dSpectra);
}


void buildMultipleScattering(
    const BakeConfig& cfg,
    const std::vector<float>& transmittance,
    const std::vector<float>& deltaRayleighSingle,
    const std::vector<float>& deltaMieSingle,
    std::vector<float>& outSkyMultiple) {

    DeviceLutDims dims{};
    dims.transmittanceR  = cfg.lut.transmittanceR;
    dims.transmittanceMu = cfg.lut.transmittanceMu;
    dims.scatteringR     = cfg.lut.irradianceR; // 最初は irradianceR を流用
    dims.skyMu           = cfg.lut.nMu;
    dims.skyMuS          = cfg.lut.nMuS;
    dims.skyNu           = cfg.lut.nNu;
    dims.irradianceR     = cfg.lut.irradianceR;
    dims.irradianceMuS   = cfg.lut.irradianceMuS;
    dims.wavelengthCount = static_cast<int>(cfg.meta.wavelengths_nm.size());

    const std::size_t transSize =
        static_cast<std::size_t>(dims.transmittanceR) *
        static_cast<std::size_t>(dims.transmittanceMu) *
        static_cast<std::size_t>(dims.wavelengthCount);

    const std::size_t scatteringSize =
        static_cast<std::size_t>(dims.scatteringR) *
        static_cast<std::size_t>(dims.skyMu) *
        static_cast<std::size_t>(dims.skyMuS) *
        static_cast<std::size_t>(dims.skyNu) *
        static_cast<std::size_t>(dims.wavelengthCount);

    const std::size_t skySize =
        static_cast<std::size_t>(dims.skyMu) *
        static_cast<std::size_t>(dims.skyMuS) *
        static_cast<std::size_t>(dims.skyNu) *
        static_cast<std::size_t>(dims.wavelengthCount);

    std::cout << "Checking buffer size..." << std::endl;
    if (transmittance.size() != transSize) {
        throw std::runtime_error("Transmittance LUT size mismatch");
    }
    if (deltaRayleighSingle.size() != scatteringSize) {
        throw std::runtime_error("DeltaRayleighSingle LUT size mismatch");
    }
    if (deltaMieSingle.size() != scatteringSize) {
        throw std::runtime_error("DeltaMieSingle LUT size mismatch");
    }

    outSkyMultiple.assign(skySize, 0.0f);

    if (!cfg.multiple.enabled || cfg.multiple.scatteringOrders <= 1) {
        return;
    }

    DeviceSpectralData dSpectra = uploadDeviceSpectralData(cfg.meta);

    DeviceAtmosphereParameters atm = makeDeviceAtmosphericParameters(cfg.meta);
    atm.wavelengths_nm        = nullptr;
    atm.solarIrradiance       = dSpectra.solarIrradiance;
    atm.rayleighScattering    = dSpectra.rayleighScattering;
    atm.mieScattering         = dSpectra.mieScattering;
    atm.mieExtinction         = dSpectra.mieExtinction;
    atm.absorptionExtinction  = dSpectra.absorptionExtinction;

    float* dTrans = nullptr;
    float* dDeltaRayleigh = nullptr;
    float* dDeltaMie = nullptr;
    float* dPrevScattering = nullptr;
    float* dScatteringDensity = nullptr;
    float* dDeltaMultiple = nullptr;
    float* dSkyAccum = nullptr;

    std::cout << "Allocating device memory..." << std::endl;
    checkCuda(cudaMalloc(&dTrans, transSize * sizeof(float)),
              "cudaMalloc transmittance LUT");
    checkCuda(cudaMalloc(&dDeltaRayleigh, scatteringSize * sizeof(float)),
              "cudaMalloc delta rayleigh single");
    checkCuda(cudaMalloc(&dDeltaMie, scatteringSize * sizeof(float)),
              "cudaMalloc delta mie single");
    checkCuda(cudaMalloc(&dPrevScattering, scatteringSize * sizeof(float)),
              "cudaMalloc prev scattering");
    checkCuda(cudaMalloc(&dScatteringDensity, scatteringSize * sizeof(float)),
              "cudaMalloc scattering density");
    checkCuda(cudaMalloc(&dDeltaMultiple, scatteringSize * sizeof(float)),
              "cudaMalloc delta multiple");
    checkCuda(cudaMalloc(&dSkyAccum, skySize * sizeof(float)),
              "cudaMalloc sky multiple accum");

    try {
        checkCuda(cudaMemcpy(dTrans,
                             transmittance.data(),
                             transSize * sizeof(float),
                             cudaMemcpyHostToDevice),
                  "cudaMemcpy transmittance LUT");
        checkCuda(cudaMemcpy(dDeltaRayleigh,
                             deltaRayleighSingle.data(),
                             scatteringSize * sizeof(float),
                             cudaMemcpyHostToDevice),
                  "cudaMemcpy delta rayleigh single");
        checkCuda(cudaMemcpy(dDeltaMie,
                             deltaMieSingle.data(),
                             scatteringSize * sizeof(float),
                             cudaMemcpyHostToDevice),
                  "cudaMemcpy delta mie single");
        checkCuda(cudaMemset(dSkyAccum, 0, skySize * sizeof(float)),
                  "cudaMemset sky multiple accum");

        const std::size_t total = scatteringElementCount(dims);
        const int blockSize = 256;
        const dim3 block(blockSize, 1, 1);
        const dim3 grid = make1DGrid(total, blockSize);

        std::cout << "Computing kernelAddSingleScattering..." << std::endl;
        kernelAddSingleScattering<<<grid, block>>>(
            atm,
            dims,
            dDeltaRayleigh,
            dDeltaMie,
            dPrevScattering);
        checkCuda(cudaGetLastError(), "kernelAddSingleScattering launch");
        checkCuda(cudaDeviceSynchronize(), "kernelAddSingleScattering sync");

        for (int order = 2; order <= cfg.multiple.scatteringOrders; ++order) {
            std::cout << "Start to compute " << order << " th scattering" << std::endl;
            std::cout << "Computing kernelBuildScatteringDensity..." << std::endl;
            kernelBuildScatteringDensity<<<grid, block>>>(
                atm,
                dims,
                dPrevScattering,
                dScatteringDensity);
            checkCuda(cudaGetLastError(), "kernelBuildScatteringDensity launch");
            checkCuda(cudaDeviceSynchronize(), "kernelBuildScatteringDensity sync");

            std::cout << "Computing kernelBuildDeltaMultipleScattering..." << std::endl;
            kernelBuildDeltaMultipleScattering<<<grid, block>>>(
                atm,
                dims,
                dTrans,
                dScatteringDensity,
                static_cast<std::uint32_t>(cfg.lut.integrationStepsView),
                dDeltaMultiple);
            checkCuda(cudaGetLastError(), "kernelBuildDeltaMultipleScattering launch");
            checkCuda(cudaDeviceSynchronize(), "kernelBuildDeltaMultipleScattering sync");

            dim3 blockSky(8, 8, 2);
            dim3 gridSky(
                static_cast<unsigned int>((dims.skyNu + blockSky.x - 1) / blockSky.x),
                static_cast<unsigned int>((dims.skyMu + blockSky.y - 1) / blockSky.y),
                static_cast<unsigned int>(((dims.skyMuS * dims.wavelengthCount) +
                                           blockSky.z - 1) / blockSky.z));

            std::cout << "Computing kernelProjectScatteringToObserverSky..." << std::endl;
            kernelProjectScatteringToObserverSky<<<gridSky, blockSky>>>(
                atm,
                dims,
                dDeltaMultiple,
                dSkyAccum);
            checkCuda(cudaGetLastError(), "kernelProjectScatteringToObserverSky launch");
            checkCuda(cudaDeviceSynchronize(), "kernelProjectScatteringToObserverSky sync");

            // 次の次数では今回の deltaMultiple を prevScattering として使う
            float* tmp = dPrevScattering;
            dPrevScattering = dDeltaMultiple;
            dDeltaMultiple = tmp;
        }

        checkCuda(cudaMemcpy(outSkyMultiple.data(),
                             dSkyAccum,
                             skySize * sizeof(float),
                             cudaMemcpyDeviceToHost),
                  "cudaMemcpy sky multiple accum");
    } catch (...) {
        cudaFree(dTrans);
        cudaFree(dDeltaRayleigh);
        cudaFree(dDeltaMie);
        cudaFree(dPrevScattering);
        cudaFree(dScatteringDensity);
        cudaFree(dDeltaMultiple);
        cudaFree(dSkyAccum);
        freeDeviceSpectralData(dSpectra);
        throw;
    }

    cudaFree(dTrans);
    cudaFree(dDeltaRayleigh);
    cudaFree(dDeltaMie);
    cudaFree(dPrevScattering);
    cudaFree(dScatteringDensity);
    cudaFree(dDeltaMultiple);
    cudaFree(dSkyAccum);
    freeDeviceSpectralData(dSpectra);
}

} // namespace bake