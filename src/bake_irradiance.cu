#include "bake_common.cuh"
#include "device_common.cuh"

#include <cuda_runtime.h>

#include <cstddef>
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

struct DeviceSolarSpectrum {
    float* solarIrradiance = nullptr;
};

void freeDeviceSolarSpectrum(DeviceSolarSpectrum& d) {
    cudaFree(d.solarIrradiance);
    d.solarIrradiance = nullptr;
}

DeviceSolarSpectrum uploadDeviceSolarSpectrum(const atmo::AtmosphereMetadata& meta) {
    DeviceSolarSpectrum d{};
    const std::size_t bytes = meta.wavelengths_nm.size() * sizeof(float);

    try {
        checkCuda(cudaMalloc(&d.solarIrradiance, bytes),
                  "cudaMalloc solarIrradiance");
        checkCuda(cudaMemcpy(d.solarIrradiance,
                             meta.solarIrradiance.data(),
                             bytes,
                             cudaMemcpyHostToDevice),
                  "cudaMemcpy solarIrradiance");
    } catch (...) {
        freeDeviceSolarSpectrum(d);
        throw;
    }

    return d;
}

ATMO_HD inline float decodeMuSFromUnit(float u, float muSMin) {
    return muSMin + atmo::clampf(u, 0.0f, 1.0f) * (1.0f - muSMin);
}

// irradiance LUT の r 軸は transmittance と同じ rho/H parameterization を使う
ATMO_HD inline float radiusFromIrradianceV(
    const DeviceAtmosphereParameters& atm,
    float v) {

    const float bottom = atm.bottomRadius_m;
    const float top    = atm.topRadius_m;

    const float H   = atmo::safeSqrt(top * top - bottom * bottom);
    const float rho = H * atmo::clampf(v, 0.0f, 1.0f);

    return atmo::safeSqrt(rho * rho + bottom * bottom);
}

__global__ void kernelBuildDirectIrradiance(
    DeviceAtmosphereParameters atm,
    DeviceLutDims dims,
    const float* transmittanceLUT,
    float* outDirectIrradiance) {

    const int iMuS    = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    const int iR      = static_cast<int>(blockIdx.y * blockDim.y + threadIdx.y);
    const int iLambda = static_cast<int>(blockIdx.z * blockDim.z + threadIdx.z);

    if (iMuS >= dims.irradianceMuS ||
        iR >= dims.irradianceR ||
        iLambda >= dims.wavelengthCount) {
        return;
    }

    const float uMuS =
        (static_cast<float>(iMuS) + 0.5f) / static_cast<float>(dims.irradianceMuS);
    const float vR =
        (static_cast<float>(iR) + 0.5f) / static_cast<float>(dims.irradianceR);

    const float muS = decodeMuSFromUnit(uMuS, atm.muSMin);
    const float r_m = radiusFromIrradianceV(atm, vR);

    const std::size_t outIdx = irradianceIndex(iR, iMuS, iLambda, dims);

    // 太陽方向が地面と交差するなら direct irradiance は 0
    if (atmo::rayIntersectsGround(atm.bottomRadius_m, r_m, muS)) {
        outDirectIrradiance[outIdx] = 0.0f;
        return;
    }

    const float T = sampleTransmittanceLUT(
        transmittanceLUT,
        dims,
        atm,
        iLambda,
        r_m,
        muS);

    // まずは位相関数の影響なし, 太陽有限角補正なし
    outDirectIrradiance[outIdx] = atm.solarIrradiance[iLambda] * T;
}

} // namespace

void buildDirectIrradiance(
    const BakeConfig& cfg,
    const std::vector<float>& transmittance,
    std::vector<float>& outDirectIrradiance) {

    DeviceLutDims dims{};
    dims.transmittanceR  = cfg.lut.transmittanceR;
    dims.transmittanceMu = cfg.lut.transmittanceMu;
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

    const std::size_t total =
        static_cast<std::size_t>(dims.irradianceR) *
        static_cast<std::size_t>(dims.irradianceMuS) *
        static_cast<std::size_t>(dims.wavelengthCount);

    outDirectIrradiance.resize(total);

    DeviceSolarSpectrum dSolar = uploadDeviceSolarSpectrum(cfg.meta);

    DeviceAtmosphereParameters atm = makeDeviceAtmosphericParameters(cfg.meta);
    atm.wavelengths_nm        = nullptr;
    atm.solarIrradiance       = dSolar.solarIrradiance;
    atm.rayleighScattering    = nullptr;
    atm.mieScattering         = nullptr;
    atm.mieExtinction         = nullptr;
    atm.absorptionExtinction  = nullptr;

    float* dTrans = nullptr;
    float* dOut   = nullptr;

    checkCuda(cudaMalloc(&dTrans, transmittance.size() * sizeof(float)),
              "cudaMalloc transmittance LUT");
    checkCuda(cudaMalloc(&dOut, total * sizeof(float)),
              "cudaMalloc direct irradiance LUT");

    try {
        checkCuda(cudaMemcpy(dTrans,
                             transmittance.data(),
                             transmittance.size() * sizeof(float),
                             cudaMemcpyHostToDevice),
                  "cudaMemcpy transmittance LUT");

        dim3 block(8, 8, 2);
        dim3 grid(
            static_cast<unsigned int>((dims.irradianceMuS + block.x - 1) / block.x),
            static_cast<unsigned int>((dims.irradianceR   + block.y - 1) / block.y),
            static_cast<unsigned int>((dims.wavelengthCount + block.z - 1) / block.z));

        kernelBuildDirectIrradiance<<<grid, block>>>(
            atm,
            dims,
            dTrans,
            dOut);

        checkCuda(cudaGetLastError(), "kernelBuildDirectIrradiance launch");
        checkCuda(cudaDeviceSynchronize(), "kernelBuildDirectIrradiance sync");

        checkCuda(cudaMemcpy(outDirectIrradiance.data(),
                             dOut,
                             total * sizeof(float),
                             cudaMemcpyDeviceToHost),
                  "cudaMemcpy direct irradiance LUT");
    } catch (...) {
        cudaFree(dTrans);
        cudaFree(dOut);
        freeDeviceSolarSpectrum(dSolar);
        throw;
    }

    cudaFree(dTrans);
    cudaFree(dOut);
    freeDeviceSolarSpectrum(dSolar);
}

} // namespace bake