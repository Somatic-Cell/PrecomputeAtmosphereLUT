#include "config_io.h"
#include "device_common.cuh"

#include <cuda_runtime.h>

#include <cmath>
#include <stdexcept>
#include <vector>

namespace bake {
namespace {

void checkCuda(cudaError_t err, const char* what) {
  if (err != cudaSuccess) {
    throw std::runtime_error(std::string(what) + ": " + cudaGetErrorString(err));
  }
}

struct DeviceOpticalCoefficients {
  float* rayleighScattering   = nullptr;
  float* mieScattering        = nullptr;
  float* mieExtinction        = nullptr;
  float* absorptionExtinction = nullptr;
};

void freeDeviceOpticalCoefficients(DeviceOpticalCoefficients& d) {
  cudaFree(d.rayleighScattering);
  cudaFree(d.mieScattering);
  cudaFree(d.mieExtinction);
  cudaFree(d.absorptionExtinction);
  d = {};
}

DeviceOpticalCoefficients uploadDeviceOpticalCoefficients(
  const atmo::AtmosphereMetadata& meta
){
  DeviceOpticalCoefficients d{};
  const std::size_t bytes = meta.wavelengths_nm.size() * sizeof(float);

  auto upload = [bytes](
    float*& dst,
    const std::vector<float>& src,
    const char* name
  ) {
    checkCuda(cudaMalloc(&dst, bytes), name);
    checkCuda(cudaMemcpy(dst, src.data(), bytes, cudaMemcpyHostToDevice), name);
  };

  try {
    upload(d.rayleighScattering,    meta.rayleighScattering,    "cudaMemcpy rayleighScattering");
    upload(d.mieScattering,         meta.mieScattering,         "cudaMemcpy mieScattering");
    upload(d.mieExtinction,         meta.mieExtinction,         "cudaMemcpy mieExtinction");
    upload(d.absorptionExtinction,  meta.absorptionExtinction,  "cudaMemcpy absorptionExtinction");
  } catch (...) {
    freeDeviceOpticalCoefficients(d);
    throw;
  }

  return d;
}

// Bruneton 2017 寄りの transmittance LUT parameterization の逆写像
// device_common.cuh 側には forward (r, mu) -> (u, v) しかないため
// bake 側では LUT 格子中心 (u, v) から (r, mu) を復元して積分する
ATMO_HD inline void transmittanceRMuFromUV(
    const DeviceAtmosphereParameters& atm,
    float u,
    float v,
    float& r_m,
    float& mu) {

    const float bottom = atm.bottomRadius_m;
    const float top    = atm.topRadius_m;

    const float H   = atmo::safeSqrt(top * top - bottom * bottom);
    const float rho = H * atmo::clampf(v, 0.0f, 1.0f);
    r_m = atmo::safeSqrt(rho * rho + bottom * bottom);

    const float dMin = top - r_m;
    const float dMax = rho + H;
    const float d    = dMin + atmo::clampf(u, 0.0f, 1.0f) * (dMax - dMin);

    if (d <= 0.0f) {
        mu = 1.0f;
        return;
    }

    // top^2 = r^2 + d^2 + 2 r d mu
    mu = (top * top - r_m * r_m - d * d) / (2.0f * r_m * d);
    mu = atmo::clampCosine(mu);
}

__global__ void kernelBuildTransmittanceLUT(
    DeviceAtmosphereParameters atm,
    DeviceLutDims dims,
    std::uint32_t integrationSteps,
    float* outLut) {

    const int iMu     = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    const int iR      = static_cast<int>(blockIdx.y * blockDim.y + threadIdx.y);
    const int iLambda = static_cast<int>(blockIdx.z * blockDim.z + threadIdx.z);

    if (iMu >= dims.transmittanceMu ||
        iR >= dims.transmittanceR ||
        iLambda >= dims.wavelengthCount) {
        return;
    }

    const float u = (static_cast<float>(iMu) + 0.5f) / static_cast<float>(dims.transmittanceMu);
    const float v = (static_cast<float>(iR) + 0.5f) / static_cast<float>(dims.transmittanceR);

    float r_m  = 0.0f;
    float mu   = 1.0f;
    transmittanceRMuFromUV(atm, u, v, r_m, mu);

    const std::size_t outIdx = transmittanceIndex(iR, iMu, iLambda, dims);

    if (atmo::rayIntersectsGround(atm.bottomRadius_m, r_m, mu)) {
        outLut[outIdx] = 0.0f;
        return;
    }

    const float sMax = atmo::distanceToTopAtmosphereBoundary(atm.topRadius_m, r_m, mu);
    const float ds   = sMax / static_cast<float>(integrationSteps);

    float opticalDepth = 0.0f;
    for (std::uint32_t i = 0; i < integrationSteps; ++i) {
        const float sMid = (static_cast<float>(i) + 0.5f) * ds;
        const float rSample = atmo::pointRadiusAlongRay(r_m, mu, sMid);
        const float altitude_m = rSample - atm.bottomRadius_m;

        opticalDepth += extinctionAt(atm, altitude_m, iLambda) * ds;
    }

    outLut[outIdx] = ::expf(-opticalDepth);
}

} // namespace

void buildTransmittanceLUT(
    const BakeConfig& cfg,
    std::vector<float>& outTransmittance) {

    DeviceLutDims dims{};
    dims.transmittanceR  = cfg.lut.transmittanceR;
    dims.transmittanceMu = cfg.lut.transmittanceMu;
    dims.wavelengthCount = static_cast<int>(cfg.meta.wavelengths_nm.size());

    const std::size_t total =
        static_cast<std::size_t>(dims.transmittanceR) *
        static_cast<std::size_t>(dims.transmittanceMu) *
        static_cast<std::size_t>(dims.wavelengthCount);

    outTransmittance.resize(total);

    DeviceOpticalCoefficients dCoeff = uploadDeviceOpticalCoefficients(cfg.meta);

    DeviceAtmosphereParameters atm = makeDeviceAtmosphericParameters(cfg.meta);

    // device 上で参照されるポインタだけ差し替える
    atm.wavelengths_nm        = nullptr; // この kernel では未使用
    atm.solarIrradiance       = nullptr; // この kernel では未使用
    atm.rayleighScattering    = dCoeff.rayleighScattering;
    atm.mieScattering         = dCoeff.mieScattering;
    atm.mieExtinction         = dCoeff.mieExtinction;
    atm.absorptionExtinction  = dCoeff.absorptionExtinction;

    float* dOut = nullptr;
    checkCuda(cudaMalloc(&dOut, total * sizeof(float)),
              "cudaMalloc transmittance LUT");

    try {
        dim3 block(8, 8, 2);
        dim3 grid(
            static_cast<unsigned int>((dims.transmittanceMu + block.x - 1) / block.x),
            static_cast<unsigned int>((dims.transmittanceR  + block.y - 1) / block.y),
            static_cast<unsigned int>((dims.wavelengthCount + block.z - 1) / block.z));

        kernelBuildTransmittanceLUT<<<grid, block>>>(
            atm,
            dims,
            static_cast<std::uint32_t>(cfg.lut.integrationStepsTransmittance),
            dOut);

        checkCuda(cudaGetLastError(), "kernelBuildTransmittanceLUT launch");
        checkCuda(cudaDeviceSynchronize(), "kernelBuildTransmittanceLUT sync");

        checkCuda(cudaMemcpy(outTransmittance.data(),
                             dOut,
                             total * sizeof(float),
                             cudaMemcpyDeviceToHost),
                  "cudaMemcpy transmittance LUT");
    } catch (...) {
        cudaFree(dOut);
        freeDeviceOpticalCoefficients(dCoeff);
        throw;
    }

    cudaFree(dOut);
    freeDeviceOpticalCoefficients(dCoeff);
}

}  // namespace bake
