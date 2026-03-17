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

struct DeviceSpectralData {
  float* solarIrradiance        = nullptr;
  float* rayleighScattering     = nullptr;
  float* mieScattering          = nullptr;
  float* mieExtinction          = nullptr;
  float* absorptionExtinction   = nullptr;
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

// 今回の最終 sky LUT は「地上から見える空」専用として [0,1] を採用
ATMO_HD inline float decodeSkyMuFromUnit(float u) {
    return -1.0f + 2.0f * atmo::clampf(u, 0.0f, 1.0f);
}

ATMO_HD inline float encodeSkyMuToUnit(float mu) {
    return 0.5f * (atmo::clampCosine(mu) + 1.0f);
}

ATMO_HD inline float decodeMuSFromUnit(float u, float muSMin) {
    return muSMin + atmo::clampf(u, 0.0f, 1.0f) * (1.0f - muSMin);
}

ATMO_HD inline float decodeNuFromUnit(float u) {
    return -1.0f + 2.0f * atmo::clampf(u, 0.0f, 1.0f);
}

__global__ void kernelBuildSkySingleScatter(
  DeviceAtmosphereParameters atm,
  DeviceLutDims dims,
  const float* transmittanceLUT,
  std::uint32_t integrationStepsView,
  float* outRayleigh,
  float* outMie
) {  

  // スレッドが担当するインデックスの決定
  const int iNu     = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
  const int iMu     = static_cast<int>(blockIdx.y * blockDim.y + threadIdx.y);
  const int packed  = static_cast<int>(blockIdx.z * blockDim.z + threadIdx.z);  
  const int iMuS    = packed % dims.skyMuS;
  const int iLambda = packed / dims.skyMuS; 
  
  // インデックスの範囲外のスレッドは終了させる
  if (iNu >= dims.skyNu ||
      iMu >= dims.skyMu ||
      iMuS >= dims.skyMuS ||
      iLambda >= dims.wavelengthCount) {
      return;
  }
  
  // インデックスを物理的な変数に変換
  const float uMu  = (static_cast<float>(iMu) + 0.5f) / static_cast<float>(dims.skyMu);
  const float uMuS = (static_cast<float>(iMuS) + 0.5f) / static_cast<float>(dims.skyMuS);
  const float uNu  = (static_cast<float>(iNu) + 0.5f) / static_cast<float>(dims.skyNu); 
  
  const float mu  = decodeSkyMuFromUnit(uMu);             // 視線方向と天頂方向がなす角の余弦 [-1, 1]
  const float muS = decodeMuSFromUnit(uMuS, atm.muSMin);  // 太陽方向と天頂方向がなす角の余弦 [munin, 1]
  const float nu  = decodeNuFromUnit(uNu);                // 視線方向と太陽方向がなす角の余弦 [-1, 1] 
  
  const float rObserver = atm.bottomRadius_m + atm.observerAltitude_m;  // 観測者がいる位置
  const std::size_t outIdx = skyIndex(iMu, iMuS, iNu, iLambda, dims); 
  
  // 視線が地面と交差するかどうかの判定
  const bool viewHitsGround = atmo::rayIntersectsGround(atm.bottomRadius_m, rObserver, mu); 
  
  // 視線方向の積分区間の決定
  // 上向きなら top atmosphere まで，下向きなら ground hit まで積分
  const float sMax = atmo::distanceToNearestAtmosphereBoundary(
    atm.bottomRadius_m,
    atm.topRadius_m,
    rObserver,
    mu,
    viewHitsGround);

  if (sMax <= 0.0f) {
    outRayleigh[outIdx] = 0.0f;
    outMie[outIdx] = 0.0f;
    return;
  }
  
  const float ds = sMax / static_cast<float>(integrationStepsView);
  
  // 視線上をレイマーチングし，各サンプリング点で寄与を計算
  float opticalDepthView = 0.0f;
  float accumRayleigh = 0.0f;
  float accumMie = 0.0f;  
  for (std::uint32_t i = 0; i < integrationStepsView; ++i) {

      // サンプリングする位置の決定
      const float sMid = (static_cast<float>(i) + 0.5f) * ds;               // 区間の中点
      const float rSample = atmo::pointRadiusAlongRay(rObserver, mu, sMid); // 中点に対応する散乱点の半径
      const float altitude_m = rSample - atm.bottomRadius_m;                // 中点に対応する散乱点の高度
      
      // そのサンプル点における消散係数 (Rayleigh, Mie もすべて)
      const float sigmaExtMid = extinctionAt(atm, altitude_m, iLambda);
      
      // 観測者からサンプル点までの透過率の計算
      const float tauToMid = opticalDepthView + 0.5f * sigmaExtMid * ds;
      const float tView = ::expf(-tauToMid);  // 透過率
      
      // 散乱点から見た太陽方向と天頂方向とのなす角の余弦
      const float muSSample =
          atmo::clampCosine((rObserver * muS + sMid * nu) / rSample); 
      
      // 散乱点から太陽までの透過率
      const float tSun = sampleTransmittanceLUT(
          transmittanceLUT,
          dims,
          atm,
          iLambda,
          rSample,
          muSSample); 
      
      // 散乱点での Rayleigh, Mie 散乱を起こす媒質の数密度
      const float rhoR = numberDensityRayleigh(atm, altitude_m);
      const float rhoM = numberDensityMie(atm, altitude_m); 
      
      // 各成分の透過率の計算
      accumRayleigh += tView * tSun * rhoR * ds;
      accumMie      += tView * tSun * rhoM * ds;  
      
      opticalDepthView += sigmaExtMid * ds;
  }
  
  // 位相関数以外の
  outRayleigh[outIdx] =
      accumRayleigh *
      atm.solarIrradiance[iLambda] *
      atm.rayleighScattering[iLambda];  
  outMie[outIdx] =
      accumMie *
      atm.solarIrradiance[iLambda] *
      atm.mieScattering[iLambda];
}

}  // namespace

void buildSkySingleScatter(
    const BakeConfig& cfg,
    const std::vector<float>& transmittance,
    std::vector<float>& outRayleigh,
    std::vector<float>& outMie) {

    // LUT のサイズの決定
    DeviceLutDims dims{};
    dims.transmittanceR  = cfg.lut.transmittanceR;
    dims.transmittanceMu = cfg.lut.transmittanceMu;
    dims.skyMu           = cfg.lut.nMu;
    dims.skyMuS          = cfg.lut.nMuS;
    dims.skyNu           = cfg.lut.nNu;
    dims.wavelengthCount = static_cast<int>(cfg.meta.wavelengths_nm.size());

    // 出力する配列のサイズの確保
    const std::size_t total =
        static_cast<std::size_t>(dims.skyMu) *
        static_cast<std::size_t>(dims.skyMuS) *
        static_cast<std::size_t>(dims.skyNu) *
        static_cast<std::size_t>(dims.wavelengthCount);

    outRayleigh.resize(total);
    outMie.resize(total);

    // metadata 内にある配列の GPU 側への転送
    DeviceSpectralData dSpectra = uploadDeviceSpectralData(cfg.meta);

    // metadata の GPU 計算向け構造体への変換
    DeviceAtmosphereParameters atm = makeDeviceAtmosphericParameters(cfg.meta);
    atm.wavelengths_nm        = nullptr;
    atm.solarIrradiance       = dSpectra.solarIrradiance;
    atm.rayleighScattering    = dSpectra.rayleighScattering;
    atm.mieScattering         = dSpectra.mieScattering;
    atm.mieExtinction         = dSpectra.mieExtinction;
    atm.absorptionExtinction  = dSpectra.absorptionExtinction;

    // 各種 LUT のデバイスメモリの確保
    float* dTrans = nullptr;      // Transmittance LUT (中間表現)
    float* dRayleigh = nullptr;   // Rayleigh 単一散乱の LUT
    float* dMie = nullptr;        // Mie 単一散乱の LUT

    checkCuda(cudaMalloc(&dTrans, transmittance.size() * sizeof(float)),
              "cudaMalloc transmittance");
    checkCuda(cudaMalloc(&dRayleigh, total * sizeof(float)),
              "cudaMalloc sky rayleigh");
    checkCuda(cudaMalloc(&dMie, total * sizeof(float)),
              "cudaMalloc sky mie");

    try {
        // transmittance データの転送
        checkCuda(cudaMemcpy(dTrans,
                             transmittance.data(),
                             transmittance.size() * sizeof(float),
                             cudaMemcpyHostToDevice),
                  "cudaMemcpy transmittance");

        // カーネルの起動
        dim3 block(8, 8, 2);
        dim3 grid(
            static_cast<unsigned int>((dims.skyNu + block.x - 1) / block.x),
            static_cast<unsigned int>((dims.skyMu + block.y - 1) / block.y),
            static_cast<unsigned int>(((dims.wavelengthCount * dims.skyMuS) + block.z - 1) / block.z));

        kernelBuildSkySingleScatter<<<grid, block>>>(
            atm,
            dims,
            dTrans,
            static_cast<std::uint32_t>(cfg.lut.integrationStepsView),
            dRayleigh,
            dMie);

        checkCuda(cudaGetLastError(), "kernelBuildSkySingleScatter launch");
        checkCuda(cudaDeviceSynchronize(), "kernelBuildSkySingleScatter sync");

        // 結果を host 側へ戻す
        checkCuda(cudaMemcpy(outRayleigh.data(),
                             dRayleigh,
                             total * sizeof(float),
                             cudaMemcpyDeviceToHost),
                  "cudaMemcpy sky rayleigh");
        checkCuda(cudaMemcpy(outMie.data(),
                             dMie,
                             total * sizeof(float),
                             cudaMemcpyDeviceToHost),
                  "cudaMemcpy sky mie");
    } catch (...) {
        cudaFree(dTrans);
        cudaFree(dRayleigh);
        cudaFree(dMie);
        freeDeviceSpectralData(dSpectra);
        throw;
    }

    // 不要になったメモリの開放など
    cudaFree(dTrans);
    cudaFree(dRayleigh);
    cudaFree(dMie);
    freeDeviceSpectralData(dSpectra);
}

}  // namespace bake
