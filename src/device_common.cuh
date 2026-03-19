#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>


#include "atmosphere_math.h"

namespace bake {

// -----------------------------------------------------------------------------
// GPU に持ち込むパラメータ
// -----------------------------------------------------------------------------

struct DeviceAtmosphereParameters {
  float bottomRadius_m      = 0.0f;
  float topRadius_m         = 0.0f;
  float observerAltitude_m  = 0.0f;
  float muSMin              = -0.2f;

  atmo::DensityProfile rayleighDensity{};
  atmo::DensityProfile mieDensity{};
  atmo::DensityProfile absorptionDensity{};

  const float* wavelengths_nm       = nullptr;
  const float* solarIrradiance      = nullptr;
  const float* rayleighScattering   = nullptr;
  const float* mieScattering        = nullptr;
  const float* mieExtinction        = nullptr;
  const float* absorptionExtinction = nullptr;

  int wavelengthCount = 0;

  float miePhaseFunctionG = 0.8f;
  std::uint32_t scatteringOrders = 7u;
};

// LUT のチャネル数だけをまとめる
struct DeviceLutDims {
  int transmittanceR = 0;
  int transmittanceMu = 0;

  int skyMu   = 0;
  int skyMuS  = 0;
  int skyNu   = 0;

  int irradianceR = 0;
  int irradianceMuS = 0;

  int wavelengthCount = 0;

  int scatteringR     = 0;
};


// -----------------------------------------------------------------------------
// Host 側: metadata -> device params 変換
// -----------------------------------------------------------------------------

inline DeviceAtmosphereParameters makeDeviceAtmosphericParameters(
  const atmo::AtmosphereMetadata& m    
) {
  const std::size_t n = m.wavelengths_nm.size();
  auto requireSameSize = [n](const std::vector<float>& v, const char* name) {
    if(v.size() != n) {
      throw std::runtime_error(
        std::string("Spectrum length mismatch for ") + name +
        ": expected " + std::to_string(n) +
        ", got " + std::to_string(v.size()) 
      );
    }
  };

  requireSameSize(m.solarIrradiance,      "solarIrradiance");
  requireSameSize(m.rayleighScattering,   "rayleighScattering");
  requireSameSize(m.mieScattering,        "mieScattering");
  requireSameSize(m.mieExtinction,        "mieExtinction");
  requireSameSize(m.absorptionExtinction, "absorptionExtinction");

  DeviceAtmosphereParameters out{};

  out.bottomRadius_m      = m.bottomRadius_m;
  out.topRadius_m         = m.topRadius_m;
  out.observerAltitude_m  = m.observerAltitude_m;
  out.muSMin              = m.muSMin;

  out.rayleighDensity     = m.rayleighDensity;
  out.mieDensity          = m.mieDensity;
  out.absorptionDensity   = m.absorptionDensity;

  out.wavelengths_nm      = m.wavelengths_nm.empty()      ? nullptr : m.wavelengths_nm.data();
  out.solarIrradiance     = m.solarIrradiance.empty()     ? nullptr : m.solarIrradiance.data();
  out.rayleighScattering  = m.rayleighScattering.empty()  ? nullptr : m.rayleighScattering.data();
  out.mieScattering       = m.mieScattering.empty()       ? nullptr : m.mieScattering.data();
  out.mieExtinction       = m.mieExtinction.empty()       ? nullptr : m.mieExtinction.data();
  out.absorptionExtinction = m.absorptionExtinction.empty() ? nullptr : m.absorptionExtinction.data();

  out.wavelengthCount   = static_cast<int>(n);
  out.miePhaseFunctionG = m.miePhaseFunctionG;
  out.scatteringOrders  = m.scatteringOrders;
  
  return out;
}

// -----------------------------------------------------------------------------
// Host/Device 共通ユーティリティ
// -----------------------------------------------------------------------------

ATMO_HD inline int clampIndex(int x, int n) {
  return x < 0 ? 0 : (x >= n ? n - 1 : x);
}

ATMO_HD inline float lerpf(float a, float b, float t) {
  return a + (b - a) * t;
}

// 中間表現用 LUT のサンプリング位置の決定 
// 3D 情報が 1D 空間上に保持
ATMO_HD inline std::size_t transmittanceIndex(
  int iR,
  int iMu,
  int iLambda,
  const DeviceLutDims& dims
) {
  return static_cast<std::size_t>(
    ((iR * dims.transmittanceMu) + iMu) * dims.wavelengthCount + iLambda 
  );
}

// 3D 情報が 1D 空間上に保持
ATMO_HD inline std::size_t irradianceIndex(
  int iR,
  int iMuS,
  int iLambda,
  const DeviceLutDims& dims
){
  return static_cast<std::size_t>(
    ((iR * dims.irradianceMuS) + iMuS) * dims.wavelengthCount + iLambda
  );
}

// 4D 情報が 1D 空間上に保持
ATMO_HD inline std::size_t skyIndex(
  int iNu,
  int iMu,
  int iMuS,
  int iLambda,
  const DeviceLutDims& dims
){
  return static_cast<std::size_t>(
    ((((iNu * dims.skyMu) + iMu) * dims.skyMuS) + iMuS) *
    dims.wavelengthCount + iLambda
  );
}

ATMO_HD inline std::size_t scatteringIndex(
    int iR,
    int iMu,
    int iMuS,
    int iNu,
    int iLambda,
    const DeviceLutDims& dims) {

      return ((((
        static_cast<std::size_t>(iR)
        * static_cast<std::size_t>(dims.skyMu)
        + static_cast<std::size_t>(iMu))
        * static_cast<std::size_t>(dims.skyMuS)
        + static_cast<std::size_t>(iMuS))
        * static_cast<std::size_t>(dims.skyNu)
        + static_cast<std::size_t>(iNu))
        * static_cast<std::size_t>(dims.wavelengthCount)
        + static_cast<std::size_t>(iLambda));

}
// -----------------------------------------------------------------------------
// 物理係数の評価
// -----------------------------------------------------------------------------

// ある高度の数密度を返す
ATMO_HD inline float numberDensityRayleigh(
    const DeviceAtmosphereParameters& atm,
    float altitude_m 
) {
    return atmo::evalProfileDensity(atm.rayleighDensity, altitude_m);
}

ATMO_HD inline float numberDensityMie(
    const DeviceAtmosphereParameters& atm,
    float altitude_m  
) {
    return atmo::evalProfileDensity(atm.mieDensity, altitude_m);
}

ATMO_HD inline float numberDensityAbsorption(
    const DeviceAtmosphereParameters& atm,
    float altitude_m  
) { 
    return atmo::evalProfileDensity(atm.absorptionDensity, altitude_m);
}

// ある高度の散乱係数を返す
ATMO_HD inline float rayleighScatteringAt(
    const DeviceAtmosphereParameters& atm,
    float altitude_m, // 高度
    int iLambda       // 波長のインデックス
) {
    return numberDensityRayleigh(atm, altitude_m) * atm.rayleighScattering[iLambda];
}

ATMO_HD inline float mieScatteringAt(
    const DeviceAtmosphereParameters& atm,
    float altitude_m, // 高度
    int iLambda       // 波長のインデックス
) {
    return numberDensityMie(atm, altitude_m) * atm.mieScattering[iLambda];
}

// ある高度の消散係数を返す
ATMO_HD inline float mieExtinctionAt(
    const DeviceAtmosphereParameters& atm,
    float altitude_m, // 高度
    int iLambda       // 波長のインデックス
) {
    return numberDensityMie(atm, altitude_m) * atm.mieExtinction[iLambda];
}

// ある高度の吸収係数を返す
ATMO_HD inline float absorptionExtinctionAt(
    const DeviceAtmosphereParameters& atm,
    float altitude_m,
    int iLambda) {
    return numberDensityAbsorption(atm, altitude_m) * atm.absorptionExtinction[iLambda];
}

// ある高度の消散係数を返す
ATMO_HD inline float extinctionAt(
    const DeviceAtmosphereParameters& atm,
    float altitude_m,
    int iLambda) {
    return rayleighScatteringAt(atm, altitude_m, iLambda) +
           mieExtinctionAt(atm, altitude_m, iLambda) +
           absorptionExtinctionAt(atm, altitude_m, iLambda);
}

// -----------------------------------------------------------------------------
// Transmittance LUT parameterization
// Bruneton 2017 寄りの (r, mu) -> (u, v) 近似マッピング
// -----------------------------------------------------------------------------


ATMO_HD inline void transmittanceUVFromRMu(
    const DeviceAtmosphereParameters& atm,
    float r_m,
    float mu,
    float& u,
    float& v) {

    const float bottom = atm.bottomRadius_m;
    const float top = atm.topRadius_m;

    const float H = atmo::safeSqrt(top * top - bottom * bottom);
    const float rho = atmo::safeSqrt(r_m * r_m - bottom * bottom);

    const float d = atmo::distanceToTopAtmosphereBoundary(top, r_m, mu);
    const float dMin = top - r_m;
    const float dMax = rho + H;

    const float xMu = (dMax > dMin) ? ((d - dMin) / (dMax - dMin)) : 0.0f;
    const float xR = (H > 0.0f) ? (rho / H) : 0.0f;

    u = atmo::clampf(xMu, 0.0f, 1.0f);
    v = atmo::clampf(xR, 0.0f, 1.0f);
}

// 1 波長ぶんだけ 2D bilinear sample する
// table は transmittanceIndex() 規約で並んでいることが前提
ATMO_HD inline float sampleTransmittanceLUT(
    const float* table,
    const DeviceLutDims& dims,
    const DeviceAtmosphereParameters& atm,
    int iLambda,
    float r_m,
    float mu) {

    float u = 0.0f;
    float v = 0.0f;
    transmittanceUVFromRMu(atm, r_m, mu, u, v);

    const float x = u * float(dims.transmittanceMu - 1);
    const float y = v * float(dims.transmittanceR - 1);

    const int x0 = clampIndex(static_cast<int>(x), dims.transmittanceMu);
    const int y0 = clampIndex(static_cast<int>(y), dims.transmittanceR);
    const int x1 = clampIndex(x0 + 1, dims.transmittanceMu);
    const int y1 = clampIndex(y0 + 1, dims.transmittanceR);

    const float tx = x - float(x0);
    const float ty = y - float(y0);

    const float f00 = table[transmittanceIndex(y0, x0, iLambda, dims)];
    const float f10 = table[transmittanceIndex(y0, x1, iLambda, dims)];
    const float f01 = table[transmittanceIndex(y1, x0, iLambda, dims)];
    const float f11 = table[transmittanceIndex(y1, x1, iLambda, dims)];

    const float fx0 = lerpf(f00, f10, tx);
    const float fx1 = lerpf(f01, f11, tx);
    return lerpf(fx0, fx1, ty);
}

}  // namespace bake
