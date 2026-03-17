#include "config_io.h"
#include "device_common.cuh"

#include <vector>

namespace bake {

inline float decodeMuSFromUnit(float u, float muSMin) {
    return muSMin + u * (1.0f - muSMin);
}

void buildSunTransmittanceFromTransmittance(
  const BakeConfig& cfg,
  const std::vector<float>& transmittance,
  std::vector<float>& outSunTransmittance
) {
  DeviceLutDims dims{};
  dims.transmittanceR   = cfg.lut.transmittanceR;
  dims.transmittanceMu  = cfg.lut.transmittanceMu;
  dims.wavelengthCount  = static_cast<int>(cfg.meta.wavelengths_nm.size());

  const DeviceAtmosphereParameters atm = makeDeviceAtmosphericParameters(cfg.meta);

  const int nLambda = dims.wavelengthCount;
  const int nMuS    = cfg.lut.nMuS;

  outSunTransmittance.resize(
    static_cast<std::size_t>(nMuS) * static_cast<std::size_t>(nLambda)
  );

  const float rObserver = cfg.meta.bottomRadius_m + cfg.meta.observerAltitude_m;
  for (int iLambda = 0; iLambda < nLambda; ++iLambda) {
    for (int iMuS = 0; iMuS < nMuS; ++iMuS) {
      const float u = (static_cast<float>(iMuS) + 0.5f) / static_cast<float>(nMuS);
      const float muS = decodeMuSFromUnit(u, cfg.meta.muSMin);

      const float T = sampleTransmittanceLUT(
        transmittance.data(),
        dims,
        atm,
        iLambda,
        rObserver,
        muS);

      outSunTransmittance[
        static_cast<std::size_t>(iMuS) *
        static_cast<std::size_t>(nLambda) +
        static_cast<std::size_t>(iLambda)] = T;
    }
  }
}


}  // namespace bake
