#pragma once

namespace atmo {

inline float encodeMu(float mu) {
  return mu;
}

inline float decodeMu(float u) {
  return u;
}

inline float encodeMuS(float mu_s, float mu_s_min) {
  return (mu_s - mu_s_min) / (1.0f - mu_s_min);
}

inline float decodeMuS(float u, float mu_s_min) {
  return mu_s_min + u * (1.0f - mu_s_min);
}

inline float encodeNu(float nu) {
  return 0.5f * (nu + 1.0f);
}

inline float decodeNu(float u) {
  return 2.0f * u - 1.0f;
}

}  // namespace atmo
