#pragma once
#include <cmath>
#include "meta_data.h"

#ifdef __CUDACC__
  #define ATMO_HD __host__ __device__
#else
  #define ATMO_HD
#endif

namespace atmo {

constexpr float kPi = 3.14159265358979323846f;

ATMO_HD inline float clampf(float x, float lo, float hi) {
  return x < lo ? lo : (x > hi ? hi : x);
}

ATMO_HD inline float sqrf(float x) { return x * x; }
ATMO_HD inline float safeSqrt(float x) { return ::sqrtf(x > 0.0f ? x : 0.0f); }
ATMO_HD inline float clampCosine(float mu) { return clampf(mu, -1.0f, 1.0f); }
ATMO_HD inline float clampDistance(float d) { return d > 0.0f ? d : 0.0f; }

ATMO_HD inline float distanceToTopAtmosphereBoundary(
    float top_radius_m, float r_m, float mu) {
  const float disc = r_m * r_m * (mu * mu - 1.0f) + top_radius_m * top_radius_m;
  return clampDistance(-r_m * mu + safeSqrt(disc));
}

ATMO_HD inline float distanceToBottomAtmosphereBoundary(
    float bottom_radius_m, float r_m, float mu) {
  const float disc = r_m * r_m * (mu * mu - 1.0f) + bottom_radius_m * bottom_radius_m;
  return clampDistance(-r_m * mu - safeSqrt(disc));
}

ATMO_HD inline bool rayIntersectsGround(
    float bottom_radius_m, float r_m, float mu) {
  return mu < 0.0f &&
         (r_m * r_m * (mu * mu - 1.0f) + bottom_radius_m * bottom_radius_m) >= 0.0f;
}

ATMO_HD inline float distanceToNearestAtmosphereBoundary(
    float bottom_radius_m, float top_radius_m, float r_m, float mu, bool intersects_ground) {
  return intersects_ground ? distanceToBottomAtmosphereBoundary(bottom_radius_m, r_m, mu)
                           : distanceToTopAtmosphereBoundary(top_radius_m, r_m, mu);
}

ATMO_HD inline float pointRadiusAlongRay(float r_m, float mu, float d_m) {
  return safeSqrt(d_m * d_m + 2.0f * r_m * mu * d_m + r_m * r_m);
}

ATMO_HD inline float evalLayerDensity(const DensityLayer& layer, float altitude_m) {
  const float density =
      layer.expTerm * ::expf(layer.expScale * altitude_m) +
      layer.linearTerm * altitude_m +
      layer.constantTerm;
  return clampf(density, 0.0f, 1.0f);
}

ATMO_HD inline float evalProfileDensity(const DensityProfile& profile, float altitude_m) {
  return altitude_m < profile.layers[0].width_m
      ? evalLayerDensity(profile.layers[0], altitude_m)
      : evalLayerDensity(profile.layers[1], altitude_m);
}

ATMO_HD inline float rayleighPhaseFunction(float nu) {
  const float k = 3.0f / (16.0f * kPi);
  return k * (1.0f + nu * nu);
}

ATMO_HD inline float miePhaseFunction(float g, float nu) {
  const float k = 3.0f / (8.0f * kPi) * (1.0f - g * g) / (2.0f + g * g);
  const float denom = 1.0f + g * g - 2.0f * g * nu;
  return k * (1.0f + nu * nu) / ::powf(denom, 1.5f);
}

}  // namespace atmo
