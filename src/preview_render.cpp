#include "preview_common.h"
#include "cie1931_2deg_1nm.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>

#define FPNG_IMPLEMENTATION
#if __has_include("../ext/fpng/src/fpng.h")
  #include "../ext/fpng/src/fpng.h"
  #define PREVIEW_HAS_FPNG 1
#else
  #define PREVIEW_HAS_FPNG 0
#endif

namespace preview {
namespace {

constexpr float kPi = 3.14159265358979323846f;

struct Vec3 {
    float x;
    float y;
    float z;
};

Vec3 makeVec3(float x, float y, float z)
{
    Vec3 v{};
    v.x = x;
    v.y = y;
    v.z = z;
    return v;
}

Vec3 operator+(const Vec3& a, const Vec3& b)
{
    return makeVec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

Vec3& operator+=(Vec3& a, const Vec3& b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}

Vec3 operator*(float s, const Vec3& v)
{
    return makeVec3(s * v.x, s * v.y, s * v.z);
}

Vec3 operator*(const Vec3& v, float s)
{
    return s * v;
}

float dot(const Vec3& a, const Vec3& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 normalize(const Vec3& v)
{
    const float len2 = dot(v, v);
    if(len2 <= 0.0f) {
        return makeVec3(0.0f, 1.0f, 0.0f);
    }
    const float invLen = 1.0f / std::sqrt(len2);
    return invLen * v;
}

float clamp01(float x)
{
    return std::max(0.0f, std::min(1.0f, x));
}

float clampSigned(float x)
{
    return std::max(-1.0f, std::min(1.0f, x));
}

float lerp(float a, float b, float t)
{
    return a + (b - a) * t;
}

std::size_t skyIndex(const PreviewMetadata& meta,
                     std::uint32_t iMu,
                     std::uint32_t iMuS,
                     std::uint32_t iNu,
                     std::uint32_t iLambda)
{
    const std::size_t lambdaCount = meta.wavelengthsNm.size();
    return ((((static_cast<std::size_t>(iMu) * meta.skyMuS + iMuS) * meta.skyNu + iNu) * lambdaCount) + iLambda);
}

float phaseRayleigh(float nu)
{
    return (3.0f / (16.0f * kPi)) * (1.0f + nu * nu);
}

float phaseCornetteShanks(float nu, float g)
{
    const float g2 = g * g;
    const float denomBase = 1.0f + g2 - 2.0f * g * nu;
    const float denom = std::pow(std::max(1.0e-6f, denomBase), 1.5f);
    const float coeff = (3.0f / (8.0f * kPi)) * ((1.0f - g2) / (2.0f + g2));
    return coeff * (1.0f + nu * nu) / denom;
}

float sampleSky4D(const std::vector<float>& data,
                  const PreviewMetadata& meta,
                  float mu,
                  float muS,
                  float nu,
                  std::uint32_t iLambda)
{
    const float muT = 0.5f * (clampSigned(mu) + 1.0f);
    const float muST = clamp01((muS - meta.muSMin) / std::max(1.0e-6f, 1.0f - meta.muSMin));
    const float nuT = clamp01(0.5f * (nu + 1.0f));

    const float x = muT * static_cast<float>(meta.skyMu - 1);
    const float y = muST * static_cast<float>(meta.skyMuS - 1);
    const float z = nuT * static_cast<float>(meta.skyNu - 1);

    const std::uint32_t x0 = static_cast<std::uint32_t>(std::floor(x));
    const std::uint32_t x1 = std::min(x0 + 1, meta.skyMu - 1);
    const std::uint32_t y0 = static_cast<std::uint32_t>(std::floor(y));
    const std::uint32_t y1 = std::min(y0 + 1, meta.skyMuS - 1);
    const std::uint32_t z0 = static_cast<std::uint32_t>(std::floor(z));
    const std::uint32_t z1 = std::min(z0 + 1, meta.skyNu - 1);

    const float fx = x - static_cast<float>(x0);
    const float fy = y - static_cast<float>(y0);
    const float fz = z - static_cast<float>(z0);

    const float c000 = data[skyIndex(meta, x0, y0, z0, iLambda)];
    const float c100 = data[skyIndex(meta, x1, y0, z0, iLambda)];
    const float c010 = data[skyIndex(meta, x0, y1, z0, iLambda)];
    const float c110 = data[skyIndex(meta, x1, y1, z0, iLambda)];
    const float c001 = data[skyIndex(meta, x0, y0, z1, iLambda)];
    const float c101 = data[skyIndex(meta, x1, y0, z1, iLambda)];
    const float c011 = data[skyIndex(meta, x0, y1, z1, iLambda)];
    const float c111 = data[skyIndex(meta, x1, y1, z1, iLambda)];

    const float c00 = lerp(c000, c100, fx);
    const float c10 = lerp(c010, c110, fx);
    const float c01 = lerp(c001, c101, fx);
    const float c11 = lerp(c011, c111, fx);
    const float c0 = lerp(c00, c10, fy);
    const float c1 = lerp(c01, c11, fy);
    return lerp(c0, c1, fz);
}

float deltaLambda(const std::vector<float>& wavelengthsNm, std::size_t i)
{
    if(wavelengthsNm.empty()) {
        return 1.0f;
    }
    if(wavelengthsNm.size() == 1) {
        return 1.0f;
    }
    if(i == 0) {
        return wavelengthsNm[1] - wavelengthsNm[0];
    }
    if(i + 1 == wavelengthsNm.size()) {
        return wavelengthsNm[i] - wavelengthsNm[i - 1];
    }
    return 0.5f * (wavelengthsNm[i + 1] - wavelengthsNm[i - 1]);
}

Vec3 xyzToLinearSrgb(const Vec3& xyz)
{
    const float r =  3.2406f * xyz.x - 1.5372f * xyz.y - 0.4986f * xyz.z;
    const float g = -0.9689f * xyz.x + 1.8758f * xyz.y + 0.0415f * xyz.z;
    const float b =  0.0557f * xyz.x - 0.2040f * xyz.y + 1.0570f * xyz.z;
    return makeVec3(r, g, b);
}

float linearToSrgb(float x)
{
    x = std::max(0.0f, x);
    if(x <= 0.0031308f) {
        return 12.92f * x;
    }
    return 1.055f * std::pow(x, 1.0f / 2.4f) - 0.055f;
}

float luminanceLinearSrgb(const Vec3& rgb)
{
    return 0.2126f * rgb.x + 0.7152f * rgb.y + 0.0722f * rgb.z;
}

float acesFitted(float x)
{
    const float a = 2.51f;
    const float b = 0.03f;
    const float c = 2.43f;
    const float d = 0.59f;
    const float e = 0.14f;
    const float y = (x * (a * x + b)) / (x * (c * x + d) + e);
    return clamp01(y);
}

Vec3 toneMapFilmicPreserveChroma(const Vec3& hdr, float exposure)
{
    const Vec3 scaled = std::max(0.0f, exposure) * hdr;
    const float Y = luminanceLinearSrgb(scaled);
    if(Y <= 1.0e-6f) {
        return makeVec3(0.0f, 0.0f, 0.0f);
    }

    const float mappedY = acesFitted(Y);
    const float scale = mappedY / Y;

    return makeVec3(
        std::max(0.0f, scaled.x * scale),
        std::max(0.0f, scaled.y * scale),
        std::max(0.0f, scaled.z * scale));
}

Vec3 computeSkyRadianceLinearSrgb(const PreviewMetadata& meta,
                                  const FinalLuts& luts,
                                  float mu,
                                  float muS,
                                  float nu,
                                  const RenderOptions& options)
{
    Vec3 xyz = makeVec3(0.0f, 0.0f, 0.0f);

    const float phaseR = phaseRayleigh(nu);
    const float phaseM = phaseCornetteShanks(nu, meta.miePhaseFunctionG);

    for(std::uint32_t iLambda = 0; iLambda < static_cast<std::uint32_t>(meta.wavelengthsNm.size()); ++iLambda) {
        float spectral = 0.0f;
        if(options.includeRayleigh) {
            spectral += phaseR * sampleSky4D(luts.skyRayleighSingle.values, meta, mu, muS, nu, iLambda);
        }
        if(options.includeMie) {
            spectral += phaseM * sampleSky4D(luts.skyMieSingle.values, meta, mu, muS, nu, iLambda);
        }
        if(options.includeMultiple) {
            spectral += sampleSky4D(luts.skyMultiple.values, meta, mu, muS, nu, iLambda);
        }

        const float wavelengthNm = meta.wavelengthsNm[iLambda];
        const float w = deltaLambda(meta.wavelengthsNm, iLambda);
        const atmo::preview::CieXyzBar cmf = atmo::preview::sampleCie1931_2deg_1nm(wavelengthNm);
        xyz.x += spectral * cmf.xBar * w;
        xyz.y += spectral * cmf.yBar * w;
        xyz.z += spectral * cmf.zBar * w;
    }

    // Bruneton 2017 に寄せて、ここでは「スペクトル radiance -> XYZ -> linear sRGB」
    // のみを行う。sky LUT 自体が既に solar irradiance と太陽方向 transmittance
    // を含んでいる設計なので、preview で SunTransmittance を再度掛けない。
    
    Vec3 rgb = xyzToLinearSrgb(xyz);
    rgb.x = std::max(0.0f, rgb.x);
    rgb.y = std::max(0.0f, rgb.y);
    rgb.z = std::max(0.0f, rgb.z);
    return rgb;
}

Vec3 directionFromLatLong(float u, float v)
{
    const float phi = (u * 2.0f - 1.0f) * kPi;
    const float theta = v * kPi;

    const float sinTheta = std::sin(theta);
    const float cosTheta = std::cos(theta);
    const float sinPhi = std::sin(phi);
    const float cosPhi = std::cos(phi);

    return normalize(makeVec3(sinTheta * sinPhi, cosTheta, sinTheta * cosPhi));
}

Vec3 sunDirection(float zenithDeg, float azimuthDeg)
{
    const float theta = zenithDeg * (kPi / 180.0f);
    const float phi = azimuthDeg * (kPi / 180.0f);

    const float sinTheta = std::sin(theta);
    const float cosTheta = std::cos(theta);
    const float sinPhi = std::sin(phi);
    const float cosPhi = std::cos(phi);

    return normalize(makeVec3(sinTheta * sinPhi, cosTheta, sinTheta * cosPhi));
}

std::vector<std::uint8_t> makeRgb8Buffer(const ImageRgb32f& image, float exposure, float gamma)
{
    (void)gamma;
    std::vector<std::uint8_t> out(static_cast<std::size_t>(image.width) * image.height * 3, 0);
    for(std::size_t i = 0; i < static_cast<std::size_t>(image.width) * image.height; ++i) {
        const Vec3 hdr = makeVec3(
            image.pixels[i * 3 + 0],
            image.pixels[i * 3 + 1],
            image.pixels[i * 3 + 2]);

        const Vec3 mapped = toneMapFilmicPreserveChroma(hdr, exposure);

        const float r = linearToSrgb(clamp01(mapped.x));
        const float g = linearToSrgb(clamp01(mapped.y));
        const float b = linearToSrgb(clamp01(mapped.z));

        out[i * 3 + 0] = static_cast<std::uint8_t>(std::round(255.0f * clamp01(r)));
        out[i * 3 + 1] = static_cast<std::uint8_t>(std::round(255.0f * clamp01(g)));
        out[i * 3 + 2] = static_cast<std::uint8_t>(std::round(255.0f * clamp01(b)));
    }
    return out;
}

std::string lowerExt(const std::filesystem::path& path)
{
    std::string ext = path.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return ext;
}

std::size_t irradianceIndex(const PreviewMetadata& meta,
                            std::uint32_t iR,
                            std::uint32_t iMuS,
                            std::uint32_t iLambda)
{
    const std::size_t lambdaCount = meta.wavelengthsNm.size();
    return (((static_cast<std::size_t>(iR) * meta.irradianceMuS + iMuS) * lambdaCount) + iLambda);
}

float sampleIrradiance2D(const std::vector<float>& data,
                         const PreviewMetadata& meta,
                         float r_m,
                         float muS)
{
    const float bottom = meta.bottomRadiusM;
    const float top    = meta.topRadiusM;

    const float H   = std::sqrt(std::max(0.0f, top * top - bottom * bottom));
    const float rho = std::sqrt(std::max(0.0f, r_m * r_m - bottom * bottom));

    const float u = clamp01((muS - meta.muSMin) / std::max(1.0e-6f, 1.0f - meta.muSMin));
    const float v = (H > 0.0f) ? clamp01(rho / H) : 0.0f;

    const float x = u * static_cast<float>(meta.irradianceMuS - 1);
    const float y = v * static_cast<float>(meta.irradianceR - 1);

    const std::uint32_t x0 = static_cast<std::uint32_t>(std::floor(x));
    const std::uint32_t x1 = std::min(x0 + 1, meta.irradianceMuS - 1);
    const std::uint32_t y0 = static_cast<std::uint32_t>(std::floor(y));
    const std::uint32_t y1 = std::min(y0 + 1, meta.irradianceR - 1);

    const float fx = x - static_cast<float>(x0);
    const float fy = y - static_cast<float>(y0);

    auto fetch = [&](std::uint32_t iR, std::uint32_t iMuS, std::uint32_t iLambda) -> float {
        return data[irradianceIndex(meta, iR, iMuS, iLambda)];
    };

    const float f00 = fetch(y0, x0, 0); // dummy, not used outside lambda loop
    (void)f00; // warning suppression only
    return 0.0f;
}

float sampleIrradiance2D(const std::vector<float>& data,
                         const PreviewMetadata& meta,
                         float r_m,
                         float muS,
                         std::uint32_t iLambda)
{
    const float bottom = meta.bottomRadiusM;
    const float top    = meta.topRadiusM;

    const float H   = std::sqrt(std::max(0.0f, top * top - bottom * bottom));
    const float rho = std::sqrt(std::max(0.0f, r_m * r_m - bottom * bottom));

    const float u = clamp01((muS - meta.muSMin) / std::max(1.0e-6f, 1.0f - meta.muSMin));
    const float v = (H > 0.0f) ? clamp01(rho / H) : 0.0f;

    const float x = u * static_cast<float>(meta.irradianceMuS - 1);
    const float y = v * static_cast<float>(meta.irradianceR - 1);

    const std::uint32_t x0 = static_cast<std::uint32_t>(std::floor(x));
    const std::uint32_t x1 = std::min(x0 + 1, meta.irradianceMuS - 1);
    const std::uint32_t y0 = static_cast<std::uint32_t>(std::floor(y));
    const std::uint32_t y1 = std::min(y0 + 1, meta.irradianceR - 1);

    const float fx = x - static_cast<float>(x0);
    const float fy = y - static_cast<float>(y0);

    const float f00 = data[irradianceIndex(meta, y0, x0, iLambda)];
    const float f10 = data[irradianceIndex(meta, y0, x1, iLambda)];
    const float f01 = data[irradianceIndex(meta, y1, x0, iLambda)];
    const float f11 = data[irradianceIndex(meta, y1, x1, iLambda)];

    const float fx0 = lerp(f00, f10, fx);
    const float fx1 = lerp(f01, f11, fx);
    return lerp(fx0, fx1, fy);
}

bool intersectGroundSphere(const PreviewMetadata& meta,
                           const Vec3& viewDir,
                           Vec3& outGroundPos,
                           Vec3& outGroundNormal)
{
    const float rObserver = meta.bottomRadiusM + meta.observerAltitudeM;
    const float R = meta.bottomRadiusM;

    const Vec3 origin = makeVec3(0.0f, rObserver, 0.0f);

    // |origin + t * dir|^2 = R^2
    const float b = dot(origin, viewDir);
    const float c = dot(origin, origin) - R * R;
    const float disc = b * b - c;
    if(disc < 0.0f) {
        return false;
    }

    const float s = std::sqrt(std::max(0.0f, disc));
    const float t0 = -b - s;
    const float t1 = -b + s;
    const float t = (t0 > 0.0f) ? t0 : t1;
    if(t <= 0.0f) {
        return false;
    }

    outGroundPos = makeVec3(
        origin.x + t * viewDir.x,
        origin.y + t * viewDir.y,
        origin.z + t * viewDir.z);

    outGroundNormal = normalize(outGroundPos);
    return true;
}

Vec3 computeGroundRadianceLinearSrgb(const PreviewMetadata& meta,
                                     const FinalLuts& luts,
                                     const Vec3& groundNormal,
                                     const Vec3& sunDir)
{
    // まずは無彩色の Lambert 地面
    constexpr float kGroundAlbedo = 0.03f;

    const float muSLocal = std::max(0.0f, dot(groundNormal, sunDir));
    if(muSLocal <= 0.0f) {
        return makeVec3(0.0f, 0.0f, 0.0f);
    }

    Vec3 xyz = makeVec3(0.0f, 0.0f, 0.0f);

    for(std::uint32_t iLambda = 0; iLambda < static_cast<std::uint32_t>(meta.wavelengthsNm.size()); ++iLambda) {
        const float E_direct =
            sampleIrradiance2D(luts.directIrradiance.values,
                               meta,
                               meta.bottomRadiusM,
                               muSLocal,
                               iLambda);

        const float spectral = E_direct * muSLocal * (kGroundAlbedo / kPi);
        
        const float wavelengthNm = meta.wavelengthsNm[iLambda];
        const float w = deltaLambda(meta.wavelengthsNm, iLambda);
        const atmo::preview::CieXyzBar cmf = atmo::preview::sampleCie1931_2deg_1nm(wavelengthNm);

        xyz.x += spectral * cmf.xBar * w;
        xyz.y += spectral * cmf.yBar * w;
        xyz.z += spectral * cmf.zBar * w;
    }

    Vec3 rgb = xyzToLinearSrgb(xyz);
    rgb.x = std::max(0.0f, rgb.x);
    rgb.y = std::max(0.0f, rgb.y);
    rgb.z = std::max(0.0f, rgb.z);
    return rgb;
}

} // namespace

void renderEquirectangularSky(const PreviewMetadata& meta,
                              const FinalLuts& luts,
                              const RenderOptions& options,
                              ImageRgb32f& outImage)
{
    outImage.width = options.width;
    outImage.height = options.height;
    outImage.pixels.assign(static_cast<std::size_t>(options.width) * options.height * 3, 0.0f);

    const Vec3 sunDir = sunDirection(options.sunZenithDeg, options.sunAzimuthDeg);
    const float muSObserver = sunDir.y;

    for(std::uint32_t y = 0; y < options.height; ++y) {
        for(std::uint32_t x = 0; x < options.width; ++x) {
            const float u = (static_cast<float>(x) + 0.5f) / static_cast<float>(options.width);
            const float v = (static_cast<float>(y) + 0.5f) / static_cast<float>(options.height);
            const Vec3 viewDir = directionFromLatLong(u, v);

            const float mu = viewDir.y;
            const float nu = dot(viewDir, sunDir);

            Vec3 rgb = computeSkyRadianceLinearSrgb(meta, luts, mu, muSObserver, nu, options);

            if(mu < 0.0f) {
                Vec3 groundPos{};
                Vec3 groundNormal{};
                if(intersectGroundSphere(meta, viewDir, groundPos, groundNormal)) {
                    rgb += computeGroundRadianceLinearSrgb(meta, luts, groundNormal, sunDir);
                }
            }

            const std::size_t idx = (static_cast<std::size_t>(y) * options.width + x) * 3;
            outImage.pixels[idx + 0] = rgb.x;
            outImage.pixels[idx + 1] = rgb.y;
            outImage.pixels[idx + 2] = rgb.z;
        }

    }
}

bool savePng8(const std::filesystem::path& path,
              const ImageRgb32f& image,
              float exposure,
              float gamma,
              std::string& outError)
{
#if PREVIEW_HAS_FPNG
    const std::vector<std::uint8_t> rgb8 = makeRgb8Buffer(image, exposure, gamma);
    fpng::fpng_init();
    if(!fpng::fpng_encode_image_to_file(path.string().c_str(),
                                        rgb8.data(),
                                        image.width,
                                        image.height,
                                        3)) {
        outError = "fpng failed to encode PNG: " + path.string();
        return false;
    }
    return true;
#else
    (void)image;
    (void)exposure;
    (void)gamma;
    outError = "savePng8() was requested, but fpng.h was not available at build time.";
    return false;
#endif
}

bool savePpm8(const std::filesystem::path& path,
              const ImageRgb32f& image,
              float exposure,
              float gamma,
              std::string& outError)
{
    std::ofstream ofs(path, std::ios::binary);
    if(!ofs) {
        outError = "Failed to open output file: " + path.string();
        return false;
    }

    ofs << "P6\n" << image.width << ' ' << image.height << "\n255\n";
    const std::vector<std::uint8_t> rgb8 = makeRgb8Buffer(image, exposure, gamma);
    ofs.write(reinterpret_cast<const char*>(rgb8.data()), static_cast<std::streamsize>(rgb8.size()));
    if(!ofs) {
        outError = "Failed to write PPM file: " + path.string();
        return false;
    }
    return true;
}

bool saveImage8(const std::filesystem::path& path,
                const ImageRgb32f& image,
                float exposure,
                float gamma,
                std::string& outError)
{
    const std::string ext = lowerExt(path);
    if(ext == ".png") {
        return savePng8(path, image, exposure, gamma, outError);
    }
    if(ext == ".ppm" || ext.empty()) {
        return savePpm8(path, image, exposure, gamma, outError);
    }

    outError = "Unsupported output image extension: " + path.extension().string();
    return false;
}

} // namespace preview
