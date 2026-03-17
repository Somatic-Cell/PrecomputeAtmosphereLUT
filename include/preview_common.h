#pragma once

#include "asset_format.h"

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace preview {

struct PreviewMetadata {
    float bottomRadiusM = 0.0f;
    float topRadiusM = 0.0f;
    float observerAltitudeM = 0.0f;
    float muSMin = -0.2f;
    float miePhaseFunctionG = 0.8f;

    std::uint32_t skyMu = 0;
    std::uint32_t skyMuS = 0;
    std::uint32_t skyNu = 0;
    std::uint32_t irradianceR;
    std::uint32_t irradianceMuS;

    std::vector<float> wavelengthsNm;
};

struct FinalTable {
    atmo::TableHeader header{};
    std::vector<float> values;
};

struct FinalLuts {
    FinalTable sunTransmittance;   // dim0 = mu_s, dim1 = lambda
    FinalTable skyRayleighSingle;  // dim0 = nu, dim1 = mu, dim2 = mu_s, dim3 = lambda
    FinalTable skyMieSingle;       // dim0 = nu, dim1 = mu, dim2 = mu_s, dim3 = lambda
    FinalTable skyMultiple;        // dim0 = nu, dim1 = mu, dim2 = mu_s, dim3 = lambda
    FinalTable directIrradiance;   // dim0 = r, dim1 = mu_s, dim2 = lambda
};

struct RenderOptions {
    std::uint32_t width = 1024;
    std::uint32_t height = 512;
    float sunZenithDeg = 30.0f;
    float sunAzimuthDeg = 0.0f;
    float exposure = 10.0f;
    float gamma = 2.2f;
    bool includeRayleigh = true;
    bool includeMie = true;
    bool includeMultiple = true;
    bool blackBelowHorizon = true;
};

struct ImageRgb32f {
    std::uint32_t width = 0;
    std::uint32_t height = 0;
    std::vector<float> pixels; // RGBRGB...
};

bool loadMetadataJson(const std::filesystem::path& path,
                      PreviewMetadata& outMeta,
                      std::string& outError);

bool loadFinalLuts(const std::filesystem::path& outDir,
                   PreviewMetadata& inOutMeta,
                   FinalLuts& outLuts,
                   std::string& outError);

void renderEquirectangularSky(const PreviewMetadata& meta,
                              const FinalLuts& luts,
                              const RenderOptions& options,
                              ImageRgb32f& outImage);

bool savePng8(const std::filesystem::path& path,
              const ImageRgb32f& image,
              float exposure,
              float gamma,
              std::string& outError);

bool savePpm8(const std::filesystem::path& path,
              const ImageRgb32f& image,
              float exposure,
              float gamma,
              std::string& outError);

bool saveImage8(const std::filesystem::path& path,
                const ImageRgb32f& image,
                float exposure,
                float gamma,
                std::string& outError);

} // namespace preview
