#pragma once
#include <string>
#include <vector>
#include "asset_format.h"
#include "meta_data.h"
#include "config_io.h"

namespace bake {


void saveMetadataJson(const BakeConfig& cfg, const char* out_dir);

void saveTableBin(const std::filesystem::path& path,
                    atmo::TableType type,
                    uint32_t d0,
                    uint32_t d1,
                    uint32_t d2,
                    uint32_t d3,
                    const std::vector<float>& data);

void buildTransmittanceLUT(const BakeConfig& cfg, std::vector<float>& out_transmittance);
void buildSunTransmittanceFromTransmittance(const BakeConfig& cfg,
                                                const std::vector<float>& transmittance,
                                                std::vector<float>& outSunTransmittance);
void buildSkySingleScatter(const BakeConfig& cfg,
                              const std::vector<float>& transmittance,
                              std::vector<float>& out_rayleigh,
                              std::vector<float>& out_mie);
void buildDirectIrradiance(
    const BakeConfig& cfg,
    const std::vector<float>& transmittance,
    std::vector<float>& outDirectIrradiance);

void buildDeltaSingleScattering(
    const BakeConfig& cfg,
    const std::vector<float>& transmittance,
    std::vector<float>& outDeltaRayleigh,
    std::vector<float>& outDeltaMie);

void buildMultipleScattering(
    const BakeConfig& cfg,
    const std::vector<float>& transmittance,
    const std::vector<float>& deltaRayleighSingle,
    const std::vector<float>& deltaMieSingle,
    std::vector<float>& outSkyMultiple);

}  // namespace bake
