#include "../src/bake_common.cuh"

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv) {
  try {
    std::filesystem::path configPath;
    std::filesystem::path outDir;

    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];
      if (arg == "--config" && i + 1 < argc) {
        configPath = argv[++i];
      } else if (arg == "--out" && i + 1 < argc) {
        outDir = argv[++i];
      }
    }

    if (configPath.empty() || outDir.empty()) {
      std::cerr << "Usage: AtmosphereBake --config config.json --out output_dir\n";
      return 1;
    }

    std::filesystem::create_directories(outDir);

    // ---------------------------------------------------------------------
    // 設定ファイルの読み込み
    // ---------------------------------------------------------------------
    bake::BakeConfig cfg = bake::loadConfigFromJson(configPath);

    // ---------------------------------------------------------------------
    // LUT の宣言
    // ---------------------------------------------------------------------
    std::vector<float> transmittance;
    std::vector<float> sunTransmittance;
    std::vector<float> skyRayleighSingle;
    std::vector<float> skyMieSingle;
    std::vector<float> directIrradiance;

    // 多重散乱の計算に使う中間表現用の LUT
    std::vector<float> deltaRayleighSingle;
    std::vector<float> deltaMieSingle;

    // 最終的な多重散乱の LUT
    std::vector<float> skyMultiple;

    // ---------------------------------------------------------------------
    // LUT の構築
    // ---------------------------------------------------------------------
    std::cout << "Building internal transmittance LUT...\n";
    bake::buildTransmittanceLUT(cfg, transmittance);

    buildDirectIrradiance(cfg, transmittance, directIrradiance);

    std::cout << "Extracting observer sun transmittance LUT...\n";
    bake::buildSunTransmittanceFromTransmittance(cfg, transmittance, sunTransmittance);

    std::cout << "Building sky single-scattering LUTs...\n";
    bake::buildSkySingleScatter(cfg, transmittance, skyRayleighSingle, skyMieSingle);

    // ---------------------------------------------------------------------
    // 多重散乱の計算
    // ---------------------------------------------------------------------
    const std::size_t finalSkySize = 
        static_cast<std::size_t>(cfg.lut.nMu) *
        static_cast<std::size_t>(cfg.lut.nMuS) *
        static_cast<std::size_t>(cfg.lut.nNu) *
        static_cast<std::size_t>(cfg.meta.wavelengths_nm.size());

    if (cfg.multiple.enabled && cfg.multiple.scatteringOrders > 1) {
        std::cout << "Building internal delta single-scattering LUTs...\n";
        bake::buildDeltaSingleScattering(
            cfg,
            transmittance,
            deltaRayleighSingle,
            deltaMieSingle);

        std::cout << "Building multiple-scattering sky LUT...\n";
        bake::buildMultipleScattering(
            cfg,
            transmittance,
            deltaRayleighSingle,
            deltaMieSingle,
            skyMultiple);
    } else {
        std::cout << "Multiple scattering disabled. Writing zero SkyMultiple LUT...\n";
        skyMultiple.assign(finalSkySize, 0.0f);
    }

    // ---------------------------------------------------------------------
    // メタデータの保存
    // ---------------------------------------------------------------------
    bake::saveMetadataJson(cfg, outDir.string().c_str());

    // ---------------------------------------------------------------------
    // 最終的な LUT の保存
    // ---------------------------------------------------------------------

    bake::saveTableBin((outDir / "transmittance_internal.bin"),
                           atmo::TableType::SunTransmittance,
                           static_cast<uint32_t>(cfg.lut.nMuS),
                           static_cast<uint32_t>(cfg.meta.wavelengths_nm.size()),
                           1,
                           1,
                           sunTransmittance);
    bake::saveTableBin(outDir / "direct_irradiance.bin",
                            atmo::TableType::DirectIrradiance, 
                            cfg.lut.irradianceR,
                            cfg.lut.irradianceMuS,
                            static_cast<std::uint32_t>(cfg.meta.wavelengths_nm.size()),
                            1u,
                            directIrradiance);

    bake::saveTableBin((outDir / "sky_rayleigh_single.bin"),
                           atmo::TableType::SkyRayleighSingle,
                           static_cast<uint32_t>(cfg.lut.nNu),
                           static_cast<uint32_t>(cfg.lut.nMu),
                           static_cast<uint32_t>(cfg.lut.nMuS),
                           static_cast<uint32_t>(cfg.meta.wavelengths_nm.size()),
                           skyRayleighSingle);
    bake::saveTableBin((outDir / "sky_mie_single.bin"),
                           atmo::TableType::SkyMieSingle,
                           static_cast<uint32_t>(cfg.lut.nNu),
                           static_cast<uint32_t>(cfg.lut.nMu),
                           static_cast<uint32_t>(cfg.lut.nMuS),
                           static_cast<uint32_t>(cfg.meta.wavelengths_nm.size()),
                           skyMieSingle);
    bake::saveTableBin((outDir / "sky_multiple.bin"),
                           atmo::TableType::SkyMultiple,
                           static_cast<uint32_t>(cfg.lut.nNu),
                           static_cast<uint32_t>(cfg.lut.nMu),
                           static_cast<uint32_t>(cfg.lut.nMuS),
                           static_cast<uint32_t>(cfg.meta.wavelengths_nm.size()),
                           skyMultiple);


    std::cout << "Done.\n";
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }
}
