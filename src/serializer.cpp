#include "bake_common.cuh"
#include "config_io.h"

#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

#include <nlohmann/json.hpp>

namespace bake {
std::uint64_t computeElementCount(
    uint32_t d0,
    uint32_t d1,
    uint32_t d2,
    uint32_t d3
)
{
    return  static_cast<std::uint64_t>(d0) * 
            static_cast<std::uint64_t>(d1) * 
            static_cast<std::uint64_t>(d2) * 
            static_cast<std::uint64_t>(d3);
}

void saveTableBin(
    const std::filesystem::path& path,
    atmo::TableType type,
    uint32_t d0,
    uint32_t d1,
    uint32_t d2,
    uint32_t d3,
    const std::vector<float>& data      // データ本体
) 
{
    // -----------------------
    // データサイズのチェックと header 部の作成
    // -----------------------

    // データの次元が一致しているかどうかを確認
    const std::uint64_t expectedCount   = computeElementCount(d0, d1, d2, d3);
    const std::uint64_t actualCount     = static_cast<std::uint64_t>(data.size());

    if(expectedCount != actualCount){
        throw std::runtime_error(
            "save_table_bin: dimension/product mismatch for file '" +
            path.string() + "' (expected: " +
            std::to_string(expectedCount) + "floats, but got " + 
            std::to_string(actualCount) + ".)"
        ); 
    }

    const std::uint64_t byteCount = actualCount * static_cast<std::uint64_t>(sizeof(float));
    if(byteCount > static_cast<std::uint64_t>(std::numeric_limits<std::streamsize>::max())){
        throw std::runtime_error(
            "saveTableBin: payload too large for std::streamsize: " + 
            path.string()
        );
    }

    // header 部の作成
    atmo::TableHeader header{};
    header.tableType = type;
    header.storageFormat = atmo::StorageFormat::FP32;
    header.dim0 = d0;
    header.dim1 = d1;
    header.dim2 = d2;
    header.dim3 = d3;

    // -----------------------
    // ファイルの作成
    // -----------------------
    std::ofstream ofs(path, std::ios::binary);
    if (!ofs) {
        throw std::runtime_error(std::string("failed to open output file: ") + path.string());
    }

    // 書き込み
    // header 部分
    ofs.write(reinterpret_cast<const char*>(&header), sizeof(header));
    if(!ofs){
        throw std::runtime_error(
            "saveTableBin: failed to write header: " + path.string()
        );
    }

    // 本体
    ofs.write(reinterpret_cast<const char*>(data.data()),
        static_cast<std::streamsize>(data.size() * sizeof(float)));
    if(!ofs){
        throw std::runtime_error(
            "saveTableBin: failed to write payload: " + path.string()
        );
    }
}

nlohmann::json densityLayerToJson(const atmo::DensityLayer& l) {
    return {
        {"width_m", l.width_m},
        {"exp_term", l.expTerm},
        {"exp_scale", l.expScale},
        {"linear_term", l.linearTerm},
        {"constant_term", l.constantTerm}
    };
}

nlohmann::json densityProfileToJson(const atmo::DensityProfile& p) {
    return {
        {"layers", {
            densityLayerToJson(p.layers[0]),
            densityLayerToJson(p.layers[1])
        }}
    };
}


void saveMetadataJson(const BakeConfig& cfg, const char* outDir) {
    namespace fs = std::filesystem;
    using json = nlohmann::json;

    const fs::path dir(outDir);
    fs::create_directories(dir);
    const fs::path path = dir / "metadata.json";

    const auto& m = cfg.meta;

    json j;
    j["metadata_version"] = m.metaDataVersion;

    j["planet"] = {
        {"bottom_radius_m", m.bottomRadius_m},
        {"top_radius_m", m.topRadius_m},
        {"observer_altitude_m", m.observerAltitude_m},
        {"mu_s_min", m.muSMin}
    };

    j["wavelengths_nm"] = m.wavelengths_nm;

    j["rayleigh_density"]  = densityProfileToJson(m.rayleighDensity);
    j["mie_density"]       = densityProfileToJson(m.mieDensity);
    j["absorption_density"] = densityProfileToJson(m.absorptionDensity);

    j["solar_irradiance"]       = m.solarIrradiance;
    j["rayleigh_scattering"]    = m.rayleighScattering;
    j["mie_scattering"]         = m.mieScattering;
    j["mie_extinction"]         = m.mieExtinction;
    j["absorption_extinction"]  = m.absorptionExtinction;

    j["mie_phase_function"] = {
        {"type", "cornette_shanks"},
        {"g", m.miePhaseFunctionG}
    };

    j["multiple_scattering"] = {
        {"scattering_orders", m.scatteringOrders}
    };

    std::ofstream ofs(path);
    if (!ofs) {
        throw std::runtime_error("Failed to open metadata output file: " + path.string());
    }
    ofs << j.dump(2) << '\n';
    if (!ofs) {
        throw std::runtime_error("Failed to write metadata output file: " + path.string());
    }
}

}  // namespace bake
