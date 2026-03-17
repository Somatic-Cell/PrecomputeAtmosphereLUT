#pragma once

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "meta_data.h"

namespace bake {

// 波長の表現方法
enum class WavelengthMode {
    UniformStep,        // min, max, step   による表現
    UniformCount,       // min, max, count  による表現
    Explicit            // 明示的な配列
};

enum class DataSourceMode {
    Builtin,
    FileCsv
};

struct CsvColumnSpec {
    std::string path;
    std::string lambdaColumn = "wavelength_nm";
    std::string valueColumn;
};

struct CsvMieSpec {
    std::string path;
    std::string lambdaColumn = "wavelength_nm";
    std::string scatteringColumn = "mie_scattering";
    std::string extinctionColumn = "mie_extinction";
};

struct WavelengthSpec {
    WavelengthMode mode = WavelengthMode::UniformStep;

    float minNm = 360.0f;
    float maxNm = 830.0f;

    float stepNm = 10.0f;
    int count = 48;
    std::vector<float> explicitValuesNm;
};

struct PlanetConfig {
    // 惑星と大気の形状
    float bottomRadius_m = 6360000.0f;  // 惑星の半径
    float topRadius_m = 6460000.0f;     // 大気の上端
    float observerAltitude_m = 0.0f;    // 観測する高度
    
    // 太陽の天頂方向の余弦が取れる最小値
    float muSMin = -0.2f;
};

// 大気の密度モデル (2017 年の改良版)
struct DensityLayerConfig {
    float width_m = 0.0f;       // この層が有効な高度の区間の幅
    float expTerm = 0.0f;       // 指数項の係数
    float expScale = 0.0f;      // 指数項のスケール
    float linearTerm = 0.0f;    // 線形項のスケール
    float constantTerm = 0.0f;  // 定数項
};

// 2区間の密度表現 (Rayleigh, Mie では 1 区間， Ozone では 2区間を使う)
struct DensityProfileConfig {
    DensityLayerConfig layers[2];
};


struct SolarIrradianceConfig {
    DataSourceMode mode = DataSourceMode::Builtin;
    std::optional<CsvColumnSpec> csv;
};

struct RayleighConfig {
    DensityProfileConfig density;
    DataSourceMode scatteringMode = DataSourceMode::Builtin;
    std::optional<CsvColumnSpec> scatteringCsv;
};

struct MiePhaseFunctionConfig {
    std::string type = "cornette_shanks";
    float g = 0.8f;
};

struct MieConfig {
    DensityProfileConfig density;
    MiePhaseFunctionConfig phase;
    DataSourceMode coefficientsMode = DataSourceMode::Builtin;
    std::optional<CsvMieSpec> coefficientsCsv;
};

struct AbsorptionConfig {
    bool enabled = false;
    DensityProfileConfig density;
    DataSourceMode extinctionMode = DataSourceMode::Builtin;
    std::optional<CsvColumnSpec> extinctionCsv;
};

struct LutConfig {
    int nMu = 64;
    int nMuS = 64;
    int nNu = 128;

    int transmittanceR = 256;
    int transmittanceMu = 128;
    
    int irradianceR = 64;
    int irradianceMuS = 16;

    int integrationStepsTransmittance = 128;
    int integrationStepsView = 128;
    int integrationStepsSun = 64;
};

struct GroundConfig {
    DataSourceMode albedoMode = DataSourceMode::Builtin;
    std::optional<CsvColumnSpec> albedoCsv;
    float constantAlbedo = 0.1f;
};

struct MultipleScatteringConfig {
    bool enabled = true;        // 多重散乱を考慮するかどうか
    int scatteringOrders = 7;   // 何次散乱まで考慮するか
};

struct BakeConfig {
    std::filesystem::path configPath;
    std::filesystem::path configDir;

    PlanetConfig            planet;
    WavelengthSpec          wavelength;
    SolarIrradianceConfig   solar;
    RayleighConfig          rayleigh;
    MieConfig               mie;
    AbsorptionConfig        absorption;
    GroundConfig            ground;
    LutConfig               lut;
    MultipleScatteringConfig multiple;

    atmo::AtmosphereMetadata meta;
};

BakeConfig loadConfigFromJson(const std::filesystem::path& path);
std::vector<float> buildWavelengthGrid_nm(const WavelengthSpec& spec);
void populateMetadataFromConfig(BakeConfig& cfg);

} // namespace bake