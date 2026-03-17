#include "config_io.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "csv_table_io.h"

namespace bake {

using json = nlohmann::json;

namespace {


// ------------------------------------------------------------
// Small helpers
// ------------------------------------------------------------

template <typename T>
T require(const json& j, const char* key) {
    if (!j.contains(key)) {
        throw std::runtime_error(std::string("Missing required key: ") + key);
    }
    try {
        return j.at(key).get<T>();
    } catch (const std::exception& e) {
        throw std::runtime_error(
            std::string("Invalid type/value for key '") + key + "': " + e.what());
    }
}

template <typename T>
T optionalValue(const json& j, const char* key, const T& defaultValue) {
    if (!j.contains(key)) {
        return defaultValue;
    }
    try {
        return j.at(key).get<T>();
    } catch (const std::exception& e) {
        throw std::runtime_error(
            std::string("Invalid type/value for key '") + key + "': " + e.what());
    }
}

const json* findObject(const json& j, const char* key) {
    if (!j.contains(key)) {
        return nullptr;
    }
    const json& child = j.at(key);
    if (!child.is_object()) {
        throw std::runtime_error(std::string("Key must be an object: ") + key);
    }
    return &child;
}

std::string requireString(const json& j, const char* key) {
    return require<std::string>(j, key);
}

float requirePositiveFloat(const json& j, const char* key) {
    const float v = require<float>(j, key);
    if (!(v > 0.0f)) {
        throw std::runtime_error(std::string("Expected positive value for key: ") + key);
    }
    return v;
}

int requirePositiveInt(const json& j, const char* key) {
    const int v = require<int>(j, key);
    if (v <= 0) {
        throw std::runtime_error(std::string("Expected positive integer for key: ") + key);
    }
    return v;
}

// ------------------------------------------------------------
// Enum parsers
// ------------------------------------------------------------

DataSourceMode parseDataSourceMode(const std::string& s) {
    if (s == "builtin")  return DataSourceMode::Builtin;
    if (s == "file_csv") return DataSourceMode::FileCsv;
    throw std::runtime_error("Unknown data source mode: " + s);
}

WavelengthMode parseWavelengthMode(const std::string& s) {
    if (s == "uniform_step")  return WavelengthMode::UniformStep;
    if (s == "uniform_count") return WavelengthMode::UniformCount;
    if (s == "explicit")      return WavelengthMode::Explicit;
    throw std::runtime_error("Unknown wavelength mode: " + s);
}

// ------------------------------------------------------------
// CSV specs
// ------------------------------------------------------------

CsvColumnSpec parseCsvColumnSpec(const json& j) {
    CsvColumnSpec out;
    out.path = requireString(j, "path");
    out.lambdaColumn = optionalValue<std::string>(j, "lambda_column", "wavelength_nm");
    out.valueColumn = requireString(j, "value_column");
    return out;
}

CsvMieSpec parseCsvMieSpec(const json& j) {
    CsvMieSpec out;
    out.path = requireString(j, "path");
    out.lambdaColumn = optionalValue<std::string>(j, "lambda_column", "wavelength_nm");
    out.scatteringColumn = optionalValue<std::string>(j, "scattering_column", "mie_scattering");
    out.extinctionColumn = optionalValue<std::string>(j, "extinction_column", "mie_extinction");
    return out;
}

// ------------------------------------------------------------
// Density profile helpers
// ------------------------------------------------------------

DensityProfileConfig makeZeroProfileConfig() {
    DensityProfileConfig p{};
    p.layers[0] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    p.layers[1] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    return p;
}

DensityProfileConfig makeExponentialProfileConfig(float scaleHeightM) {
    if (!(scaleHeightM > 0.0f)) {
        throw std::runtime_error("scale_height_m must be positive");
    }

    DensityProfileConfig p{};
    p.layers[0] = {
        1.0e9f,
        1.0f,
        -1.0f / scaleHeightM,
        0.0f,
        0.0f
    };
    p.layers[1] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    return p;
}

DensityProfileConfig makeBuiltinOzoneProfileConfig() {
    DensityProfileConfig p{};
    p.layers[0] = {
        25000.0f,
        0.0f,
        0.0f,
        1.0f / 15000.0f,
        -2.0f / 3.0f
    };
    p.layers[1] = {
        1.0e9f,
        0.0f,
        0.0f,
        -1.0f / 15000.0f,
        8.0f / 3.0f
    };
    return p;
}

DensityLayerConfig parseDensityLayerConfig(const json& j) {
    DensityLayerConfig L;
    L.width_m      = optionalValue<float>(j, "width_m", 0.0f);
    L.expTerm      = optionalValue<float>(j, "exp_term", 0.0f);
    L.expScale     = optionalValue<float>(j, "exp_scale", 0.0f);
    L.linearTerm   = optionalValue<float>(j, "linear_term", 0.0f);
    L.constantTerm = optionalValue<float>(j, "constant_term", 0.0f);
    return L;
}

DensityProfileConfig parseDensityProfileConfig(
    const json& j,
    const DensityProfileConfig& fallback) {

    if (j.contains("layers")) {
        const json& layers = j.at("layers");
        if (!layers.is_array()) {
            throw std::runtime_error("'layers' must be an array");
        }
        if (layers.empty() || layers.size() > 2) {
            throw std::runtime_error("'layers' must contain 1 or 2 entries");
        }

        DensityProfileConfig out{};
        out.layers[0] = parseDensityLayerConfig(layers.at(0));
        out.layers[1] = (layers.size() >= 2)
            ? parseDensityLayerConfig(layers.at(1))
            : DensityLayerConfig{};
        return out;
    }

    if (j.contains("type")) {
        const std::string type = requireString(j, "type");
        if (type == "exp") {
            const float H = requirePositiveFloat(j, "scale_height_m");
            return makeExponentialProfileConfig(H);
        }
        if (type == "zero") {
            return makeZeroProfileConfig();
        }
        if (type == "builtin_ozone") {
            return makeBuiltinOzoneProfileConfig();
        }
        throw std::runtime_error("Unknown density profile type: " + type);
    }

    return fallback;
}

atmo::DensityProfile toAtmoDensityProfile(const DensityProfileConfig& src) {
    atmo::DensityProfile dst{};

    for (int i = 0; i < 2; ++i) {
        dst.layers[i].width_m       = src.layers[i].width_m;
        dst.layers[i].expTerm      = src.layers[i].expTerm;
        dst.layers[i].expScale     = src.layers[i].expScale;
        dst.layers[i].linearTerm   = src.layers[i].linearTerm;
        dst.layers[i].constantTerm = src.layers[i].constantTerm;
    }

    return dst;
}

// ------------------------------------------------------------
// Wavelength helpers
// ------------------------------------------------------------

void validateStrictlyIncreasing(const std::vector<float>& xs, const char* name) {
    if (xs.empty()) {
        throw std::runtime_error(std::string(name) + " must not be empty");
    }
    for (size_t i = 1; i < xs.size(); ++i) {
        if (!(xs[i] > xs[i - 1])) {
            throw std::runtime_error(std::string(name) + " must be strictly increasing");
        }
    }
}

// ------------------------------------------------------------
// Builtin spectra
// ------------------------------------------------------------

float planckRelative(float wavelengthNm, float temperatureK) {
    constexpr double h = 6.62607015e-34;
    constexpr double c = 2.99792458e8;
    constexpr double k = 1.380649e-23;

    const double lambda = static_cast<double>(wavelengthNm) * 1.0e-9;
    const double c1 = 2.0 * h * c * c;
    const double c2 = h * c / k;
    const double denom = std::exp(c2 / (lambda * temperatureK)) - 1.0;
    const double B = c1 / (std::pow(lambda, 5.0) * denom);
    return static_cast<float>(B);
}

std::vector<float> makeBuiltinSolarIrradiance(const std::vector<float>& wavelengthsNm) {
    std::vector<float> out;
    out.reserve(wavelengthsNm.size());

    float maxV = 0.0f;
    for (float lam : wavelengthsNm) {
        const float v = planckRelative(lam, 5778.0f);
        out.push_back(v);
        maxV = std::max(maxV, v);
    }

    if (maxV > 0.0f) {
        for (float& v : out) {
            v /= maxV;
        }
    }
    return out;
}

std::vector<float> makeBuiltinRayleighScattering(const std::vector<float>& wavelengthsNm) {
    constexpr float lambda0Nm = 680.0f;
    constexpr float beta0 = 5.802e-6f;

    std::vector<float> out;
    out.reserve(wavelengthsNm.size());

    for (float lam : wavelengthsNm) {
        out.push_back(beta0 * std::pow(lambda0Nm / lam, 4.0f)); // 波長の4乗分の1で消散係数が変化
    }
    return out;
}

void makeBuiltinMieCoefficients(
    const std::vector<float>& wavelengthsNm,
    std::vector<float>& scattering,
    std::vector<float>& extinction) {

    constexpr float betaExt550 = 3.996e-6f;
    constexpr float singleScatteringAlbedo = 0.9f;
    constexpr float alpha = 0.0f;

    scattering.clear();
    extinction.clear();
    scattering.reserve(wavelengthsNm.size());
    extinction.reserve(wavelengthsNm.size());

    for (float lam : wavelengthsNm) {
        const float betaExt = betaExt550 * std::pow(550.0f / lam, alpha);
        extinction.push_back(betaExt);
        scattering.push_back(singleScatteringAlbedo * betaExt);
    }
}

std::vector<float> makeBuiltinAbsorptionExtinction(
    const std::vector<float>& wavelengthsNm,
    bool enabled) {
    std::vector<float> out(wavelengthsNm.size(), 0.0f);

    if (!enabled) {
        return out;
    }

    for (size_t i = 0; i < wavelengthsNm.size(); ++i) {
        const float lam = wavelengthsNm[i];
        const float d = (lam - 600.0f) / 80.0f;
        out[i] = 1.5e-6f * std::exp(-0.5f * d * d);
    }

    return out;
}

// ------------------------------------------------------------
// Validation
// ------------------------------------------------------------

void validateCommonRanges(const BakeConfig& cfg) {
    if (!(cfg.planet.topRadius_m > cfg.planet.bottomRadius_m)) {
        throw std::runtime_error("planet.top_radius_m must be greater than planet.bottom_radius_m");
    }

    if (cfg.planet.observerAltitude_m < 0.0f) {
        throw std::runtime_error("planet.observer_altitude_m must be >= 0");
    }

    if (cfg.planet.bottomRadius_m + cfg.planet.observerAltitude_m >= cfg.planet.topRadius_m) {
        throw std::runtime_error("Observer must lie inside the atmosphere shell");
    }

    if (cfg.lut.nMu <= 0 || cfg.lut.nMuS <= 0 || cfg.lut.nNu <= 0) {
        throw std::runtime_error("Final sky LUT dimensions must be positive");
    }

    if (cfg.lut.transmittanceR <= 0 || cfg.lut.transmittanceMu <= 0) {
        throw std::runtime_error("Transmittance LUT dimensions must be positive");
    }

    if (cfg.lut.integrationStepsTransmittance <= 0 ||
        cfg.lut.integrationStepsView <= 0 ||
        cfg.lut.integrationStepsSun <= 0) {
        throw std::runtime_error("Integration step counts must be positive");
    }

    if (cfg.lut.irradianceR <= 0 || cfg.lut.irradianceMuS <= 0) {
        throw std::runtime_error("Irradiance LUT dimensions must be positive");
    }

    if (!(cfg.mie.phase.g > -1.0f && cfg.mie.phase.g < 1.0f)) {
        throw std::runtime_error("mie.phase_function.g must be in (-1, 1)");
    }

    if (cfg.multiple.enabled && cfg.multiple.scatteringOrders < 2) {
        throw std::runtime_error(
            "multiple.scatteringOrders must be >= 2 when multiple scattering is enabled");
    }

    validateStrictlyIncreasing(buildWavelengthGrid_nm(cfg.wavelength), "wavelength grid");
}

void validateSpectralLengths(const atmo::AtmosphereMetadata& m) {
    const size_t n = m.wavelengths_nm.size();

    auto check = [n](const std::vector<float>& v, const char* name) {
        if (v.size() != n) {
            throw std::runtime_error(
                std::string("Spectrum length mismatch for ") + name +
                ": expected " + std::to_string(n) +
                ", got " + std::to_string(v.size()));
        }
    };

    check(m.solarIrradiance, "solarIrradiance");
    check(m.rayleighScattering, "rayleighScattering");
    check(m.mieScattering, "mieScattering");
    check(m.mieExtinction, "mieExtinction");
    check(m.absorptionExtinction, "absorptionExtinction");
}

} // namespace

std::vector<float> buildWavelengthGrid_nm(const WavelengthSpec& spec) {
    switch (spec.mode) {
    case WavelengthMode::UniformStep: {
        if (!(spec.maxNm > spec.minNm)) {
            throw std::runtime_error("wavelength.max_nm must be greater than min_nm");
        }
        if (!(spec.stepNm > 0.0f)) {
            throw std::runtime_error("wavelength.step_nm must be positive");
        }

        std::vector<float> out;
        const float eps = 1e-4f;
        for (float lam = spec.minNm; lam <= spec.maxNm + eps; lam += spec.stepNm) {
            out.push_back(lam);
        }
        if (out.empty()) {
            throw std::runtime_error("Generated wavelength grid is empty");
        }
        return out;
    }

    case WavelengthMode::UniformCount: {
        if (!(spec.maxNm > spec.minNm)) {
            throw std::runtime_error("wavelength.max_nm must be greater than min_nm");
        }
        if (spec.count < 2) {
            throw std::runtime_error("wavelength.count must be >= 2");
        }

        std::vector<float> out;
        out.reserve(static_cast<size_t>(spec.count));
        const float step = (spec.maxNm - spec.minNm) / float(spec.count - 1);
        for (int i = 0; i < spec.count; ++i) {
            out.push_back(spec.minNm + step * float(i));
        }
        return out;
    }

    case WavelengthMode::Explicit:
        validateStrictlyIncreasing(spec.explicitValuesNm, "wavelength.values_nm");
        return spec.explicitValuesNm;
    }

    throw std::runtime_error("Unhandled WavelengthMode");
}

void populateMetadataFromConfig(BakeConfig& cfg) {
    auto& m = cfg.meta;

    m.metaDataVersion = 1u;

    m.bottomRadius_m = cfg.planet.bottomRadius_m;
    m.topRadius_m = cfg.planet.topRadius_m;
    m.observerAltitude_m = cfg.planet.observerAltitude_m;
    m.muSMin = cfg.planet.muSMin;

    m.wavelengths_nm = buildWavelengthGrid_nm(cfg.wavelength);

    m.rayleighDensity = toAtmoDensityProfile(cfg.rayleigh.density);
    m.mieDensity      = toAtmoDensityProfile(cfg.mie.density);

    if (cfg.absorption.enabled) {
        m.absorptionDensity = toAtmoDensityProfile(cfg.absorption.density);
    } else {
        m.absorptionDensity = toAtmoDensityProfile(makeZeroProfileConfig());
    }

    if (cfg.solar.mode == DataSourceMode::Builtin) {
        m.solarIrradiance = makeBuiltinSolarIrradiance(m.wavelengths_nm);
    } else {
        if (!cfg.solar.csv.has_value()) {
            throw std::runtime_error("solar_irradiance.mode=file_csv but csv spec is missing");
        }
        const auto csvPath = resolvePath(cfg.configDir, cfg.solar.csv->path);
        m.solarIrradiance = loadCsvSpectrumResampled(
            csvPath,
            cfg.solar.csv->lambdaColumn,
            cfg.solar.csv->valueColumn,
            m.wavelengths_nm);
    }

    if (cfg.rayleigh.scatteringMode == DataSourceMode::Builtin) {
        m.rayleighScattering = makeBuiltinRayleighScattering(m.wavelengths_nm);
    } else {
        if (!cfg.rayleigh.scatteringCsv.has_value()) {
            throw std::runtime_error("rayleigh.scattering.mode=file_csv but csv spec is missing");
        }
        const auto csvPath = resolvePath(cfg.configDir, cfg.rayleigh.scatteringCsv->path);
        m.rayleighScattering = loadCsvSpectrumResampled(
            csvPath,
            cfg.rayleigh.scatteringCsv->lambdaColumn,
            cfg.rayleigh.scatteringCsv->valueColumn,
            m.wavelengths_nm);
    }

    if (cfg.mie.phase.type != "cornette_shanks") {
        throw std::runtime_error(
            "Only mie.phase_function.type = 'cornette_shanks' is supported currently");
    }
    m.miePhaseFunctionG = cfg.mie.phase.g;

    if (cfg.mie.coefficientsMode == DataSourceMode::Builtin) {
        makeBuiltinMieCoefficients(
            m.wavelengths_nm,
            m.mieScattering,
            m.mieExtinction);
    } else {
        if (!cfg.mie.coefficientsCsv.has_value()) {
            throw std::runtime_error("mie.coefficients.mode=file_csv but csv spec is missing");
        }
        const auto csvPath = resolvePath(cfg.configDir, cfg.mie.coefficientsCsv->path);
        loadCsvMieCoefficientsResampled(
            csvPath,
            cfg.mie.coefficientsCsv->lambdaColumn,
            cfg.mie.coefficientsCsv->scatteringColumn,
            cfg.mie.coefficientsCsv->extinctionColumn,
            m.wavelengths_nm,
            m.mieScattering,
            m.mieExtinction);
    }

    if (!cfg.absorption.enabled) {
        m.absorptionExtinction.assign(m.wavelengths_nm.size(), 0.0f);
        m.absorptionDensity = toAtmoDensityProfile(makeZeroProfileConfig());
    } else if (cfg.absorption.extinctionMode == DataSourceMode::Builtin) {
        m.absorptionExtinction =
            makeBuiltinAbsorptionExtinction(m.wavelengths_nm, true);
    } else {
        if (!cfg.absorption.extinctionCsv.has_value()) {
            throw std::runtime_error("absorption.extinction.mode=file_csv but csv spec is missing");
        }
        const auto csvPath = resolvePath(cfg.configDir, cfg.absorption.extinctionCsv->path);
        m.absorptionExtinction = loadCsvSpectrumResampled(
            csvPath,
            cfg.absorption.extinctionCsv->lambdaColumn,
            cfg.absorption.extinctionCsv->valueColumn,
            m.wavelengths_nm);
    }

    m.scatteringOrders = static_cast<std::uint32_t>(
        cfg.multiple.enabled ? cfg.multiple.scatteringOrders : 1);

    validateCommonRanges(cfg);
    validateSpectralLengths(m);
}

BakeConfig loadConfigFromJson(const std::filesystem::path& path) {
    json root;

    try {
        std::ifstream ifs(path);
        if (!ifs) {
            throw std::runtime_error("Failed to open config file: " + path.string());
        }
        root = json::parse(ifs);
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to parse config JSON: " + std::string(e.what()));
    }

    BakeConfig cfg;
    cfg.configPath = std::filesystem::absolute(path);
    cfg.configDir = cfg.configPath.parent_path();

    if (const json* j = findObject(root, "planet")) {
        cfg.planet.bottomRadius_m =
            optionalValue<float>(*j, "bottom_radius_m", cfg.planet.bottomRadius_m);
        cfg.planet.topRadius_m =
            optionalValue<float>(*j, "top_radius_m", cfg.planet.topRadius_m);
        cfg.planet.observerAltitude_m =
            optionalValue<float>(*j, "observer_altitude_m", cfg.planet.observerAltitude_m);
        cfg.planet.muSMin =
            optionalValue<float>(*j, "mu_s_min", cfg.planet.muSMin);
    }

    if (const json* j = findObject(root, "wavelength")) {
        cfg.wavelength.mode =
            parseWavelengthMode(optionalValue<std::string>(*j, "mode", "uniform_step"));

        switch (cfg.wavelength.mode) {
        case WavelengthMode::UniformStep:
            cfg.wavelength.minNm = require<float>(*j, "min_nm");
            cfg.wavelength.maxNm = require<float>(*j, "max_nm");
            cfg.wavelength.stepNm = requirePositiveFloat(*j, "step_nm");
            break;

        case WavelengthMode::UniformCount:
            cfg.wavelength.minNm = require<float>(*j, "min_nm");
            cfg.wavelength.maxNm = require<float>(*j, "max_nm");
            cfg.wavelength.count = requirePositiveInt(*j, "count");
            break;

        case WavelengthMode::Explicit:
            cfg.wavelength.explicitValuesNm = require<std::vector<float>>(*j, "values_nm");
            break;
        }
    }

    if (const json* j = findObject(root, "solar_irradiance")) {
        cfg.solar.mode =
            parseDataSourceMode(optionalValue<std::string>(*j, "mode", "builtin"));
        if (cfg.solar.mode == DataSourceMode::FileCsv) {
            cfg.solar.csv = parseCsvColumnSpec(*j);
        }
    }

    if (const json* j = findObject(root, "rayleigh")) {
        if (const json* d = findObject(*j, "density")) {
            cfg.rayleigh.density =
                parseDensityProfileConfig(*d, makeExponentialProfileConfig(8000.0f));
        } else {
            cfg.rayleigh.density = makeExponentialProfileConfig(8000.0f);
        }

        if (const json* s = findObject(*j, "scattering")) {
            cfg.rayleigh.scatteringMode =
                parseDataSourceMode(optionalValue<std::string>(*s, "mode", "builtin"));
            if (cfg.rayleigh.scatteringMode == DataSourceMode::FileCsv) {
                cfg.rayleigh.scatteringCsv = parseCsvColumnSpec(*s);
            }
        }
    } else {
        cfg.rayleigh.density = makeExponentialProfileConfig(8000.0f);
    }

    if (const json* j = findObject(root, "mie")) {
        if (const json* d = findObject(*j, "density")) {
            cfg.mie.density =
                parseDensityProfileConfig(*d, makeExponentialProfileConfig(1200.0f));
        } else {
            cfg.mie.density = makeExponentialProfileConfig(1200.0f);
        }

        if (const json* p = findObject(*j, "phase_function")) {
            cfg.mie.phase.type =
                optionalValue<std::string>(*p, "type", cfg.mie.phase.type);
            cfg.mie.phase.g =
                optionalValue<float>(*p, "g", cfg.mie.phase.g);
        }

        if (const json* c = findObject(*j, "coefficients")) {
            cfg.mie.coefficientsMode =
                parseDataSourceMode(optionalValue<std::string>(*c, "mode", "builtin"));
            if (cfg.mie.coefficientsMode == DataSourceMode::FileCsv) {
                cfg.mie.coefficientsCsv = parseCsvMieSpec(*c);
            }
        }
    } else {
        cfg.mie.density = makeExponentialProfileConfig(1200.0f);
    }

    if (const json* j = findObject(root, "absorption")) {
        cfg.absorption.enabled =
            optionalValue<bool>(*j, "enabled", cfg.absorption.enabled);

        if (const json* d = findObject(*j, "density")) {
            cfg.absorption.density =
                parseDensityProfileConfig(*d, makeBuiltinOzoneProfileConfig());
        } else {
            cfg.absorption.density = makeBuiltinOzoneProfileConfig();
        }

        if (const json* e = findObject(*j, "extinction")) {
            cfg.absorption.extinctionMode =
                parseDataSourceMode(optionalValue<std::string>(*e, "mode", "builtin"));
            if (cfg.absorption.extinctionMode == DataSourceMode::FileCsv) {
                cfg.absorption.extinctionCsv = parseCsvColumnSpec(*e);
            }
        }
    } else {
        cfg.absorption.enabled = false;
        cfg.absorption.density = makeZeroProfileConfig();
    }

    if (const json* j = findObject(root, "lut")) {
        cfg.lut.nMu = optionalValue<int>(*j, "n_mu", cfg.lut.nMu);
        cfg.lut.nMuS = optionalValue<int>(*j, "n_mu_s", cfg.lut.nMuS);
        cfg.lut.nNu = optionalValue<int>(*j, "n_nu", cfg.lut.nNu);

        cfg.lut.transmittanceR =
            optionalValue<int>(*j, "transmittance_r", cfg.lut.transmittanceR);
        cfg.lut.transmittanceMu =
            optionalValue<int>(*j, "transmittance_mu", cfg.lut.transmittanceMu);

        cfg.lut.integrationStepsTransmittance =
            optionalValue<int>(*j, "integration_steps_transmittance",
                               cfg.lut.integrationStepsTransmittance);
        cfg.lut.integrationStepsView =
            optionalValue<int>(*j, "integration_steps_view", cfg.lut.integrationStepsView);
        cfg.lut.integrationStepsSun =
            optionalValue<int>(*j, "integration_steps_sun", cfg.lut.integrationStepsSun);

        cfg.lut.irradianceR =
            optionalValue<int>(*j, "irradiance_r", cfg.lut.irradianceR);
        cfg.lut.irradianceMuS =
            optionalValue<int>(*j, "irradiance_mu_s", cfg.lut.irradianceMuS);
    }

    if (const json* j = findObject(root, "multiple_scattering")) {
        cfg.multiple.enabled =
            optionalValue<bool>(*j, "enabled", cfg.multiple.enabled);
        cfg.multiple.scatteringOrders =
            optionalValue<int>(*j, "scattering_orders", cfg.multiple.scatteringOrders);
    }

    populateMetadataFromConfig(cfg);
    return cfg;
}

} // namespace bake