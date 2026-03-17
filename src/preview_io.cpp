#include "preview_common.h"

#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace preview {
namespace {

constexpr std::uint32_t kAtmoMagic = 0x4F4D5441u; // 'ATMO'
constexpr std::uint32_t kSupportedVersion = 1u;

bool readTextFile(const std::filesystem::path& path, std::string& outText, std::string& outError)
{
    std::ifstream ifs(path, std::ios::binary);
    if(!ifs) {
        outError = "Failed to open file: " + path.string();
        return false;
    }
    std::ostringstream oss;
    oss << ifs.rdbuf();
    outText = oss.str();
    return true;
}

std::size_t skipWs(const std::string& s, std::size_t pos)
{
    while(pos < s.size() && std::isspace(static_cast<unsigned char>(s[pos])) != 0) {
        ++pos;
    }
    return pos;
}

bool extractObjectText(const std::string& text, const std::string& key, std::string& outObject)
{
    const std::string needle = "\"" + key + "\"";
    const std::size_t keyPos = text.find(needle);
    if(keyPos == std::string::npos) {
        return false;
    }

    std::size_t pos = text.find(':', keyPos + needle.size());
    if(pos == std::string::npos) {
        return false;
    }
    pos = skipWs(text, pos + 1);
    if(pos >= text.size() || text[pos] != '{') {
        return false;
    }

    int depth = 0;
    const std::size_t begin = pos;
    for(; pos < text.size(); ++pos) {
        if(text[pos] == '{') {
            ++depth;
        }
        else if(text[pos] == '}') {
            --depth;
            if(depth == 0) {
                outObject.assign(text.begin() + static_cast<std::ptrdiff_t>(begin),
                                 text.begin() + static_cast<std::ptrdiff_t>(pos + 1));
                return true;
            }
        }
    }
    return false;
}

bool extractArrayText(const std::string& text, const std::string& key, std::string& outArray)
{
    const std::string needle = "\"" + key + "\"";
    const std::size_t keyPos = text.find(needle);
    if(keyPos == std::string::npos) {
        return false;
    }

    std::size_t pos = text.find(':', keyPos + needle.size());
    if(pos == std::string::npos) {
        return false;
    }
    pos = skipWs(text, pos + 1);
    if(pos >= text.size() || text[pos] != '[') {
        return false;
    }

    int depth = 0;
    const std::size_t begin = pos;
    for(; pos < text.size(); ++pos) {
        if(text[pos] == '[') {
            ++depth;
        }
        else if(text[pos] == ']') {
            --depth;
            if(depth == 0) {
                outArray.assign(text.begin() + static_cast<std::ptrdiff_t>(begin),
                                text.begin() + static_cast<std::ptrdiff_t>(pos + 1));
                return true;
            }
        }
    }
    return false;
}

bool parseLeadingFloat(const std::string& text, std::size_t pos, float& outValue)
{
    pos = skipWs(text, pos);
    if(pos >= text.size()) {
        return false;
    }

    const char* begin = text.c_str() + pos;
    char* end = nullptr;
    const float value = std::strtof(begin, &end);
    if(end == begin) {
        return false;
    }
    if(value == HUGE_VALF || value == -HUGE_VALF || !std::isfinite(value)) {
        return false;
    }
    outValue = value;
    return true;
}

bool extractNumberByKeys(const std::string& text,
                         const std::initializer_list<const char*>& keys,
                         float& outValue)
{
    for(const char* key : keys) {
        const std::string needle = std::string("\"") + key + "\"";
        const std::size_t keyPos = text.find(needle);
        if(keyPos == std::string::npos) {
            continue;
        }
        const std::size_t colonPos = text.find(':', keyPos + needle.size());
        if(colonPos == std::string::npos) {
            continue;
        }
        if(parseLeadingFloat(text, colonPos + 1, outValue)) {
            return true;
        }
    }
    return false;
}

bool extractNumberByKeys(const std::string& text,
                         const std::initializer_list<const char*>& keys,
                         std::uint32_t& outValue)
{
    float tmp = 0.0f;
    if(!extractNumberByKeys(text, keys, tmp)) {
        return false;
    }
    if(tmp < 0.0f || tmp > static_cast<float>(std::numeric_limits<std::uint32_t>::max())) {
        return false;
    }
    outValue = static_cast<std::uint32_t>(tmp);
    return true;
}

bool parseFloatArray(const std::string& arrayText, std::vector<float>& outValues)
{
    outValues.clear();
    if(arrayText.size() < 2 || arrayText.front() != '[' || arrayText.back() != ']') {
        return false;
    }

    std::size_t pos = 1;
    while(pos + 1 < arrayText.size()) {
        pos = skipWs(arrayText, pos);
        if(pos >= arrayText.size() - 1) {
            break;
        }
        if(arrayText[pos] == ',') {
            ++pos;
            continue;
        }

        float value = 0.0f;
        if(!parseLeadingFloat(arrayText, pos, value)) {
            return false;
        }

        const char* begin = arrayText.c_str() + pos;
        char* end = nullptr;
        std::strtof(begin, &end);
        pos += static_cast<std::size_t>(end - begin);
        outValues.push_back(value);
    }
    return !outValues.empty();
}

std::uint64_t computeElementCount(const atmo::TableHeader& header)
{
    return static_cast<std::uint64_t>(header.dim0) *
           static_cast<std::uint64_t>(header.dim1) *
           static_cast<std::uint64_t>(header.dim2) *
           static_cast<std::uint64_t>(header.dim3);
}

bool readTableWithHeader(const std::filesystem::path& path,
                         atmo::TableType expectedType,
                         FinalTable& outTable,
                         std::string& outError)
{
    std::ifstream ifs(path, std::ios::binary | std::ios::ate);
    if(!ifs) {
        outError = "Failed to open LUT file: " + path.string();
        return false;
    }

    const std::streamsize fileSize = ifs.tellg();
    if(fileSize < static_cast<std::streamsize>(sizeof(atmo::TableHeader))) {
        outError = "LUT file is too small to contain TableHeader: " + path.string();
        return false;
    }

    ifs.seekg(0, std::ios::beg);
    if(!ifs.read(reinterpret_cast<char*>(&outTable.header), sizeof(outTable.header))) {
        outError = "Failed to read TableHeader: " + path.string();
        return false;
    }

    if(outTable.header.magic != kAtmoMagic) {
        outError = "Invalid table magic in: " + path.string();
        return false;
    }
    if(outTable.header.version != kSupportedVersion) {
        outError = "Unsupported table version in: " + path.string();
        return false;
    }
    if(outTable.header.tableType != expectedType) {
        outError = "Unexpected table type in: " + path.string();
        return false;
    }
    if(outTable.header.storageFormat != atmo::StorageFormat::FP32) {
        outError = "Unsupported storage format in: " + path.string();
        return false;
    }

    const std::uint64_t count = computeElementCount(outTable.header);
    if(count == 0) {
        outError = "Table has zero-sized dimension(s): " + path.string();
        return false;
    }

    const std::uint64_t expectedPayloadBytes = count * static_cast<std::uint64_t>(sizeof(float));
    const std::uint64_t actualPayloadBytes =
        static_cast<std::uint64_t>(fileSize) - static_cast<std::uint64_t>(sizeof(atmo::TableHeader));

    if(actualPayloadBytes != expectedPayloadBytes) {
        std::ostringstream oss;
        oss << "Payload size mismatch in " << path.string()
            << ": header expects " << expectedPayloadBytes
            << " bytes, but file contains " << actualPayloadBytes << " bytes.";
        outError = oss.str();
        return false;
    }

    outTable.values.resize(static_cast<std::size_t>(count));
    if(!ifs.read(reinterpret_cast<char*>(outTable.values.data()),
                 static_cast<std::streamsize>(expectedPayloadBytes))) {
        outError = "Failed to read LUT payload: " + path.string();
        return false;
    }

    return true;
}

bool validateSkyHeaderCompatibility(const atmo::TableHeader& ref,
                                    const atmo::TableHeader& other,
                                    const char* otherName,
                                    std::string& outError)
{
    if(ref.dim0 != other.dim0 || ref.dim1 != other.dim1 || ref.dim2 != other.dim2 || ref.dim3 != other.dim3) {
        std::ostringstream oss;
        oss << "Sky LUT dimension mismatch between sky_rayleigh_single.bin and "
            << otherName << ": ("
            << ref.dim0 << ", " << ref.dim1 << ", " << ref.dim2 << ", " << ref.dim3 << ") vs ("
            << other.dim0 << ", " << other.dim1 << ", " << other.dim2 << ", " << other.dim3 << ").";
        outError = oss.str();
        return false;
    }
    return true;
}

bool validateOrAssignDim(const char* label,
                         std::uint32_t currentValue,
                         std::uint32_t headerValue,
                         std::uint32_t& outValue,
                         std::string& outError)
{
    if(currentValue != 0 && currentValue != headerValue) {
        std::ostringstream oss;
        oss << "Dimension mismatch for " << label << ": metadata/CLI says "
            << currentValue << ", but LUT header says " << headerValue << ".";
        outError = oss.str();
        return false;
    }
    outValue = headerValue;
    return true;
}

} // namespace

bool loadMetadataJson(const std::filesystem::path& path,
                      PreviewMetadata& outMeta,
                      std::string& outError)
{
    std::string text;
    if(!readTextFile(path, text, outError)) {
        return false;
    }

    std::string planetText;
    if(extractObjectText(text, "planet", planetText)) {
        extractNumberByKeys(planetText, {"bottom_radius_m", "bottomRadiusM"}, outMeta.bottomRadiusM);
        extractNumberByKeys(planetText, {"top_radius_m", "topRadiusM"}, outMeta.topRadiusM);
        extractNumberByKeys(planetText, {"observer_altitude_m", "observerAltitudeM"}, outMeta.observerAltitudeM);
        extractNumberByKeys(planetText, {"mu_s_min", "muSMin"}, outMeta.muSMin);
    }

    std::string mieText;
    if(extractObjectText(text, "mie_phase_function", mieText)) {
        extractNumberByKeys(mieText, {"g", "miePhaseFunctionG"}, outMeta.miePhaseFunctionG);
    }
    else {
        extractNumberByKeys(text, {"mie_phase_function_g", "miePhaseFunctionG"}, outMeta.miePhaseFunctionG);
    }

    std::string lutText;
    if(extractObjectText(text, "lut", lutText)) {
        extractNumberByKeys(lutText, {"sky_mu", "n_mu", "skyMu"}, outMeta.skyMu);
        extractNumberByKeys(lutText, {"sky_mu_s", "n_mu_s", "skyMuS"}, outMeta.skyMuS);
        extractNumberByKeys(lutText, {"sky_nu", "n_nu", "skyNu"}, outMeta.skyNu);

        extractNumberByKeys(lutText, {"irradiance_r", "n_irradiance_r", "irradianceR"}, outMeta.irradianceR);
        extractNumberByKeys(lutText, {"irradiance_mu_s", "n_irradiance_mu_s", "irradianceMuS"}, outMeta.irradianceMuS);
    }

    std::string wavelengthsText;
    if(!extractArrayText(text, "wavelengths_nm", wavelengthsText)) {
        outError = "metadata.json does not contain wavelengths_nm array.";
        return false;
    }
    if(!parseFloatArray(wavelengthsText, outMeta.wavelengthsNm)) {
        outError = "Failed to parse wavelengths_nm array in metadata.json.";
        return false;
    }

    return true;
}

bool loadFinalLuts(const std::filesystem::path& outDir,
                   PreviewMetadata& inOutMeta,
                   FinalLuts& outLuts,
                   std::string& outError)
{
    if(inOutMeta.wavelengthsNm.empty()) {
        outError = "wavelengths_nm is empty. metadata.json is required for spectral preview.";
        return false;
    }

    if(!readTableWithHeader(outDir / "transmittance_internal.bin",
                            atmo::TableType::SunTransmittance,
                            outLuts.sunTransmittance,
                            outError)) {
        return false;
    }
    if(!readTableWithHeader(outDir / "sky_rayleigh_single.bin",
                            atmo::TableType::SkyRayleighSingle,
                            outLuts.skyRayleighSingle,
                            outError)) {
        return false;
    }
    if(!readTableWithHeader(outDir / "sky_mie_single.bin",
                            atmo::TableType::SkyMieSingle,
                            outLuts.skyMieSingle,
                            outError)) {
        return false;
    }
    if(!readTableWithHeader(outDir / "sky_multiple.bin",
                            atmo::TableType::SkyMultiple,
                            outLuts.skyMultiple,
                            outError)) {
        return false;
    }
    if(!readTableWithHeader(outDir / "direct_irradiance.bin",
                            atmo::TableType::DirectIrradiance,
                            outLuts.directIrradiance,
                            outError)) {
        return false;
    }

    const atmo::TableHeader& skyHdr = outLuts.skyRayleighSingle.header;
    if(!validateSkyHeaderCompatibility(skyHdr, outLuts.skyMieSingle.header, "sky_mie_single.bin", outError)) {
        return false;
    }
    if(!validateSkyHeaderCompatibility(skyHdr, outLuts.skyMultiple.header, "sky_multiple.bin", outError)) {
        return false;
    }

    if(skyHdr.dim0 == 0 || skyHdr.dim1 == 0 || skyHdr.dim2 == 0 || skyHdr.dim3 == 0) {
        outError = "Sky LUT header contains zero dimension(s).";
        return false;
    }

    std::uint32_t resolvedSkyNu = 0;
    std::uint32_t resolvedSkyMu = 0;
    std::uint32_t resolvedSkyMuS = 0;

    if(!validateOrAssignDim("skyNu", inOutMeta.skyNu, skyHdr.dim0, resolvedSkyNu, outError)) {
        return false;
    }
    if(!validateOrAssignDim("skyMu", inOutMeta.skyMu, skyHdr.dim1, resolvedSkyMu, outError)) {
        return false;
    }
    if(!validateOrAssignDim("skyMuS", inOutMeta.skyMuS, skyHdr.dim2, resolvedSkyMuS, outError)) {
        return false;
    }

    const std::uint32_t lambdaCount = skyHdr.dim3;
    if(inOutMeta.wavelengthsNm.size() != static_cast<std::size_t>(lambdaCount)) {
        std::ostringstream oss;
        oss << "wavelength count mismatch: metadata.json has "
            << inOutMeta.wavelengthsNm.size() << " wavelengths, but sky LUT header says "
            << lambdaCount << ".";
        outError = oss.str();
        return false;
    }

    const atmo::TableHeader& sunHdr = outLuts.sunTransmittance.header;
    if(sunHdr.dim0 != resolvedSkyMuS) {
        std::ostringstream oss;
        oss << "sun_transmittance.bin dim0 (mu_s count) is " << sunHdr.dim0
            << ", but sky LUT dim2 (mu_s count) is " << resolvedSkyMuS << ".";
        outError = oss.str();
        return false;
    }
    if(sunHdr.dim1 != lambdaCount) {
        std::ostringstream oss;
        oss << "sun_transmittance.bin dim1 (lambda count) is " << sunHdr.dim1
            << ", but sky LUT dim3 (lambda count) is " << lambdaCount << ".";
        outError = oss.str();
        return false;
    }

    const atmo::TableHeader& irrHdr = outLuts.directIrradiance.header;
    if(irrHdr.dim0 == 0 || irrHdr.dim1 == 0 || irrHdr.dim2 == 0 || irrHdr.dim3 == 0) {
        outError = "direct_irradiance.bin header contains zero dimension(s).";
        return false;
    }

    std::uint32_t resolvedIrradianceR = 0;
    std::uint32_t resolvedIrradianceMuS = 0;

    if(!validateOrAssignDim("irradianceR", inOutMeta.irradianceR, irrHdr.dim0, resolvedIrradianceR, outError)) {
        return false;
    }
    if(!validateOrAssignDim("irradianceMuS", inOutMeta.irradianceMuS, irrHdr.dim1, resolvedIrradianceMuS, outError)) {
        return false;
    }

    if(irrHdr.dim2 != lambdaCount) {
        std::ostringstream oss;
        oss << "direct_irradiance.bin dim2 (lambda count) is " << irrHdr.dim2
            << ", but sky LUT dim3 (lambda count) is " << lambdaCount << ".";
        outError = oss.str();
        return false;
    }

    if(irrHdr.dim3 != 1u) {
        std::ostringstream oss;
        oss << "direct_irradiance.bin dim3 is " << irrHdr.dim3
            << ", expected 1.";
        outError = oss.str();
        return false;
    }

    inOutMeta.skyNu = resolvedSkyNu;
    inOutMeta.skyMu = resolvedSkyMu;
    inOutMeta.skyMuS = resolvedSkyMuS;
    inOutMeta.irradianceR = resolvedIrradianceR;
    inOutMeta.irradianceMuS = resolvedIrradianceMuS;
    return true;
}

} // namespace preview
