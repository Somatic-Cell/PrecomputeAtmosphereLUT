#include "csv_table_io.h"

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace bake {

namespace {

std::string trim(std::string s) {
    auto notSpace = [](unsigned char c) { return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), notSpace));
    s.erase(std::find_if(s.rbegin(), s.rend(), notSpace).base(), s.end());
    return s;
}

void stripUtf8Bom(std::string& s) {
    if (s.size() >= 3 &&
        static_cast<unsigned char>(s[0]) == 0xEF &&
        static_cast<unsigned char>(s[1]) == 0xBB &&
        static_cast<unsigned char>(s[2]) == 0xBF) {
        s.erase(0, 3);
    }
}

bool isCommentOrEmpty(const std::string& line) {
    const std::string t = trim(line);
    return t.empty() || (!t.empty() && t[0] == '#');
}

std::vector<std::string> splitCsvLine(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    bool inQuotes = false;

    for (size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];

        if (c == '"') {
            if (inQuotes && i + 1 < line.size() && line[i + 1] == '"') {
                cur.push_back('"');
                ++i;
            } else {
                inQuotes = !inQuotes;
            }
        } else if (c == ',' && !inQuotes) {
            out.push_back(trim(cur));
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    
    if (inQuotes) {
        throw std::runtime_error("Unterminated quoted field in CSV line");
    }

    out.push_back(trim(cur));
    return out;
}

float parseFloatStrict(const std::string& s, const std::string& what) {
    try {
        size_t pos = 0;
        float v = std::stof(s, &pos);
        if (pos != s.size()) {
            throw std::runtime_error("trailing characters");
        }
        if (!std::isfinite(v)){
            throw std::runtime_error("non-finite value");
        }
        return v;
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to parse float for " + what + ": '" + s + "' (" + e.what() + ")");
    }
}

void validateStrictlyIncreasing(
    const std::vector<float>& xs,
    const std::string& name) {
    if (xs.empty()) {
        throw std::runtime_error(name + " must not be empty");
    }
    for (size_t i = 1; i < xs.size(); ++i) {
        if (!(xs[i] > xs[i - 1])) {
            throw std::runtime_error(name + " must be strictly increasing");
        }
    }
}

} // namespace

bool CsvTable::hasColumn(const std::string& name) const {
    return headerToIndex.find(name) != headerToIndex.end();
}

size_t CsvTable::columnIndex(const std::string& name) const {
    auto it = headerToIndex.find(name);
    if (it == headerToIndex.end()) {
        throw std::runtime_error("CSV column not found: " + name);
    }
    return it->second;
}

std::filesystem::path resolvePath(
    const std::filesystem::path& baseDir,
    const std::string& pathLike) {
    std::filesystem::path p(pathLike);
    if (p.is_absolute()) {
        return p;
    }
    return baseDir / p;
}

CsvTable loadCsvTable(const std::filesystem::path& path) {
    std::ifstream ifs(path);
    if (!ifs) {
        throw std::runtime_error("Failed to open CSV file: " + path.string());
    }

    CsvTable table;
    std::string line;
    bool headerRead = false;
    size_t lineNumber = 0;

    while (std::getline(ifs, line)) {
        ++lineNumber;

        if (!headerRead) {
            if (isCommentOrEmpty(line)) {
                continue;
            }

            stripUtf8Bom(line);
            table.headers = splitCsvLine(line);
            if (table.headers.empty()) {
                throw std::runtime_error("CSV header is empty: " + path.string());
            }

            for (size_t i = 0; i < table.headers.size(); ++i) {
                const std::string name = trim(table.headers[i]);
                if (name.empty()) {
                    throw std::runtime_error("CSV header contains empty column name: " + 
                        std::to_string(lineNumber) + ": " + path.string());
                }
                if (table.headerToIndex.count(name) != 0) {
                    throw std::runtime_error("Duplicate CSV column name '" + name +
                        "' at line " + std::to_string(lineNumber) + 
                        "' in " + path.string());
                }
                table.headers[i] = name;
                table.headerToIndex[name] = i;
            }

            headerRead = true;
            continue;
        }

        if (isCommentOrEmpty(line)) {
            continue;
        }

        auto fields = splitCsvLine(line);
        if (fields.size() != table.headers.size()) {
            throw std::runtime_error(
                "CSV row has wrong number of columns at line " +
                std::to_string(lineNumber) + " in " + path.string() +
                " (expected " + std::to_string(table.headers.size()) +
                ", got " + std::to_string(fields.size()) + ")");
        }

        table.rows.push_back(std::move(fields));
    }

    if (!headerRead) {
        throw std::runtime_error("CSV file does not contain a header: " + path.string());
    }

    return table;
}

std::vector<float> extractNumericColumn(
    const CsvTable& table,
    const std::string& columnName,
    const std::filesystem::path& sourcePath) {
    const size_t idx = table.columnIndex(columnName);

    std::vector<float> out;
    out.reserve(table.rows.size());

    for (size_t row = 0; row < table.rows.size(); ++row) {
        const std::string& s = table.rows[row][idx];
        out.push_back(parseFloatStrict(
            s,
            sourcePath.string() + ":" + columnName + " row " + std::to_string(row + 2)));
    }

    return out;
}

std::vector<float> resampleLinearSpectrum(
    const std::vector<float>& srcLambdaNm,
    const std::vector<float>& srcValues,
    const std::vector<float>& dstLambdaNm,
    const std::string& debugName) {
    if (srcLambdaNm.size() != srcValues.size()) {
        throw std::runtime_error(debugName + ": wavelength/value size mismatch");
    }

    validateStrictlyIncreasing(srcLambdaNm, debugName + " source wavelengths");
    validateStrictlyIncreasing(dstLambdaNm, debugName + " destination wavelengths");

    const float eps = 1e-4f;
    const float srcMin = srcLambdaNm.front();
    const float srcMax = srcLambdaNm.back();

    std::vector<float> out;
    out.reserve(dstLambdaNm.size());

    for (float lam : dstLambdaNm) {
        if (lam < srcMin - eps || lam > srcMax + eps) {
            throw std::runtime_error(
                debugName + ": destination wavelength out of source range: " +
                std::to_string(lam) + " not in [" +
                std::to_string(srcMin) + ", " +
                std::to_string(srcMax) + "]");
        }

        if (lam <= srcMin + eps) {
            out.push_back(srcValues.front());
            continue;
        }

        if (lam >= srcMax - eps) {
            out.push_back(srcValues.back());
            continue;
        }

        auto it = std::upper_bound(srcLambdaNm.begin(), srcLambdaNm.end(), lam);
        const size_t i1 = static_cast<size_t>(it - srcLambdaNm.begin());
        const size_t i0 = i1 - 1;

        const float x0 = srcLambdaNm[i0];
        const float x1 = srcLambdaNm[i1];
        const float y0 = srcValues[i0];
        const float y1 = srcValues[i1];

        const float t = (lam - x0) / (x1 - x0);
        out.push_back((1.0f - t) * y0 + t * y1);
    }

    return out;
}

std::vector<float> loadCsvSpectrumResampled(
    const std::filesystem::path& path,
    const std::string& lambdaColumn,
    const std::string& valueColumn,
    const std::vector<float>& dstLambdaNm) {
    const CsvTable table = loadCsvTable(path);
    const auto srcLambda = extractNumericColumn(table, lambdaColumn, path);
    const auto srcValues = extractNumericColumn(table, valueColumn, path);

    return resampleLinearSpectrum(
        srcLambda,
        srcValues,
        dstLambdaNm,
        path.string() + " [" + valueColumn + "]");
}

void loadCsvMieCoefficientsResampled(
    const std::filesystem::path& path,
    const std::string& lambdaColumn,
    const std::string& scatteringColumn,
    const std::string& extinctionColumn,
    const std::vector<float>& dstLambdaNm,
    std::vector<float>& outScattering,
    std::vector<float>& outExtinction) {
    const CsvTable table = loadCsvTable(path);

    const auto srcLambda = extractNumericColumn(table, lambdaColumn, path);
    const auto srcScattering = extractNumericColumn(table, scatteringColumn, path);
    const auto srcExtinction = extractNumericColumn(table, extinctionColumn, path);

    auto scattering = resampleLinearSpectrum(
        srcLambda, srcScattering, dstLambdaNm,
        path.string() + " [" + scatteringColumn + "]");

    auto extinction = resampleLinearSpectrum(
        srcLambda, srcExtinction, dstLambdaNm,
        path.string() + " [" + extinctionColumn + "]");

    outScattering = std::move(scattering);
    outExtinction = std::move(extinction);
}

} // namespace bake