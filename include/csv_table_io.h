#pragma once

#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>
#include <cstddef>

namespace bake {

struct CsvTable {
    std::vector<std::string> headers;
    std::unordered_map<std::string, size_t> headerToIndex;
    std::vector<std::vector<std::string>> rows;

    bool hasColumn(const std::string& name) const;
    
    // 列が存在しない場合は例外を投げる
    std::size_t columnIndex(const std::string& name) const;
};

std::filesystem::path resolvePath(
    const std::filesystem::path& baseDir,
    const std::string& pathLike);

// UTF-8 BOM を許容し，ヘッダの重複や不正な行は例外
CsvTable loadCsvTable(const std::filesystem::path& path);

// 指定列を float として読みだす．変換不能な値や列の欠落は例外
std::vector<float> extractNumericColumn(
    const CsvTable& table,
    const std::string& columnName,
    const std::filesystem::path& sourcePath);

// srcLambdaNm と srcValues は同じ長さで、srcLambdaNm は昇順を仮定。
// 範囲外の dstLambdaNm をどう扱うかは実装コメントで明確化する。
std::vector<float> resampleLinearSpectrum(
    const std::vector<float>& srcLambdaNm,
    const std::vector<float>& srcValues,
    const std::vector<float>& dstLambdaNm,
    const std::string& debugName);

std::vector<float> loadCsvSpectrumResampled(
    const std::filesystem::path& path,
    const std::string& lambdaColumn,
    const std::string& valueColumn,
    const std::vector<float>& dstLambdaNm);

void loadCsvMieCoefficientsResampled(
    const std::filesystem::path& path,
    const std::string& lambdaColumn,
    const std::string& scatteringColumn,
    const std::string& extinctionColumn,
    const std::vector<float>& dstLambdaNm,
    std::vector<float>& outScattering,
    std::vector<float>& outExtinction);

} // namespace bake