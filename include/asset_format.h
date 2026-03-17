#pragma once
#include <cstdint>

namespace atmo {

enum class TableType : uint32_t {
  Transmittance2D = 0,    // 2D の透過率用 LUT (半径方向r と 方向余弦 mu を想定)
  SkyRayleighSingle = 1,  // Rayleigh 単一散乱の位相関数を含めない LUT
  SkyMieSingle = 2,       // Mie 単一散乱の位相関数を含めない LUT
  SkyMultiple = 3,        // 多重散乱の LUT
  SunTransmittance = 4,   // 固定観測点から見た太陽方向の透過率 LUT
  DirectIrradiance = 5    // 
};

enum class StorageFormat : uint32_t {
  FP32 = 0
};

struct TableHeader {
  uint32_t magic = 0x4F4D5441u;   // 'ATMO'
  uint32_t version = 1u;          // ファイルフォーマットのバージョン
  TableType tableType = TableType::SkyRayleighSingle;      // ファイルの中身
  StorageFormat storageFormat = StorageFormat::FP32;
  
  // Dimension semantics depend on table_type.
  // Transmittance2D                        : dim0 = n_mu, dim1 = n_r
  // SkyRayleighSingle/MieSingle /Multiple  : dim0 = n_nu, dim1 = n_mu, dim2 = n_mu_s, dim3 = n_lambda
  // SunTransmittance                       : dim0 = n_mu_s, dim1 = n_lambda
  uint32_t dim0 = 0u;
  uint32_t dim1 = 0u;
  uint32_t dim2 = 0u;
  uint32_t dim3 = 0u;
};

}  // namespace atmo
