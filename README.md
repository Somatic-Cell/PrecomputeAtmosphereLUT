 AtmosphereBake

`AtmosphereBake` は，固定観測高度・sky-only・spectral な大気散乱 LUT を
オフラインで生成するための CUDA ベースの前計算ツールです。

現時点では主に以下を対象とします。

- 観測者高度固定
- 天球（sky dome / sky sphere）のみ
- spectral な Rayleigh / Mie / absorption
- single scattering を中心とした前計算
- レンダラ起動時ではなく，事前に LUT を生成して保存する運用

## 生成する主な出力

- `transmittance.bin`
  - 大気内の透過率 LUT
- `sun_transmittance.bin`
  - 観測点から見た太陽方向の透過率 LUT
- `sky_rayleigh.bin`
  - Rayleigh single scattering の phase-free LUT
- `sky_mie.bin`
  - Mie single scattering の phase-free LUT
- `metadata.json`
  - LUT の解釈に必要なパラメータ群
- `preview.ppm`（任意）
  - 魚眼形式の全天プレビュー画像

## ディレクトリ構成

```text
AtmosphereBake/
  CMakeLists.txt
  README.md

  include/
    asset_format.h
    metadata.h
    sampling.h
    atmosphere_math.h

  src/
    config_io.h
    config_io.cpp
    csv_table_io.h
    csv_table_io.cpp
    serializer.cpp

    bake_common.cuh
    bake_transmittance.cu
    bake_sun_transmittance.cu
    bake_sky_single_scatter.cu

  app/
    main.cpp
    preview.cpp

  data/
    solar_am0.csv
    mie_coefficients.csv
```
