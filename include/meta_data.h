#pragma once
#include <vector>
#include <cstdint>

namespace atmo {

// 大気の密度モデル (2017 年の改良版)
struct DensityLayer {
    float width_m = 0.0f;       // この層が有効な高度の区間の幅
    float expTerm = 0.0f;       // 指数項の係数
    float expScale = 0.0f;      // 指数項のスケール
    float linearTerm = 0.0f;    // 線形項のスケール
    float constantTerm = 0.0f;  // 定数項
};

// 2区間の密度表現 (Rayleigh, Mie では 1 区間， Ozone では 2区間を使う)
struct DensityProfile {
    DensityLayer layers[2];
};

// アセット全体の物理パラメータの記録
struct AtmosphereMetadata {
    // metadata.json のバージョン
    uint32_t metaDataVersion = 1u;

    // 惑星，大気の形状
    float bottomRadius_m = 6360000.0f;  // 地球の半径
    float topRadius_m = 6460000.0f;     // 大気の上限
    float observerAltitude_m = 0.0f;    // 観測者の地面からの高度

    // 太陽の天頂方向の余弦が取れる最小値
    float muSMin = -0.2f;

    // 波長方向のチャネル
    std::vector<float> wavelengths_nm;

    // 高度 h [m] における大気の密度
    DensityProfile rayleighDensity;
    DensityProfile mieDensity;
    DensityProfile absorptionDensity;

    // 波長ごとのプロファイル
    std::vector<float> solarIrradiance;       // 大気上端での太陽光源のスペクトル分布
    std::vector<float> rayleighScattering;    // 波長毎のRayleigh 散乱係数
    std::vector<float> mieScattering;         // 波長毎のMie 散乱係数
    std::vector<float> mieExtinction;         // 波長毎のMie 消散係数
    std::vector<float> absorptionExtinction;  // 波長毎の吸収係数

    // 位相関数のパラメータ
    float miePhaseFunctionG = 0.8f;           // Cornette-Shanks 関数のパラメータ

    // 多重散乱で何次散乱まで考慮するか
    uint32_t scatteringOrders = 7u;
};

} // namespace atmo