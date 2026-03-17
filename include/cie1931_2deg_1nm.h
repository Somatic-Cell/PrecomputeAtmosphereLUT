#pragma once

#include <cstddef>

namespace atmo::preview {

struct CieXyzBar
{
    float xBar;
    float yBar;
    float zBar;
};

// Source:
//   CIE 2019, "Colour-matching functions of CIE 1931 standard colorimetric observer"
//   DOI: 10.25039/CIE.DS.xvudnb9b
// Dataset:
//   CIE_xyz_1931_2deg.csv
// Wavelength range:
//   360 nm to 830 nm at 1 nm increments
// Notes:
//   - Linear interpolation inside the tabulated range
//   - Zero extrapolation outside the tabulated range
constexpr std::size_t kCie1931_2deg_1nm_Count = 471;

float cie1931_2deg_1nm_wavelength_nm(std::size_t index);

// Linear interpolation, zero outside [360, 830].
CieXyzBar sampleCie1931_2deg_1nm(float wavelengthNm);

} // namespace atmo::preview
