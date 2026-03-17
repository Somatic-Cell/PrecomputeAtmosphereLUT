cmake_minimum_required(VERSION 3.24)
project(AtmosphereBake VERSION 0.1.0 LANGUAGES CXX CUDA)

find_package(CUDAToolkit REQUIRED)

add_subdirectory(${PROJECT_SOURCE_DIR}/ext/json)
add_subdirectory(${PROJECT_SOURCE_DIR}/ext/fpng)


# ----------------------------
# 共有ヘッダだけを公開する target
# ----------------------------
add_library(AtmoHeaders INTERFACE)

target_include_directories(AtmoHeaders
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(AtmoHeaders INTERFACE cxx_std_20)

# ----------------------------
# 前計算の実装本体
# ----------------------------
add_library(AtmoBakeCore STATIC
  src/config_io.cpp
  src/csv_table_io.cpp
  src/serializer.cpp
  src/bake_transmittance.cu
  src/bake_sun_transmittance.cu
  src/bake_sky_single_scatter.cu
  src/bake_irradiance.cu
  src/bake_scattering.cu
)

target_include_directories(AtmoBakeCore
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(AtmoBakeCore
  PUBLIC
    AtmoHeaders
  PRIVATE
    CUDA::cudart
    nlohmann_json::nlohmann_json
)

target_compile_features(AtmoBakeCore PUBLIC cxx_std_20)

set_target_properties(AtmoBakeCore PROPERTIES
  CUDA_STANDARD 17
  CUDA_STANDARD_REQUIRED ON
  CXX_STANDARD 20
  CXX_STANDARD_REQUIRED ON
  CUDA_ARCHITECTURES native
)

# ----------------------------
# Preview 実装本体
# ----------------------------
add_library(AtmoPreviewCore STATIC
  src/preview_io.cpp
  src/preview_render.cpp
  src/cie1931_2deg_1nm.cpp

)

target_link_libraries(AtmoPreviewCore
  PUBLIC
    AtmoHeaders
)

target_compile_features(AtmoPreviewCore PUBLIC cxx_std_20)

set_target_properties(AtmoPreviewCore PROPERTIES
  CXX_STANDARD 20
  CXX_STANDARD_REQUIRED ON
)


# ----------------------------
# CLI 実行ファイル
# ----------------------------
add_executable(AtmosphereBake
  app/main.cpp
)

target_link_libraries(AtmosphereBake
  PRIVATE
    AtmoBakeCore
)

target_compile_features(AtmosphereBake PRIVATE cxx_std_20)

add_executable(AtmosphereBakePreview
  "app/preview.cpp"
  "${PROJECT_SOURCE_DIR}/ext/fpng/src/fpng.h"
  "${PROJECT_SOURCE_DIR}/ext/fpng/src/fpng.cpp"
)

target_link_libraries(AtmosphereBakePreview
  PRIVATE
    AtmoPreviewCore
)

target_compile_features(AtmosphereBakePreview PRIVATE cxx_std_20)

# ----------------------------
# install / export
# ----------------------------
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(TARGETS
  AtmoHeaders
  AtmoBakeCore
  AtmoPreviewCore
  AtmosphereBake
  AtmosphereBakePreview
  EXPORT AtmosphereBakeTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT AtmosphereBakeTargets
  FILE AtmosphereBakeTargets.cmake
  NAMESPACE Atmo::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AtmosphereBake
)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/AtmosphereBakeConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/AtmosphereBakeConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AtmosphereBake
)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/AtmosphereBakeConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/AtmosphereBakeConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/AtmosphereBakeConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AtmosphereBake
)