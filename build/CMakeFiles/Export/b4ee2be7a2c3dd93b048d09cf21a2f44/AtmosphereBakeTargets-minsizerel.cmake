#----------------------------------------------------------------
# Generated CMake target import file for configuration "MinSizeRel".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Atmo::AtmoBakeCore" for configuration "MinSizeRel"
set_property(TARGET Atmo::AtmoBakeCore APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(Atmo::AtmoBakeCore PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_MINSIZEREL "CUDA;CXX"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/lib/AtmoBakeCore.lib"
  )

list(APPEND _cmake_import_check_targets Atmo::AtmoBakeCore )
list(APPEND _cmake_import_check_files_for_Atmo::AtmoBakeCore "${_IMPORT_PREFIX}/lib/AtmoBakeCore.lib" )

# Import target "Atmo::AtmoPreviewCore" for configuration "MinSizeRel"
set_property(TARGET Atmo::AtmoPreviewCore APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(Atmo::AtmoPreviewCore PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_MINSIZEREL "CXX"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/lib/AtmoPreviewCore.lib"
  )

list(APPEND _cmake_import_check_targets Atmo::AtmoPreviewCore )
list(APPEND _cmake_import_check_files_for_Atmo::AtmoPreviewCore "${_IMPORT_PREFIX}/lib/AtmoPreviewCore.lib" )

# Import target "Atmo::AtmosphereBake" for configuration "MinSizeRel"
set_property(TARGET Atmo::AtmosphereBake APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(Atmo::AtmosphereBake PROPERTIES
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/bin/AtmosphereBake.exe"
  )

list(APPEND _cmake_import_check_targets Atmo::AtmosphereBake )
list(APPEND _cmake_import_check_files_for_Atmo::AtmosphereBake "${_IMPORT_PREFIX}/bin/AtmosphereBake.exe" )

# Import target "Atmo::AtmosphereBakePreview" for configuration "MinSizeRel"
set_property(TARGET Atmo::AtmosphereBakePreview APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(Atmo::AtmosphereBakePreview PROPERTIES
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/bin/AtmosphereBakePreview.exe"
  )

list(APPEND _cmake_import_check_targets Atmo::AtmosphereBakePreview )
list(APPEND _cmake_import_check_files_for_Atmo::AtmosphereBakePreview "${_IMPORT_PREFIX}/bin/AtmosphereBakePreview.exe" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
