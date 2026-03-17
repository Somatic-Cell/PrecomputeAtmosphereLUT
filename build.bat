rmdir /S build
mkdir build
cd build
@REM cmake .. -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake
cmake ..
cmake --build . --config Release --verbose 
cd ..