SET vcpkg=C:\vcpkg
SET triplet=x64-windows
SETLOCAL

rmdir /S /Q build
mkdir build
pushd build
cmake .. -DCMAKE_TOOLCHAIN_FILE=%vcpkg%\scripts\buildsystems\vcpkg.cmake -DVCPKG_TARGET_TRIPLET=%triplet%
popd