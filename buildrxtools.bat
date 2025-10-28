SET vcpkg=C:\vcpkg
SET triplet=x64-windows
SET maindir=F:\RxTools
SETLOCAL

rmdir /S /Q %~1\build
mkdir %~1\build
pushd %~1\build
cmake %maindir%\%~1 -DCMAKE_TOOLCHAIN_FILE=%vcpkg%\scripts\buildsystems\vcpkg.cmake -DVCPKG_TARGET_TRIPLET=%triplet% -DPROJ_DB_PATH=%vcpkg%\installed\%triplet%\share\proj\proj.db -DLAPISGISCMAKE_PATH=%maindir%\LapisGis\LapisGis.cmake
popd