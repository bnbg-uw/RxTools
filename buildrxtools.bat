SET vcpkg=C:\vcpkg
SET triplet=x64-windows
SETLOCAL

IF EXIST src\lapisgis\ (
	pushd src\lapisgis\
	git pull
	popd
) ELSE (
	pushd src
	git clone --recurse-submodules https://github.com/jontkane/LapisGis
	popd
)

IF EXIST src\lico\ (
	pushd src\lico\
	git pull
	popd
) ELSE (
	pushd src
	git clone https://github.com/jontkane/lico
	popd
)

IF EXIST src\processedfolder\ (
	pushd src\processedfolder\
	git pull
	popd
) ELSE (
	pushd src
	git clone https://github.com/bnbg-uw/ProcessedFolder
	popd
)

rmdir /S /Q build
mkdir build
pushd build
cmake .. -DCMAKE_TOOLCHAIN_FILE=%vcpkg%\scripts\buildsystems\vcpkg.cmake -DVCPKG_TARGET_TRIPLET=%triplet%
popd

build\rxtools.sln