include_guard(GLOBAL)

set(RXTOOLS_DIR ${CMAKE_CURRENT_LIST_DIR})

file(GLOB RXTOOLS_SOURCES
	${RXTOOLS_DIR}/src/*.hpp
	${RXTOOLS_DIR}/src/*.cpp)

add_library(RxTools STATIC ${RXTOOLS_SOURCES})

find_package(Eigen3 REQUIRED)

set(LAPISGISCMAKE_PATH "${RXTOOLS_DIR}/src/gis/LapisGis.cmake" CACHE PATH "Path to LapisGis.cmake")
include(${LAPISGISCMAKE_PATH})

set(RXTOOLS_EXTERNAL_INCLUDES
	${Eigen_INCLUDE_DIRS}
	${LAPISGIS_INCLUDES}
	)

set(RXTOOLS_EXTERNAL_LINKS
	${Eigen_LIBRARIES}
	${LAPISGIS_LINKS}
	)

target_include_directories(RxTools PRIVATE ${RXTOOLS_EXTERNAL_INCLUDES})
target_precompile_headers(RxTools PRIVATE ${RXTOOLS_DIR}/src/rxtools_pch.hpp)
#Set that up

set(RXTOOLS_INCLUDES
	${RXTOOLS_EXTERNAL_INCLUDES}
	${RXTOOLS_DIR}/src
	)
set(RXTOOLS_LINKS
	${RXTOOLS_EXTERNAL_LINKS}
	RxTools
	)

if (MSVC)
	target_compile_options(RxTools PRIVATE /W3 /WX)
else()
	target_compile_options(RxTools PRIVATE -Wall -Wextra -Werror)
endif()
