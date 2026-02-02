include_guard(GLOBAL)

set(RXTOOLS_DIR ${CMAKE_CURRENT_LIST_DIR})

file(GLOB RXTOOLS_SOURCES
	${RXTOOLS_DIR}/src/*.hpp
	${RXTOOLS_DIR}/src/*.cpp)

if(EXISTS "${RXTOOLS_DIR}/src/LapisGis/LapisGis.cmake")
	include("${RXTOOLS_DIR}/src/LapisGis/LapisGis.cmake")
	include("${RXTOOLS_DIR}/src/lico/lico.cmake")
	include("${RXTOOLS_DIR}/src/ProcessedFolder/ProcessedFolder.cmake")
else()
	include("${RXTOOLS_DIR}/../LapisGis/LapisGis.cmake")
	include("${RXTOOLS_DIR}/../lico/lico.cmake")
	include("${RXTOOLS_DIR}/../ProcessedFolder/ProcessedFolder.cmake")
endif()

find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options)

set(RXTOOLS_EXTERNAL_INCLUDES
	${Eigen_INCLUDE_DIRS}
	${Boost_INCLUDE_DIRS}
	${LAPISGIS_INCLUDES}
	${LICO_INCLUDES}
	${PROCESSEDFOLDER_INCLUDES}
	)

set(RXTOOLS_EXTERNAL_LINKS
	${Eigen_LIBRARIES}
	${Boost_LIBRARIES}
	${LAPISGIS_LINKS}
	${LICO_LINKS}
	${PROCESSEDFOLDER_LINKS}
	)

add_library(RxTools STATIC ${RXTOOLS_SOURCES})
target_include_directories(RxTools PRIVATE ${RXTOOLS_EXTERNAL_INCLUDES})
target_precompile_headers(RxTools PRIVATE ${RXTOOLS_DIR}/src/rxtools_pch.hpp)

set(RXTOOLS_INCLUDES
	${RXTOOLS_EXTERNAL_INCLUDES}
	${RXTOOLS_DIR}/src
	)
set(RXTOOLS_LINKS
	${RXTOOLS_EXTERNAL_LINKS}
	RxTools
	)
