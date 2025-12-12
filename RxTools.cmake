include_guard(GLOBAL)

set(RXTOOLS_DIR ${CMAKE_CURRENT_LIST_DIR})

file(GLOB RXTOOLS_SOURCES
	${RXTOOLS_DIR}/src/*.hpp
	${RXTOOLS_DIR}/src/*.cpp)

if(EXISTS "${RXTOOLS_DIR}/src/LapisGis/LapisGis.cmake")
	include("${RXTOOLS_DIR}/src/LapisGis/LapisGis.cmake")
	include("${RXTOOLS_DIR}/src/lico/lico.cmake")
else()
	include("${RXTOOLS_DIR}/../LapisGis/LapisGis.cmake")
	include("${RXTOOLS_DIR}/../lico/lico.cmake")
endif()

find_package(Eigen3 REQUIRED)

set(RXTOOLS_EXTERNAL_INCLUDES
	${Eigen_INCLUDE_DIRS}
	${LAPISGIS_INCLUDES}
	)

set(RXTOOLS_EXTERNAL_LINKS
	${Eigen_LIBRARIES}
	${LAPISGIS_LINKS}
	)

set(RXTOOLS_INCLUDES
	${RXTOOLS_EXTERNAL_INCLUDES}
	${RXTOOLS_DIR}/src
	)
set(RXTOOLS_LINKS
	${RXTOOLS_EXTERNAL_LINKS}
	RxTools
	)
