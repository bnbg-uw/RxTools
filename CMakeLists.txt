cmake_minimum_required(VERSION 3.24)
project(RxTools VERSION 1.0)

if(" ${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "
FATAL: In-source builds are not allowed.
       You should create a separate directory for build files.
")
endif()
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR})

include(RxTools.cmake)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/rxtools.cmake
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/.gitignore
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/.gitmodules
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/README.md
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (MSVC)
	target_compile_options(RxTools PRIVATE /W3 /WX)
else()
	target_compile_options(RxTools PRIVATE -Wall -Wextra -Werror)
endif()