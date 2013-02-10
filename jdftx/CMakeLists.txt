project(JDFTx)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake-Modules/")

#Package configuration:
set(CPACK_PACKAGE_NAME "JDFTx")
set(CPACK_PACKAGE_DESCRIPTION "The playground for density-functional theory")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Electronic and classical density-functional theory software with a modern multi-hardware C++11 code-base")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "preAlpha")
set(VERSION_STRING "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
