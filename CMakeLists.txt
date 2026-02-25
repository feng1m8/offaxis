cmake_minimum_required(VERSION 3.26)

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project("offaxis"
    LANGUAGES CXX
    VERSION "2.6"
)

file(GENERATE OUTPUT ".gitignore" CONTENT "*")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory("src/relxill" EXCLUDE_FROM_ALL)
add_subdirectory("src/YNOGK-CXX" EXCLUDE_FROM_ALL)

find_library(DL dl REQUIRED)
find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(CCfits REQUIRED IMPORTED_TARGET CCfits)
pkg_check_modules(chealpix REQUIRED IMPORTED_TARGET chealpix)

add_library(
    "offaxis_cxx" SHARED
    "src/offaxis/offaxline/envs.cxx"
    "src/offaxis/offaxline/kbhtables.cxx"
    "src/offaxis/offaxline/offaxline.cxx"
    "src/offaxis/offaxline/raytracing.cxx"
    "src/offaxis/offaxline/sphere.cxx"
    "src/offaxis/offaxconv/convolve.cxx"
    "src/offaxis/offaxconv/offaxconv.cxx"
    "src/offaxis/offaxxill/spectrum.cxx"
    "src/offaxis/offaxxill/offaxxill.cxx"
    "src/offaxis/offaxxill/doppler.cxx"
    "src/offaxis/offaxxill/xilltable.cxx"
    "src/offaxis/offaxxill/emission.cxx"
    "src/offaxis/xspecmodel.cxx"
)

target_include_directories(
    "offaxis_cxx"
    PUBLIC "src"
    PRIVATE "src/offaxis"
    PUBLIC "include"
)

target_link_libraries(
    "offaxis_cxx"
    PRIVATE ${DL}
    PRIVATE Relxill
    PRIVATE ynogk_cxx
    PRIVATE GSL::gsl
    PRIVATE OpenMP::OpenMP_CXX
    PRIVATE PkgConfig::chealpix
    PRIVATE PkgConfig::CCfits
)

find_package(
    Python REQUIRED
    Interpreter Development.SABIModule
)

list(APPEND CMAKE_MODULE_PATH "${Python_SITEARCH}/cython_cmake/cmake")
find_package(Cython MODULE REQUIRED)
include(UseCython)

cython_transpile(
    "python/offaxis.pyx"
    LANGUAGE CXX
    OUTPUT_VARIABLE offaxis_pyx
)

Python_add_library(
    "offaxis" MODULE
    USE_SABI 3.8 WITH_SOABI
    ${offaxis_pyx}
)

set_target_properties("offaxis" PROPERTIES INSTALL_RPATH "$ORIGIN")
target_link_libraries("offaxis" PRIVATE "offaxis_cxx")

install(
    TARGETS "offaxis_cxx"
    RUNTIME DESTINATION "offaxis"
    LIBRARY DESTINATION "offaxis"
)

install(
    TARGETS "offaxis"
    DESTINATION "offaxis"
)

install(
    FILES "python/__init__.py"
    FILES "python/initpackage.py"
    FILES "include/offaxis/offaxis.dat"
    DESTINATION "offaxis"
)
