cmake_minimum_required(VERSION 3.26)

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project("offaxis"
    LANGUAGES CXX
)

file(GENERATE OUTPUT ".gitignore" CONTENT "*")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory("src/relxill" EXCLUDE_FROM_ALL)
add_subdirectory("subprojects/YNOGK-CXX" EXCLUDE_FROM_ALL)

find_library(DL dl REQUIRED)
find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(CCfits REQUIRED IMPORTED_TARGET CCfits)
pkg_check_modules(chealpix REQUIRED IMPORTED_TARGET chealpix)

add_library(
    "offaxis_cxx" SHARED
    "src/offaxis/envs.cxx"
    "src/offaxis/kyn.cxx"
    "src/offaxis/offaxline.cxx"
    "src/offaxis/raytracing.cxx"
    "src/offaxis/sphere.cxx"
    "src/offaxis/convolve/convolve.cxx"
    "src/offaxis/convolve/offaxconv.cxx"
    "src/offaxis/reflection/spectrum.cxx"
    "src/offaxis/reflection/offaxxill.cxx"
    "src/offaxis/reflection/primary.cxx"
    "src/offaxis/reflection/initialize.cxx"
    "src/offaxis/reflection/emission.cxx"
    "src/offaxis/xspecmodel/xspecmodel.cxx"
)

target_include_directories(
    "offaxis_cxx"
    PRIVATE "src"
    PUBLIC "include"
)

target_link_libraries(
    "offaxis_cxx"
    PRIVATE ${DL} Relxill ynogk_cxx GSL::gsl OpenMP::OpenMP_CXX PkgConfig::chealpix PkgConfig::CCfits
    INTERFACE $<TARGET_OBJECTS:ynogk_cxx>
)

find_package(
    Python REQUIRED
    Interpreter Development.SABIModule
)

list(APPEND CMAKE_MODULE_PATH "${Python_SITEARCH}/cython_cmake/cmake")
find_package(Cython MODULE REQUIRED)
include(UseCython)

cython_transpile(
    "src/offaxis/xspecmodel/offaxis.pyx"
    LANGUAGE CXX
    OUTPUT_VARIABLE offaxis_pyx
)

Python_add_library(
    "offaxis" MODULE
    USE_SABI 3.7 WITH_SOABI
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
    FILES "offaxis/__init__.py"
    FILES "offaxis/initpackage.py"
    FILES "offaxis/offaxis.dat"
    DESTINATION "offaxis"
)
