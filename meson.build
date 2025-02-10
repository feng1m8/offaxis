project(
    'offaxis', 'cpp', 'c',
    meson_version: '>=1.3.0',
    license: 'GPL-3.0-or-later',
    license_files: [
        'LICENSE',
        'src/relxill/LICENSE',
    ],
    default_options: {
        'b_lto': true,
        'b_ndebug': 'if-release',
        'buildtype': 'release',
        'cpp_std': 'c++17',
        'c_std': 'c11',
    },
)

python = import('python').find_installation(pure: false)

headas = run_command(
    python, '-c', 'import os; print(os.environ["HEADAS"])',
    check: false,
).stdout().strip()
headas = include_directories(headas + '/include')

cxx = meson.get_compiler('cpp')
cxx.has_header('xsTypes.h', include_directories: headas, required: true)

ynogk_cxx = subproject('YNOGK-CXX').get_variable('ynogk_cxx')

fftw3 = dependency('fftw3')
cfitsio = dependency('cfitsio')
ccfits = dependency('CCfits')
chealpix = dependency('chealpix')
gsl = dependency('gsl')
openmp = dependency('openmp')

libdl = cxx.find_library('dl')

configure_file(
    input: 'src/relxill/src/config.h.in',
    output: 'config.h',
    copy: true,
)

custom_target(
    output: [
        'xspec_wrapper_lmodels.cpp',
        'xspec_wrapper_lmodels.h',
    ],
    input: [
        'src/relxill/src/create_wrapper_xspec.py',
        'src/relxill/src/modelfiles/lmodel_relxill_public.dat',
    ],
    command: [
        python, '@INPUT@', '@OUTPUT0@',
    ],
    build_by_default: true,
)

relxill = static_library(
    'relxill',
    'src/relxill/src/relutility.c',
    'src/relxill/src/xilltable.c',
    'src/relxill/src/reltable.c',
    'src/relxill/src/writeOutfiles.c',
    'src/relxill/src/donthcomp.c',
    'src/relxill/src/Xillspec.cpp',
    'src/relxill/src/Relphysics.cpp',
    'src/relxill/src/Relcache.cpp',
    'src/relxill/src/Relprofile.cpp',
    'src/relxill/src/Rellp.cpp',
    'src/relxill/src/ModelParams.cpp',
    'src/relxill/src/Relbase.cpp',
    'src/relxill/src/Relreturn_Corona.cpp',
    'src/relxill/src/Relreturn_Datastruct.cpp',
    'src/relxill/src/Relreturn_Table.cpp',
    dependencies: [
        fftw3,
        cfitsio,
    ],
    include_directories: [
        headas,
        'src/relxill/src',
    ],
    build_by_default: false,
)

relxill = declare_dependency(
    objects: [
        relxill.extract_all_objects(recursive: false),
    ],
    include_directories: [
        headas,
        'src/relxill/src',
    ],
    dependencies: [
        cfitsio,
        fftw3,
    ],
)

offaxis = shared_library(
    'offaxis',
    'src/offaxis/envs.cxx',
    'src/offaxis/kyn.cxx',
    'src/offaxis/offaxline.cxx',
    'src/offaxis/raytracing.cxx',
    'src/offaxis/sphere.cxx',
    'src/offaxis/convolve/convolve.cxx',
    'src/offaxis/convolve/offaxconv.cxx',
    'src/offaxis/reflection/spectrum.cxx',
    'src/offaxis/reflection/offaxxill.cxx',
    'src/offaxis/reflection/primary.cxx',
    'src/offaxis/reflection/initialize.cxx',
    'src/offaxis/reflection/emission.cxx',
    'src/offaxis/xspecmodel/xspecmodel.cxx',
    include_directories: [
        'src',
        'include',
    ],
    dependencies: [
        ynogk_cxx,
        relxill,
        ccfits,
        gsl,
        chealpix,
        libdl,
        openmp,
    ],
    install: true,
    install_dir: python.get_install_dir() / 'offaxis',
)

offaxis = declare_dependency(
    link_with: [
        offaxis,
    ],
    include_directories: [
        'include',
    ],
)

add_languages('cython', native: false)

python.extension_module(
    'offaxis',
    'src/offaxis/xspecmodel/offaxis.pyx',
    dependencies: [
        offaxis,
    ],
    cpp_args: [
        '-DCYTHON_LIMITED_API',
    ],
    override_options: [
        'cython_language=cpp',
    ],
    limited_api: '3.7',
    install: true,
    subdir: 'offaxis',
    install_rpath: '$ORIGIN',
)

python.install_sources(
    'src/offaxis/xspecmodel/__init__.py',
    'src/offaxis/xspecmodel/initpackage.py',
    'src/offaxis/xspecmodel/offaxis.dat',
    subdir: 'offaxis',
)
