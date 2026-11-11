# CMake build for EddyFlow's engine executables
#
#   cmake --preset release
#   cmake --build --preset release
#
# builds eddyflow_rp and eddyflow_fcc into bin/<os>/, where the GUI, the
# release scripts and tests/regression/run.sh look for them. The `checks`
# preset builds the same two with the subscript and runtime checks compiled
# in, into the same place, as `make debug` did: a checks build cannot then be
# mistaken for a release one sitting beside it.
#
# The test fixtures (pwbref, pwbcmp, ccfcheck, covmaxcheck, langfordcheck,
# pwbtablecheck) are not part of the default build. Name them:
#
#   cmake --build --preset release --target pwbref
#
# They land in obj/<os>/, not bin/: bin/ holds what a user runs.
#
# CMake reads every `use` and `include` itself, and Ninja rebuilds whatever a
# changed module or include file reaches, so there is no dependency list to
# keep. The old Makefile's generated one existed because an object with no rule
# was silently linked from its previous compile; that cannot happen here.
#
# Note the executables are NOT byte-reproducible: the linker stamps a
# timestamp into the PE header. Compare object code, or behaviour, when
# checking whether a change is inert.

cmake_minimum_required(VERSION 3.21)

project(EddyFlowEngine LANGUAGES Fortran)

if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR "The engine is built with gfortran; found ${CMAKE_Fortran_COMPILER_ID}.")
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Release" FORCE)
endif()

# Every flag the compiler sees is spelled out below. CMake's per-configuration
# defaults (-O3, -O2 -g, -DNDEBUG ...) are cleared so that nothing is added
# behind them - the fft4 sources in particular are compiled with no
# optimisation at all, and must stay that way.
foreach(config DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
    set(CMAKE_Fortran_FLAGS_${config} "")
endforeach()
# Nor any Windows system import library the engine does not name.
set(CMAKE_Fortran_STANDARD_LIBRARIES "")

# -fbounds-check is NOT in the release flags. It traps a subscript outside its
# array, which is a debugging aid and not something a released binary needs:
# it cannot change any value, only turn a wrong one into a stop. It costs real
# time on the per-sample loops, so it lives in EDDYFLOW_CHECKS instead.
#
# Before it was moved, the whole regression sweep was run with it still on and
# came back 25 passed, 0 failed - so nothing in the shipped paths is relying on
# it to catch anything. Build the `checks` preset and run the sweep again after
# any change that indexes raw data, and keep it that way.
#
# Budget for it: a checks build is about 1.75x slower than release, not the
# 6-7 % that -fbounds-check alone costs. -fcheck=all is the expensive part - it
# adds array temporary, pointer and recursion checks on top of the subscripts.
# Measured on one day of CH-LAE planar fit: release 50 s, -fbounds-check alone
# 54 s, checks 89 s. So the regression sweep, normally around 50 minutes, takes
# about 90 against a checks build.
option(EDDYFLOW_CHECKS "Compile in the subscript and runtime checks" OFF)

# -ffp-contract=off: never fuse a multiply and an add into one FMA. A fused
# product is rounded once instead of twice, so allowing it would move results
# - and the blocked cross-correlation in m_pwb_core is bit-identical to the
# loop it replaced only because each product is rounded exactly as before.
# Today it changes nothing: the default x86-64 target has no FMA instruction.
# It is here so that adding -march one day cannot silently change the numbers.
set(EDDYFLOW_FFLAGS
    -Wno-surprising -Wall -pedantic -fall-intrinsics -Wextra
    -Wcharacter-truncation -Wunderflow -Wno-compare-reals
    -std=f2008 -O3 -ffp-contract=off)
set(EDDYFLOW_CHECK_FLAGS -fbounds-check -fcheck=all -g -fbacktrace)
if(EDDYFLOW_CHECKS)
    list(APPEND EDDYFLOW_FFLAGS ${EDDYFLOW_CHECK_FLAGS})
endif()

# libgfortran and libquadmath linked in. On Windows that is what lets
# src/src_os/win/setlocale_override.f90 replace the setlocale they call around
# every formatted READ and WRITE, whose saved locale string the UCRT has freed
# by the time it is restored (see that file); everywhere it means the engine
# needs no Fortran runtime installed next to it.
#
# A cache string, so a platform that needs otherwise can say so: the macOS
# release script passes its own.
set(EDDYFLOW_LDFLAGS "-shared-libgcc -static-libgfortran -static-libquadmath" CACHE STRING
    "Link flags for eddyflow_rp and eddyflow_fcc")
separate_arguments(EDDYFLOW_LDFLAGS_LIST UNIX_COMMAND "${EDDYFLOW_LDFLAGS}")

# The one module that needs a different source per operating system, and the
# folders the binaries go to.
if(WIN32)
    set(OS_FAMILY win)
    set(OS_DIR win)
elseif(APPLE)
    set(OS_FAMILY posix)
    set(OS_DIR mac)
else()
    set(OS_FAMILY posix)
    set(OS_DIR linux)
endif()
message(STATUS "Compiling for ${OS_DIR} (src/src_os/${OS_FAMILY})")

# Overridable, for eddyflow-build-script, which builds into its own folders.
set(EDDYFLOW_BIN_DIR "${PROJECT_SOURCE_DIR}/bin/${OS_DIR}" CACHE PATH
    "Where eddyflow_rp and eddyflow_fcc go")
set(EDDYFLOW_FIXTURE_DIR "${PROJECT_SOURCE_DIR}/obj/${OS_DIR}" CACHE PATH
    "Where the test drivers go; the static checks look in obj/<os>")

# All modules in one folder, as gfortran's -J put them before; every target
# reads them from there.
set(EDDYFLOW_MODULE_DIR "${CMAKE_BINARY_DIR}/mod")

# Sources: every file in each folder, as the Makefile's $(wildcard) took them.
# CONFIGURE_DEPENDS re-globs at build time, so a new file needs no re-run of
# cmake by hand.
set(SRC "${PROJECT_SOURCE_DIR}/src")
file(GLOB COMMON_SOURCES CONFIGURE_DEPENDS "${SRC}/src_common/*.f90")
file(GLOB OS_SOURCES     CONFIGURE_DEPENDS "${SRC}/src_os/${OS_FAMILY}/*.f90")
file(GLOB RP_SOURCES     CONFIGURE_DEPENDS "${SRC}/src_rp/*.f90")
file(GLOB FFT4_SOURCES   CONFIGURE_DEPENDS "${SRC}/src_rp/fft4/*.F")
file(GLOB FCC_SOURCES    CONFIGURE_DEPENDS "${SRC}/src_fcc/*.f90")
set(RP_MAIN  "${SRC}/src_rp/eddyflow-rp_main.f90")
set(FCC_MAIN "${SRC}/src_fcc/eddyflow-fcc_main.f90")
list(REMOVE_ITEM RP_SOURCES "${RP_MAIN}")
list(REMOVE_ITEM FCC_SOURCES "${FCC_MAIN}")

# The free-form sources are not preprocessed, and never were. Left to itself,
# CMake's Ninja generator runs cpp over every source before compiling it, and
# cpp predefines names such as `linux` and `unix` on those systems - so a
# variable of that name would quietly become a 1.
function(eddyflow_fortran_target target)
    set_target_properties(${target} PROPERTIES
        Fortran_MODULE_DIRECTORY "${EDDYFLOW_MODULE_DIR}"
        Fortran_PREPROCESS OFF)
    target_compile_options(${target} PRIVATE ${EDDYFLOW_FFLAGS})
endfunction()

# Shared by both executables: compiled once, linked twice.
add_library(engine_common OBJECT ${COMMON_SOURCES} ${OS_SOURCES})
eddyflow_fortran_target(engine_common)

add_library(engine_rp_objs OBJECT ${RP_SOURCES})
eddyflow_fortran_target(engine_rp_objs)
target_link_libraries(engine_rp_objs PUBLIC engine_common)

# FFTPACK, fixed form, compiled as it always was: warnings off and no
# optimisation flags at all.
add_library(engine_fft4 OBJECT ${FFT4_SOURCES})
target_compile_options(engine_fft4 PRIVATE -w)

add_library(engine_fcc_objs OBJECT ${FCC_SOURCES})
eddyflow_fortran_target(engine_fcc_objs)
target_link_libraries(engine_fcc_objs PUBLIC engine_common)
# One writer packs fixed-width FLUXNET fields and truncates by design.
set_source_files_properties("${SRC}/src_fcc/write_out_fluxnet_fcc.f90"
    PROPERTIES COMPILE_OPTIONS "-Wno-character-truncation")

function(eddyflow_executable target output_name)
    set_target_properties(${target} PROPERTIES
        OUTPUT_NAME ${output_name}
        RUNTIME_OUTPUT_DIRECTORY "${EDDYFLOW_BIN_DIR}")
    target_link_options(${target} PRIVATE ${EDDYFLOW_LDFLAGS_LIST})
endfunction()

add_executable(eddyflow_rp "${RP_MAIN}"
    $<TARGET_OBJECTS:engine_rp_objs>
    $<TARGET_OBJECTS:engine_fft4>
    $<TARGET_OBJECTS:engine_common>)
eddyflow_fortran_target(eddyflow_rp)
add_dependencies(eddyflow_rp engine_rp_objs)
eddyflow_executable(eddyflow_rp eddyflow_rp)

add_executable(eddyflow_fcc "${FCC_MAIN}"
    $<TARGET_OBJECTS:engine_fcc_objs>
    $<TARGET_OBJECTS:engine_common>)
eddyflow_fortran_target(eddyflow_fcc)
add_dependencies(eddyflow_fcc engine_fcc_objs)
eddyflow_executable(eddyflow_fcc eddyflow_fcc)

# Test fixtures. None is built by default, each lands in obj/<os>/, and each is
# linked -static: it is run by a test harness that has no reason to know where
# libgfortran lives, and a driver that reports the wrong numbers is a worse
# failure than one that is a megabyte larger.
#
# Most link against a named handful of the engine's objects and nothing else,
# picked out of engine_common by name. That is the point of them: pwbref, for
# one, proves the pre-whitening chain does not depend on engine state, and it
# stops linking the moment that stops being true.
function(eddyflow_fixture target output_name main)
    add_executable(${target} EXCLUDE_FROM_ALL "${SRC}/src_tools/${main}" ${ARGN})
    eddyflow_fortran_target(${target})
    add_dependencies(${target} engine_common)
    set_target_properties(${target} PROPERTIES
        OUTPUT_NAME ${output_name}
        RUNTIME_OUTPUT_DIRECTORY "${EDDYFLOW_FIXTURE_DIR}")
    target_link_options(${target} PRIVATE -static)
endfunction()

# The objects of the named engine_common sources, and only those.
function(eddyflow_common_objects out)
    list(JOIN ARGN "|" names)
    set(${out} "$<FILTER:$<TARGET_OBJECTS:engine_common>,INCLUDE,[/\\\\](${names})\\.f90\\.o(bj)?$>" PARENT_SCOPE)
endfunction()

# The PWB reference driver: the pre-whitening chain alone.
# static_checks/test_pwb_reference_static.py builds and runs it.
eddyflow_common_objects(pwb_core_objects m_numeric_kinds m_pwb_core)
eddyflow_fixture(pwbref pwb_reference pwb_reference_main.f90 ${pwb_core_objects})

# The PWB comparison driver: the engine's whole full-rate chain for one
# period - pre-whitening, the four block bootstraps, mode, HDI, choice of
# combination - on block starts it is handed, so tests/pwb_dyco/compare.py can
# give dyco the same resamples and compare every number.
eddyflow_fixture(pwbcmp pwb_compare pwb_compare_main.f90 ${pwb_core_objects})

# The cross-correlation bit-exactness check: the pre-rewrite ComputeCcfWindow
# frozen, compared on raw bits with the shipped m_pwb_core.
# static_checks/test_pwb_ccf_blocked_static.py runs it.
eddyflow_fixture(ccfcheck ccf_bitexact ccf_bitexact_main.f90 ${pwb_core_objects})

# The covariance-maximisation bit-exactness check, built like ccfcheck: the
# old per-lag path frozen, compared on raw bits with the shipped module.
# static_checks/test_covmax_series_static.py runs it.
eddyflow_common_objects(covmax_core_objects m_numeric_kinds m_covmax_core)
eddyflow_fixture(covmaxcheck covmax_bitexact covmax_bitexact_main.f90 ${covmax_core_objects})

# The Langford et al. (2015) far-lag noise check: FarLagCovarianceNoise on
# synthetic series with known answers. static_checks/test_langford_lod_static.py
# runs it.
eddyflow_fixture(langfordcheck langford_check langford_check_main.f90 ${covmax_core_objects})

# The PWB table check: store a year of rows through the shipped routines and
# confirm every lookup returns the row a scan from row 1 would. Unlike the
# others it needs the engine's own state, so it links every RP object but the
# main program. static_checks/test_pwb_table_static.py runs it.
eddyflow_fixture(pwbtablecheck pwb_table_check pwb_table_check_main.f90
    $<TARGET_OBJECTS:engine_rp_objs>
    $<TARGET_OBJECTS:engine_fft4>
    $<TARGET_OBJECTS:engine_common>)
add_dependencies(pwbtablecheck engine_rp_objs)
