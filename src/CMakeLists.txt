# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/h4nn3s/apps/lambda/CMakeLists.txt
#
# CMakeLists.txt file for lambda.
# ===========================================================================

# ----------------------------------------------------------------------------
# App version
# ----------------------------------------------------------------------------

# change this after every release
set (SEQAN_APP_VERSION_MAJOR "2")
set (SEQAN_APP_VERSION_MINOR "0")
set (SEQAN_APP_VERSION_PATCH "0")

# don't change the following
set (SEQAN_APP_VERSION "${SEQAN_APP_VERSION_MAJOR}.${SEQAN_APP_VERSION_MINOR}.${SEQAN_APP_VERSION_PATCH}")

# adapt when necessary
set (MINIMUM_SEQAN_VERSION "2.3.1")

# ----------------------------------------------------------------------------
# Dependencies (continued)
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
find_package(OpenMP QUIET)
find_package(ZLIB   QUIET)
find_package(BZip2  QUIET)
find_package(SeqAn  QUIET REQUIRED CONFIG)

message(STATUS "These dependencies were found:")
message(   "     OPENMP     ${OPENMP_FOUND}      ${OpenMP_CXX_FLAGS}")
message(   "     ZLIB       ${ZLIB_FOUND}      ${ZLIB_VERSION_STRING}")
message(   "     BZIP2      ${BZIP2_FOUND}      ${BZIP2_VERSION_STRING}")
message(   "     SEQAN      ${SEQAN_FOUND}      ${SEQAN_VERSION_STRING}")

# Warn if OpenMP was not found.
if (NOT OPENMP_FOUND)
    message (WARNING "WARNING WARNING WARNING\nWARNING: OpenMP not found. Lambda will be built without multi-threading! "
    "This is probably not what you want! Use GCC >= 4.9.1, Clang >= 3.8.0 or ICC >= 16.0.2\nWARNING WARNING WARNING")
endif (NOT OPENMP_FOUND)

# Warn if Zlib was not found.
if (NOT ZLIB_FOUND)
  message (WARNING "WARNING: Zlib not found. Building lambda without support for gzipped input and output (this includes support for .bam).")
endif (NOT ZLIB_FOUND)

# Warn if BZip2 was not found.
if (NOT BZIP2_FOUND)
  message (WARNING "WARNING: BZip2 not found. Building lambda without support for bzipped input and output.")
endif (NOT BZIP2_FOUND)

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9.1)
        message (FATAL_ERROR "Your GCC version is too old. Minimum version is GCC-4.9.1!")
        return ()
    endif ()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.6)
        message (FATAL_ERROR "Your Clang version is too old. Please upgrade to 3.8.0 or use GCC.")
        return()
    elseif (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.8)
        message (WARNING "Your Clang version is too old, you will not have parallelism! Upgrade to 3.8.0 or newer.")
    endif ()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 17)
        message(WARNING "Your Intel Compiler version is too old. Please upgrade to 17.0.0 or newer!")
    endif ()
else ()
    message(WARNING "Unknown compiler, you are own your own!")
endif ()

if (SEQAN_VERSION_STRING VERSION_LESS "${MINIMUM_SEQAN_VERSION}")
    message (FATAL_ERROR "The minimum SeqAn verison required is ${MINIMUM_SEQAN_VERSION}!")
    return ()
endif ()

message(STATUS "The requirements were met.")

# ----------------------------------------------------------------------------
# App-Level Configuration
# ----------------------------------------------------------------------------

message ("\n${ColourBold}Build configuration${ColourReset}")

message (STATUS "LAMBDA version is: ${SEQAN_APP_VERSION}")

option (LAMBDA_FASTBUILD                "Build only blastp and blastx modes (speeds up build)."                                     OFF)
option (LAMBDA_LINGAPS_OPT              "Add optimized codepaths for linear gap costs (increases bin size and compile time)."       OFF)
option (LAMBDA_LONG_PROTEIN_SUBJ_SEQS   "Make max protein sequence length == 4.3billion instead of 65,535. INVALIDATES INDEXS!"     OFF)
option (LAMBDA_MMAPPED_DB               "Use mmapped access to the database."                                                       OFF)
option (LAMBDA_NATIVE_BUILD             "Architecture-specific optimizations, i.e. g++ -march=native."                              ON)
option (LAMBDA_MULTIOPT_BUILD           "Build and install multiple binaries with different optimisation levels + a dispatcher."    OFF)
option (LAMBDA_STATIC_BUILD             "Include all libraries in the binaries."                                                    OFF)

if (LAMBDA_FASTBUILD)
    add_definitions (-DFASTBUILD=1)
endif (LAMBDA_FASTBUILD)

if (LAMBDA_NATIVE_BUILD)
    add_definitions (-DLAMBDA_NATIVE_BUILD=1)
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -march=native")
    if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -xHOST -ipo -no-prec-div -fp-model fast=2")
    endif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
endif (LAMBDA_NATIVE_BUILD)

if (LAMBDA_STATIC_BUILD)
    add_definitions (-DLAMBDA_STATIC_BUILD=1)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    # apple does not support fully static builds, but at least libgcc and libstdc++
    if (APPLE)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
        message (WARNING "WARNING: Builds on Mac are never fully static.")
    else (APPLE)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
    endif (APPLE)
    # on linux cmake adds -rdynamic automatically which clang can't handle in static builds
    if (CMAKE_SYSTEM_NAME MATCHES "Linux")
        SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")
    endif (CMAKE_SYSTEM_NAME MATCHES "Linux")
endif (LAMBDA_STATIC_BUILD)

if (LAMBDA_MMAPPED_DB)
    add_definitions (-DLAMBDA_MMAPPED_DB=1)
endif (LAMBDA_MMAPPED_DB)

if (LAMBDA_LINGAPS_OPT)
    add_definitions (-DLAMBDA_LINGAPS_OPT=1)
endif ()

if (LAMBDA_LONG_PROTEIN_SUBJ_SEQS)
    add_definitions (-DLAMBDA_LONG_PROTEIN_SUBJ_SEQS=1)
endif ()

message(STATUS "The following options are selected for the build:")
message(   "     LAMBDA_FASTBUILD                ${LAMBDA_FASTBUILD}")
message(   "     LAMBDA_LINGAPS_OPT              ${LAMBDA_LINGAPS_OPT}")
message(   "     LAMBDA_LONG_PROTEIN_SUBJ_SEQS   ${LAMBDA_LONG_PROTEIN_SUBJ_SEQS}")
message(   "     LAMBDA_MMAPPED_DB               ${LAMBDA_MMAPPED_DB}")
message(   "     LAMBDA_NATIVE_BUILD             ${LAMBDA_NATIVE_BUILD}")
message(   "     LAMBDA_MULTIOPT_BUILD           ${LAMBDA_MULTIOPT_BUILD}")
message(   "     LAMBDA_STATIC_BUILD             ${LAMBDA_STATIC_BUILD}")
message(STATUS "Run 'cmake -LH' to get a comment on each option.")
message(STATUS "Remove CMakeCache.txt and re-run cmake with -DOPTIONNAME=ON|OFF to change an option.")

# deactivate the version check on broken seqan releases
if (SEQAN_VERSION_STRING VERSION_LESS "2.3.2")
    add_definitions (-DSEQAN_DISABLE_VERSION_CHECK="YES")
endif ()

# ----------------------------------------------------------------------------
# Compiler specifics
# ----------------------------------------------------------------------------

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -ftemplate-depth-1024")

    # do not warn for variable length arrays
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -Wno-vla-extension")
endif ()

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    # do not warn for variable length arrays
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -Wno-vla")

    # parallelize parts of build even for one translation unit
    if (NOT DEFINED LAMBDA_COMPILE_THREADS)
        include(ProcessorCount)
        ProcessorCount(LAMBDA_COMPILE_THREADS)
    endif ()
    if (LAMBDA_COMPILE_THREADS GREATER 1)
        set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto=${LAMBDA_COMPILE_THREADS}")
    endif()

    # strip binaries to make them smaller
    set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -s")
endif ()

if (NOT CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    # -Wextra not met by current SeqAn on intel
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -Wextra")
endif ()

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Enable global exception handler for all seqan apps.
set (SEQAN_DEFINITIONS ${SEQAN_DEFINITIONS} -DSEQAN_GLOBAL_EXCEPTION_HANDLER=1)

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")
add_definitions (-DCMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE}")

# Set the right output directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} -Wall -pedantic")

# Update the list of file names below if you add source files to your application.
set (LAMBDA_SOURCE_FILES
                lambda.cpp
                    shared_definitions.hpp
                    shared_misc.hpp
                    shared_options.hpp
                search.hpp
                    search_algo.hpp
                    search_datastructures.hpp
                    search_misc.hpp
                    search_output.hpp
                    search_options.hpp
                mkindex.hpp
                    mkindex_algo.hpp
                    mkindex_misc.hpp
                    mkindex_options.hpp
                    mkindex_saca.hpp)

add_executable (lambda2 ${LAMBDA_SOURCE_FILES})
target_link_libraries (lambda2 ${SEQAN_LIBRARIES})

if (LAMBDA_MULTIOPT_BUILD)
    add_executable (lambda2-sse4 ${LAMBDA_SOURCE_FILES})
    target_link_libraries (lambda2-sse4 ${SEQAN_LIBRARIES})
    set_target_properties (lambda2-sse4 PROPERTIES COMPILE_FLAGS "-mmmx -msse -msse2 -msse3 -mssse3 -msse4 -mpopcnt")
endif ()

# ----------------------------------------------------------------------------
# Man-pages
# ----------------------------------------------------------------------------

# Umbrella man-page
add_custom_command (OUTPUT lambda2.1
                    COMMAND lambda2 --export-help=man > lambda2.1
                    DEPENDS lambda2)
# searchn subcommand
add_custom_command (OUTPUT lambda2-searchn.1
                    COMMAND lambda2 searchn --export-help=man > lambda2-searchn.1
                    DEPENDS lambda2)
# searchp subcommand
add_custom_command (OUTPUT lambda2-searchp.1
                    COMMAND lambda2 searchp --export-help=man > lambda2-searchp.1
                    DEPENDS lambda2)
# mkindexn subcommand
add_custom_command (OUTPUT lambda2-mkindexn.1
                    COMMAND lambda2 mkindexn --export-help=man > lambda2-mkindexn.1
                    DEPENDS lambda2)
# mkindexp subcommand
add_custom_command (OUTPUT lambda2-mkindexp.1
                    COMMAND lambda2 mkindexp --export-help=man > lambda2-mkindexp.1
                    DEPENDS lambda2)

add_custom_target (manual
                   ALL
                   DEPENDS lambda2.1 lambda2-searchn.1 lambda2-searchp.1 lambda2-mkindexn.1 lambda2-mkindexp.1)

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Adapt to system paths
include (GNUInstallDirs)

if (LAMBDA_MULTIOPT_BUILD)
    # Install lambda binaries into LIBEXECDIR
    install (TARGETS lambda2 lambda2-sse4
             DESTINATION ${CMAKE_INSTALL_LIBEXECDIR})

    # Install dispatcher script into BINDIR
    configure_file(../bin/lambda2.in
                   ${CMAKE_CURRENT_BINARY_DIR}/lambda2-dispatch
                   @ONLY)
    install (FILES ${CMAKE_CURRENT_BINARY_DIR}/lambda2-dispatch
             DESTINATION ${CMAKE_INSTALL_BINDIR}
             RENAME lambda2
             PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
elseif ()
    # Install lambda into BINDIR, usually ${PREFIX}/bin
    install (TARGETS lambda2
             DESTINATION ${CMAKE_INSTALL_BINDIR})
endif ()

# Install non-binary files for the package to DOCDIR, usually ${PREFIX}/share/doc/lambda2
install (FILES ../LICENSE.rst
               ../LICENSE-BSD.rst
               ../LICENSE-AGPL3.rst
               ../README.rst
         DESTINATION ${CMAKE_INSTALL_DOCDIR})

# Man pages into MANDIR, usually ${PREFIX}/share/man/man1 (or without share)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/lambda2.1
               ${CMAKE_CURRENT_BINARY_DIR}/lambda2-searchn.1
               ${CMAKE_CURRENT_BINARY_DIR}/lambda2-searchp.1
               ${CMAKE_CURRENT_BINARY_DIR}/lambda2-mkindexn.1
               ${CMAKE_CURRENT_BINARY_DIR}/lambda2-mkindexp.1
         DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

# Information
set (CPACK_PACKAGE_NAME "lambda2")
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "lambda -- the local aligner for massive bioligical data")
set (CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/../README.rst")
set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE.rst")
set (CPACK_PACKAGE_VENDOR "Hannes Hauswedell <hannes.hauswedell@fu-berlin.de>")
set (CPACK_PACKAGE_CONTACT "${CPACK_PACKAGE_VENDOR}")
set (CPACK_PACKAGE_VERSION_MAJOR "${SEQAN_APP_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${SEQAN_APP_VERSION_MINOR}")
set (CPACK_PACKAGE_VERSION_PATCH "${SEQAN_APP_VERSION_PATCH}")
set (CPACK_PACKAGE_VERSION "${SEQAN_APP_VERSION}")
set (CPACK_PACKAGE_INSTALL_DIRECTORY "${CPACK_PACKAGE_NAME} ${CPACK_PACKAGE_VERSION}")

# Package format(s)
if (CMAKE_SYSTEM_NAME MATCHES "Windows")
    set(CPACK_GENERATOR "ZIP;NSIS")
elseif (CMAKE_SYSTEM_NAME MATCHES "Darwin")
    set(CPACK_GENERATOR "ZIP;DragNDrop")
elseif (CMAKE_VERSION VERSION_LESS "3.1") # TXZ support since 3.1
    set(CPACK_GENERATOR "TBZ2")
else()
    set(CPACK_GENERATOR "TXZ")
endif ()

if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    set(CPACK_GENERATOR "${CPACK_GENERATOR};DEB;RPM")
endif ()

# Package architecture
if (CMAKE_SYSTEM_PROCESSOR STREQUAL "amd64")
    set(CMAKE_SYSTEM_PROCESSOR "x86_64")
    set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
endif ()
if (CMAKE_CXX_FLAGS MATCHES "avx2")
    set (CMAKE_SYSTEM_PROCESSOR "${CMAKE_SYSTEM_PROCESSOR}_avx2")
elseif (CMAKE_CXX_FLAGS MATCHES "sse4")
    set (CMAKE_SYSTEM_PROCESSOR "${CMAKE_SYSTEM_PROCESSOR}_sse4")
endif()

# Include architecture in package name
if (NOT DEFINED CPACK_SYSTEM_NAME)
  set(CPACK_SYSTEM_NAME "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
endif (NOT DEFINED CPACK_SYSTEM_NAME)

include (CPack)
