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
set (SEQAN_APP_VERSION_MAJOR "3")
set (SEQAN_APP_VERSION_MINOR "0")
set (SEQAN_APP_VERSION_PATCH "0")

# don't change the following
set (SEQAN_APP_VERSION "${SEQAN_APP_VERSION_MAJOR}.${SEQAN_APP_VERSION_MINOR}.${SEQAN_APP_VERSION_PATCH}")

# adapt when necessary
set (MINIMUM_SEQAN_VERSION "2.3.1")
set (MINIMUM_SEQAN3_VERSION "3.0.0")

# ----------------------------------------------------------------------------
# App-Level Configuration
# ----------------------------------------------------------------------------

message ("\n${ColourBold}Build configuration${ColourReset}")

message (STATUS "LAMBDA version is: ${SEQAN_APP_VERSION}")

option (LAMBDA_STATIC_BUILD "Include all libraries in the binaries." OFF)

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


message(STATUS "The following options are selected for the build:")
message(   "     LAMBDA_STATIC_BUILD             ${LAMBDA_STATIC_BUILD}")
message(STATUS "Run 'cmake -LH' to get a comment on each option.")
message(STATUS "Remove CMakeCache.txt and re-run cmake with -DOPTIONNAME=ON|OFF to change an option.")

# deactivate the version check on broken seqan releases
if (SEQAN_VERSION_STRING VERSION_LESS "2.3.2")
    add_definitions (-DSEQAN_DISABLE_VERSION_CHECK="YES")
endif ()

# ----------------------------------------------------------------------------
# Dependencies (continued)
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
find_package(OpenMP QUIET)
find_package(ZLIB   QUIET)
find_package(BZip2  QUIET)
find_package(SeqAn  QUIET REQUIRED CONFIG)
find_package(SeqAn3 QUIET REQUIRED CONFIG)

# SHARG
find_package (sharg QUIET REQUIRED HINTS "${CMAKE_SOURCE_DIR}/submodules/sharg-parser/build_system")

message(STATUS "These dependencies were found:")
message(   "     OPENMP     ${OPENMP_FOUND}      ${OpenMP_CXX_FLAGS}")
message(   "     ZLIB       ${ZLIB_FOUND}      ${ZLIB_VERSION_STRING}")
message(   "     BZIP2      ${BZIP2_FOUND}      ${BZIP2_VERSION_STRING}")
message(   "     SEQAN      ${SEQAN_FOUND}      ${SEQAN_VERSION_STRING}")
message(   "     SEQAN3     ${SEQAN3_FOUND}      ${SEQAN3_VERSION}")
message(   "     SHARG      ${SHARG_FOUND}      ${SHARG_VERSION}")

# Warn if Zlib was not found.
if (NOT ZLIB_FOUND)
  message (WARNING "WARNING: Zlib not found. Building lambda without support for gzipped input and output (this includes support for .bam).")
endif (NOT ZLIB_FOUND)

# Warn if BZip2 was not found.
if (NOT BZIP2_FOUND)
  message (WARNING "WARNING: BZip2 not found. Building lambda without support for bzipped input and output.")
endif (NOT BZIP2_FOUND)

if (SEQAN_VERSION_STRING VERSION_LESS "${MINIMUM_SEQAN_VERSION}")
    message (FATAL_ERROR "The minimum SeqAn version required is ${MINIMUM_SEQAN_VERSION}!")
    return ()
endif ()

if (SEQAN3_VERSION VERSION_LESS "${MINIMUM_SEQAN3_VERSION}")
    message (FATAL_ERROR "The minimum SeqAn3 version required is ${MINIMUM_SEQAN3_VERSION}!")
    return ()
endif ()

message(STATUS "The requirements were met.")

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
set (LAMBDA_SOURCE_FILES lambda.cpp  search.cpp mkindex.cpp)

add_executable (lambda3 ${LAMBDA_SOURCE_FILES})
target_link_libraries (lambda3 ${SEQAN_LIBRARIES} seqan3::seqan3 fmindex-collection)
set_target_properties (lambda3 PROPERTIES COMPILE_FLAGS "-mmmx -msse -msse2 -msse3 -mssse3 -msse4 -mpopcnt")

# ----------------------------------------------------------------------------
# Man-pages
# ----------------------------------------------------------------------------

# Umbrella man-page
add_custom_command (OUTPUT lambda3.1
                    COMMAND lambda3 --export-help man > lambda3.1
                    DEPENDS lambda3)
# searchn subcommand
add_custom_command (OUTPUT lambda3-searchn.1
                    COMMAND lambda3 searchn --export-help man > lambda3-searchn.1
                    DEPENDS lambda3)
# searchp subcommand
add_custom_command (OUTPUT lambda3-searchp.1
                    COMMAND lambda3 searchp --export-help man > lambda3-searchp.1
                    DEPENDS lambda3)
# mkindexn subcommand
add_custom_command (OUTPUT lambda3-mkindexn.1
                    COMMAND lambda3 mkindexn --export-help man > lambda3-mkindexn.1
                    DEPENDS lambda3)
# mkindexp subcommand
add_custom_command (OUTPUT lambda3-mkindexp.1
                    COMMAND lambda3 mkindexp --export-help man > lambda3-mkindexp.1
                    DEPENDS lambda3)

add_custom_target (manual
                   ALL
                   DEPENDS lambda3.1 lambda3-searchn.1 lambda3-searchp.1 lambda3-mkindexn.1 lambda3-mkindexp.1)

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Adapt to system paths
include (GNUInstallDirs)

# Install lambda binaries into LIBEXECDIR
install (TARGETS lambda3
         DESTINATION ${CMAKE_INSTALL_LIBEXECDIR})

# Install dispatcher script into BINDIR
configure_file(../bin/lambda3.in
               ${CMAKE_CURRENT_BINARY_DIR}/lambda3-dispatch
               @ONLY)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/lambda3-dispatch
         DESTINATION ${CMAKE_INSTALL_BINDIR}
         RENAME lambda3
         PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# Install non-binary files for the package to DOCDIR, usually ${PREFIX}/share/doc/lambda3
install (FILES ../LICENSE.rst
               ../LICENSE-BSD.rst
               ../LICENSE-AGPL3.rst
               ../README.rst
         DESTINATION ${CMAKE_INSTALL_DOCDIR})

# Man pages into MANDIR, usually ${PREFIX}/share/man/man1 (or without share)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/lambda3.1
               ${CMAKE_CURRENT_BINARY_DIR}/lambda3-searchn.1
               ${CMAKE_CURRENT_BINARY_DIR}/lambda3-searchp.1
               ${CMAKE_CURRENT_BINARY_DIR}/lambda3-mkindexn.1
               ${CMAKE_CURRENT_BINARY_DIR}/lambda3-mkindexp.1
         DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

# Information
set (CPACK_PACKAGE_NAME "lambda3")
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "lambda -- the local aligner for massive bioligical data")
set (CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/../README.rst")
set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE.rst")
set (CPACK_PACKAGE_VENDOR "Hannes Hauswedell <hannes.hauswedell@decode.is>")
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

# Include architecture in package name
if (NOT DEFINED CPACK_SYSTEM_NAME)
  set(CPACK_SYSTEM_NAME "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
endif (NOT DEFINED CPACK_SYSTEM_NAME)

include (CPack)
