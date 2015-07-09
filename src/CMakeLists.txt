# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/h4nn3s/apps/lambda/CMakeLists.txt
#
# CMakeLists.txt file for lambda.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (lambda)
message (STATUS "Configuring lambda")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES OpenMP ZLIB BZip2)
find_package(SeqAn REQUIRED)
find_package(CXX11 REQUIRED)
find_package(CXX14)

if (NOT CXX11_FOUND)
    message (FATAL_ERROR "  C++11 not found: Not building lambda.")
    return ()
endif (NOT CXX11_FOUND)

# Warn if OpenMP was not found.
if (NOT SEQAN_HAS_OPENMP)
    message (WARNING "WARNING: OpenMP not found. Lambda will be built without multi-threading! This is probably not what you want!")
endif (NOT SEQAN_HAS_OPENMP)

# Warn if Zlib was not found.
if (NOT SEQAN_HAS_ZLIB)
  message (WARNING "WARNING: Zlib not found. Building lambda without gzipped input and output.")
endif (NOT SEQAN_HAS_ZLIB)

# Warn if BZip2 was not found.
if (NOT SEQAN_HAS_BZIP2)
  message (WARNING "WARNING: BZip2 not found. Building lambda without bzipped input and output.")
endif (NOT SEQAN_HAS_BZIP2)

if (COMPILER_IS_CLANG)
    set (CXX14_CXX_FLAGS "${CXX14_CXX_FLAGS} -ftemplate-depth-1024")
endif (COMPILER_IS_CLANG)

if(CMAKE_COMPILER_IS_GNUCXX)
  if (491 GREATER _GCC_VERSION)
    message (FATAL_ERROR "Your GCC version is too old. Minimum version is GCC-4.9.1!")
    return ()
  endif(491 GREATER _GCC_VERSION)
endif(CMAKE_COMPILER_IS_GNUCXX)

# ----------------------------------------------------------------------------
# App-Level Configuration
# ----------------------------------------------------------------------------

set (SEQAN_APP_VERSION "0.9.0")

option (LAMBDA_FASTBUILD "Set to ON to skip certain rarely used codepaths (speeds up build)." OFF)
if (LAMBDA_FASTBUILD)
    add_definitions (-DFASTBUILD=1)
endif (LAMBDA_FASTBUILD)

option (LAMBDA_NATIVE_BUILD "Set to OFF to disable architecture-specific optimizations, e.g. g++ -march=native." ON)
if (LAMBDA_NATIVE_BUILD)
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -march=native")
endif (LAMBDA_NATIVE_BUILD)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")

# Update the list of file names below if you add source files to your application.
add_executable (lambda lambda.cpp
                lambda.hpp
                misc.hpp
                match.hpp
                options.hpp
                alph.hpp
                holders.hpp)
add_executable (lambda_indexer lambda_indexer.cpp
                lambda_indexer.hpp
                options.hpp
                misc.hpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (lambda ${SEQAN_LIBRARIES})
target_link_libraries (lambda_indexer ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ${CXX11_CXX_FLAGS} ${CXX14_CXX_FLAGS} -Wno-vla")

set(BUILD_SHARED_LIBS OFF)

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (lambda)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install lambda in ${PREFIX}/bin directory
install (TARGETS lambda lambda_indexer
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/lambda for SeqAn release builds.
install (FILES ../LICENSE
               ../README.md
               ../CHANGELOG.md
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
# install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (lambda)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:lambda")
  set (CPACK_PACKAGE_NAME "lambda")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "lambda")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Your Name <your.email@example.net>")
  set (CPACK_PACKAGE_VENDOR "Your Name <your.email@example.net>")

  seqan_configure_cpack_app (lambda "lambda")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:lambda")

