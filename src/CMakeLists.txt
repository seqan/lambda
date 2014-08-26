# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/h4nn3s/apps/lambda/CMakeLists.txt
#
# CMakeLists.txt file for lambda.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_h4nn3s_apps_lambda)
message (STATUS "Configuring extras/apps/lambda")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES OpenMP)
find_package(SeqAn REQUIRED)
find_package(CXX11 REQUIRED)
find_package(CXX14)

if (NOT CXX11_FOUND)
    message (STATUS "  C++11 not found: Not building lambda.")
    return ()
endif (NOT CXX11_FOUND)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set (CXX14_CXX_FLAGS "${CXX14_CXX_FLAGS} -ftemplate-depth-1024")
endif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
add_executable (lambda lambda.cpp
                lambda.hpp
                misc.hpp
                match.hpp
                options.hpp
                alph.hpp
                holders.hpp
                trans.hpp)
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
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
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

