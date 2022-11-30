# ===========================================================================
#                  Lambda
# ===========================================================================

cmake_minimum_required (VERSION 3.4.0)
string(ASCII 27 Esc)
set(ColourBold  "${Esc}[1m")
set(ColourReset "${Esc}[m")
set(ColourRed   "${Esc}[31m")

message ("${ColourBold}Compiler Detection${ColourReset}")

project (lambda3 CXX)

# ----------------------------------------------------------------------------
# Make "Release" the default cmake build type
# ----------------------------------------------------------------------------

if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING
         "Choose the type of build, options are: Debug Release RelWithDebInfo"
         FORCE)
endif ()

# ----------------------------------------------------------------------------
# Begin of dependency detection
# ----------------------------------------------------------------------------

message ("\n${ColourBold}Dependency detection${ColourReset}")

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/submodules/seqan/include/seqan/version.h")
    set (CMAKE_INCLUDE_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/submodules/seqan/include
        ${CMAKE_INCLUDE_PATH})
    set (CMAKE_PREFIX_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/submodules/seqan/util/cmake
        ${CMAKE_PREFIX_PATH})
    set (CMAKE_MODULE_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/submodules/seqan/util/cmake
        ${CMAKE_MODULE_PATH})
    message (STATUS "Found a local SeqAn library provided with the Lambda source code.")
    message (    "   This will be preferred over system global headers.")
endif ()

# CEREAL
if (EXISTS ${CMAKE_SOURCE_DIR}/submodules/cereal)
    set(BUILD_SANDBOX OFF)
    set(BUILD_DOC OFF)
    set(BUILD_TESTS OFF)
    set(SKIP_PERFORMANCE_COMPARISON ON)
    set(CEREAL_INSTALL OFF)
    add_subdirectory(submodules/cereal)
endif ()

# SGS FMindex
add_subdirectory(submodules/fmindex-collection)

# ----------------------------------------------------------------------------
# Add Lambda targets
# ----------------------------------------------------------------------------

add_subdirectory(src)

# ----------------------------------------------------------------------------
# Warn if cmake build type is not "Release"
# ----------------------------------------------------------------------------

if (NOT CMAKE_BUILD_TYPE STREQUAL Release)
    message (STATUS "${ColourRed}CMAKE_BUILD_TYPE is not \"Release\", your binaries will be slow.${ColourReset}")
endif ()

# ----------------------------------------------------------------------------
# Add Tests
# ----------------------------------------------------------------------------

message ("\n${ColourBold}Setting up unit tests${ColourReset}")

enable_testing ()
add_subdirectory (test EXCLUDE_FROM_ALL)

