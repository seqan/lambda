# ===========================================================================
#                  Lambda
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
string(ASCII 27 Esc)
set(ColourBold  "${Esc}[1m")
set(ColourReset "${Esc}[m")
message (STATUS "${ColourBold}Compiler Detection...${ColourReset}")

project (lambda CXX)

# ----------------------------------------------------------------------------
# Load SeqAn build system
# ----------------------------------------------------------------------------

message (STATUS "${ColourBold}Configuring SeqAn...${ColourReset}")
set (CMAKE_INCLUDE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/include
    ${CMAKE_INCLUDE_PATH})
set (CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/util/cmake/FindTBB
    ${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/util/cmake
    ${CMAKE_MODULE_PATH})
include (SeqAnContribs)

# Build in release mode by default
if (CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release")
endif()

set (SEQAN_USE_SEQAN_BUILD_SYSTEM TRUE CACHE INTERNAL "Use SeqAn build system." FORCE)
include (SeqAnBuildSystem)
set (SEQAN_ROOT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/" CACHE INTERNAL "Root source directory." FORCE)
seqan_build_system_init ()
seqan_get_repository_info ()
seqan_setup_library ()
include (package)
include (SeqAnCtdSetup)

# ----------------------------------------------------------------------------
# Add Lambda targets
# ----------------------------------------------------------------------------

message (STATUS "${ColourBold}Configuring LAMBDA...${ColourReset}")
add_subdirectory(src)

# ----------------------------------------------------------------------------
# Add Tests
# ----------------------------------------------------------------------------

message (STATUS "${ColourBold}Configuring LAMBDA Tests...${ColourReset}")
add_subdirectory(tests)
