# ===========================================================================
#                  Lambda
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
string(ASCII 27 Esc)
set(ColourBold  "${Esc}[1m")
set(ColourReset "${Esc}[m")
message ("${ColourBold}Compiler Detection${ColourReset}")

project (lambda CXX)

# ----------------------------------------------------------------------------
# Begin of dependency detection
# ----------------------------------------------------------------------------

message ("\n${ColourBold}Dependency detection${ColourReset}")

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/include/seqan/version.h")
    set (CMAKE_INCLUDE_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/include
        ${CMAKE_INCLUDE_PATH})
    set (CMAKE_MODULE_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/util/cmake
        ${CMAKE_MODULE_PATH})
    message (STATUS "Found a local SeqAn library provided with the Lambda source code.")
    message (    "   This will be preferred over system global headers.")
endif ()

# ----------------------------------------------------------------------------
# Add Lambda targets
# ----------------------------------------------------------------------------

add_subdirectory(src)

# ----------------------------------------------------------------------------
# Add Tests
# ----------------------------------------------------------------------------

message ("\n${ColourBold}Setting up unit tests${ColourReset}")
add_subdirectory(tests)
