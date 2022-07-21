# ===========================================================================
#  Lambda tests
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)

enable_testing ()
include (CTest)

## only subset of tests
set (PROGS blastn blastp blastx tblastn tblastx blastn_bs)

## basic indexer tests
foreach(PROG ${PROGS})
    foreach(DI fm bifm)
        add_test (NAME test_mkindex_${PROG}_${DI}
                COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/maintests.sh
            "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" ${PROG} ${DI} "MKINDEX" " ")
    endforeach()
endforeach()

## basic search tests
foreach(PROG ${PROGS})
    foreach(DI fm bifm)
        if(NOT CMAKE_SYSTEM_NAME STREQUAL Darwin)
            foreach(FF m0 m8 m9 sam bam m9.gz sam.bz2)
                add_test (NAME test_search_${PROG}_${DI}_fullSIMD_${FF}
                          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/maintests.sh
                          "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" ${PROG} ${DI} "SEARCH" ${FF})
            endforeach()
        else()
            foreach(FF m0 m8 m9 sam bam sam.bz2)
                add_test (NAME test_search_${PROG}_${DI}_fullSIMD_${FF}
                          COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/maintests.sh
                          "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" ${PROG} ${DI} "SEARCH" ${FF})
            endforeach()
        endif()
    endforeach()
endforeach()
