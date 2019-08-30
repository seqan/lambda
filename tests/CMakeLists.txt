# ===========================================================================
#  Lambda tests
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)

enable_testing ()
include (CTest)

## only subset of tests
if (LAMBDA_FASTBUILD)
    set (PROGS blastp blastx)
else ()
    set (PROGS blastn blastp blastx tblastn tblastx)
endif ()

## basic indexer tests
foreach(PROG ${PROGS})
    add_test (NAME test_mkindex_${PROG}
              COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/maintests.sh
                      "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" ${PROG} ${DI} "MKINDEX" " ")
endforeach()

## basic search tests
foreach(PROG ${PROGS})
    foreach(FF m0 m8 m9 sam bam m9.gz sam.bz2)
        add_test (NAME test_search_${PROG}_${FF}
                  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/maintests.sh
                          "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" ${PROG} ${DI} "SEARCH" ${FF})
    endforeach()
endforeach()
