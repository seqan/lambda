# ===========================================================================
#  Lambda tests
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)

enable_testing ()
include (CTest)

## only subset of tests
set (PROGS blastn blastp blastx tblastn tblastx)
set (MODES fullSerial fullSIMD)

## basic indexer tests
foreach(PROG ${PROGS})
    foreach(DI fm bifm)
        add_test (NAME test_mkindex_${PROG}_${DI}
                COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/maintests.sh
            "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" ${PROG} ${DI} " " "MKINDEX" " ")
    endforeach()
endforeach()

## basic search tests
foreach(PROG ${PROGS})
    foreach(DI fm bifm)
        foreach(MODE ${MODES})
            if(NOT CMAKE_SYSTEM_NAME STREQUAL Darwin)
                foreach(FF m0 m8 m9 sam bam m9.gz sam.bz2)
                    add_test (NAME test_search_${PROG}_${DI}_${MODE}_${FF}
                              COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/maintests.sh
                              "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" ${PROG} ${DI} ${MODE} "SEARCH" ${FF})
                endforeach()
            else()
                foreach(FF m0 m8 m9 sam bam sam.bz2)
                    add_test (NAME test_search_${PROG}_${DI}_${MODE}_${FF}
                              COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/maintests.sh
                              "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" ${PROG} ${DI} ${MODE} "SEARCH" ${FF})
                endforeach()
            endif()
        endforeach()
    endforeach()
endforeach()

## Add separate bisulfite tests
foreach(DI fm bifm)
    add_test (NAME test_mkindex_blastn_bs_${DI}
              COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/maintests.sh
          "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" blastn_bs ${DI} " " "MKINDEX" " ")
endforeach()

foreach(DI fm bifm)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL Darwin)
        foreach(FF m0 m8 m9 sam bam m9.gz sam.bz2)
            add_test (NAME test_search_blastn_bs_${DI}_fullSerial_${FF}
                      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/maintests.sh
                      "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" blastn_bs ${DI} fullSerial "SEARCH" ${FF})
        endforeach()
    else()
        foreach(FF m0 m8 m9 sam bam sam.bz2)
            add_test (NAME test_search_blastn_bs_${DI}_fullSerial_${FF}
                      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/maintests.sh
                      "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" blastn_bs ${DI} fullSerial "SEARCH" ${FF})
        endforeach()
    endif()
endforeach()
