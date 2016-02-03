
## basic indexer tests
foreach(PROG blastn blastp blastx tblastn tblastx)
    foreach(DI sa fm)
        add_test (NAME test_mkindex_${PROG}_${DI}
                  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/maintests.sh
                          "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" ${PROG} ${DI} "MKINDEX" " ")
    endforeach()
endforeach()

## basic search tests
foreach(PROG blastn blastp blastx tblastn tblastx)
    foreach(DI sa fm)
        foreach(FF m0 m8 m9 sam bam m9.gz sam.bz2)
            add_test (NAME test_search_${PROG}_${DI}_${FF}
                      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/maintests.sh
                             "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" ${PROG} ${DI} "SEARCH" ${FF})
        endforeach()
    endforeach()
endforeach()

## TODO remove files afterwards, remove folder hierarchy
## only .fasta.gz and md5sums files are required