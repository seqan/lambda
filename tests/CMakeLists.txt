

cmake_minimum_required (VERSION 2.8.2)

message (STATUS "Configuring tests")

foreach(PROG blastn blastp blastx tblastn tblastx)
    foreach(DI sa fm)
        add_test (NAME lambda_indexer_test_${PROG}_${DI}
                  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/mkindex.sh
                          "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" ${PROG} ${DI})
    endforeach()
endforeach()

## TODO remove files afterwards, remove folder hierarchy
## only .fasta.gz and md5sums files are required