#!/bin/sh

errorout()
{
    echo $1 #> /dev/stderr
    [ "$MYTMP" = "" ] || rm -r "${MYTMP}"
    exit 1
}

[ $# -ne 7 ] && exit 1

SRCDIR=$1
BINDIR=$2
PROG=$3
DI=$4
MODE=$5
TOOL=$6
EXTENSION=$7

# check existence of commands
which openssl gunzip mktemp diff cat zcat zgrep > /dev/null
[ $? -eq 0 ] || errorout "Not all required programs found. Needs: openssl gunzip mktemp diff cat zcat zgrep"

SALPH=prot      # actual subject alph
QALPHIN=prot    # query input file alph
SALPHIN=prot    # subject input file alph

MKINDEX=mkindexp # protein mode by default
SEARCH=searchp   # protein mode by default

case "$PROG" in "blastn")
    QALPHIN=nucl
    SALPH=nucl
    SALPHIN=nucl
    MKINDEX=mkindexn
    SEARCH=searchn
    ;;
"blastp")
    ;;
"blastx")
    QALPHIN=nucl
    ;;
"tblastn")
    SALPH=trans
    SALPHIN=nucl
    ;;
"tblastx")
    SALPH=trans
    QALPHIN=nucl
    SALPHIN=nucl
    ;;
esac

MYTMP="$(mktemp -q -d -t "$(basename "$0").XXXXXX" 2>/dev/null || mktemp -q -d)"
[ $? -eq 0 ] || errorout "Could not create tmp"

cd "$MYTMP"
[ $? -eq 0 ] || errorout "Could not cd to tmp"

gunzip < "${SRCDIR}/tests/db_${SALPHIN}.fasta.gz" > db.fasta
[ $? -eq 0 ] || errorout "Could not unzip database file"

${BINDIR}/bin/lambda3 ${MKINDEX} -d db.fasta -i db_${SALPH}_${DI}.fasta.gz.lba --db-index-type ${DI}
[ $? -eq 0 ] || errorout "Could not run the indexer"

[ "$(openssl md5 db_${SALPH}_${DI}.fasta.gz.lba)" = \
"$(zgrep "(db_${SALPH}_${DI}.fasta.gz.lba)" "${SRCDIR}/tests/index_test_outfile.md5sums.gz")" ] || errorout "MD5 mismatch of index file"

# INDEXER tests end here
if [ "$TOOL" = "MKINDEX" ]; then
    rm -r "${MYTMP}"
    exit 0
fi

gunzip < "${SRCDIR}/tests/queries_${QALPHIN}.fasta.gz" > queries.fasta
[ $? -eq 0 ] || errorout "Could not unzip queries.fasta"

${BINDIR}/bin/lambda3 ${SEARCH} -i db_${SALPH}_${DI}.fasta.gz.lba -q queries.fasta -t 1 --version-to-outputfile 0 --seed-offset 7 --seed-length 14 --adaptive-seeding 1 -m $MODE \
-o output_${PROG}_${DI}_${MODE}.${EXTENSION}
[ $? -eq 0 ] || errorout "Search failed."

[ "$(openssl md5 output_${PROG}_${DI}_${MODE}.${EXTENSION})" = \
"$(zgrep "(output_${PROG}_${DI}_${MODE}.${EXTENSION})" "${SRCDIR}/tests/search_test_outfile.md5sums.gz")" ] || errorout "MD5 mismatch of output file"

rm -r "${MYTMP}"
