#!/bin/sh

errorout()
{
    echo $1 #> /dev/stderr
    [ "$MYTMP" = "" ] || rm -r "${MYTMP}"
    exit 1
}

[ $# -ne 6 ] && exit 1

SRCDIR=$1
BINDIR=$2
PROG=$3
DI=$4
TOOL=$5
EXTENSION=$6

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
"blastn_bs")
    QALPHIN=nucl_bs
    SALPH=nucl_bs
    SALPHIN=nucl_bs
    MKINDEX=mkindexn
    SEARCH=searchn
    ;;
esac

MYTMP="$(mktemp -q -d -t "$(basename "$0").XXXXXX" 2>/dev/null || mktemp -q -d)"
[ $? -eq 0 ] || errorout "Could not create tmp"

cd "$MYTMP"
[ $? -eq 0 ] || errorout "Could not cd to tmp"

gunzip < "${SRCDIR}/tests/db_${SALPHIN}.fasta.gz" > db.fasta
[ $? -eq 0 ] || errorout "Could not unzip database file"

if [ "$PROG" = "blastn_bs" ]; then
    ${BINDIR}/bin/lambda3 ${MKINDEX} -d db.fasta -i db_${SALPH}_${DI}.fasta.gz.lba --db-index-type ${DI} -r dna3bs
else
    ${BINDIR}/bin/lambda3 ${MKINDEX} -d db.fasta -i db_${SALPH}_${DI}.fasta.gz.lba --db-index-type ${DI}
fi
[ $? -eq 0 ] || errorout "Could not run the indexer"

# INDEXER tests end here
if [ "$TOOL" = "MKINDEX" ]; then
    [ "$(openssl md5 db_${SALPH}_${DI}.fasta.gz.lba)" = \
    "$(zgrep "(db_${SALPH}_${DI}.fasta.gz.lba)" "${SRCDIR}/tests/index_test_outfile.md5sums.gz")" ] || errorout "MD5 mismatch of index file"
    rm -r "${MYTMP}"
    exit 0
fi

gunzip < "${SRCDIR}/tests/queries_${QALPHIN}.fasta.gz" > queries.fasta
[ $? -eq 0 ] || errorout "Could not unzip queries.fasta"

${BINDIR}/bin/lambda3 ${SEARCH} -i db_${SALPH}_${DI}.fasta.gz.lba -q queries.fasta -t 2 --version-to-outputfile 0 \
-o output_${PROG}_${DI}_fullSIMD.${EXTENSION}
[ $? -eq 0 ] || errorout "Search failed."

[ "$(openssl md5 output_${PROG}_${DI}_fullSIMD.${EXTENSION})" = \
"$(zgrep "(output_${PROG}_${DI}_fullSIMD.${EXTENSION})" "${SRCDIR}/tests/search_test_outfile.md5sums.gz")" ] || errorout "MD5 mismatch of output file"

rm -r "${MYTMP}"
