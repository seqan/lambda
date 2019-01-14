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
MODE=$5
EXTENSION=$6

# check existence of commands
which openssl gunzip mktemp diff cat zcat zgrep > /dev/null
[ $? -eq 0 ] || errorout "Not all required programs found. Needs: openssl gunzip mktemp diff cat zcat zgrep"

SALPH=prot      # actual subject alph
QALPHIN=prot    # query input file alph
SALPHIN=prot    # subject input file alph
INDEXER=mkindexp
SEARCHER=searchp
case "$PROG" in "blastn")
    INDEXER=mkindexn
    SEARCHER=searchn
    QALPHIN=nucl
    SALPH=nucl
    SALPHIN=nucl
    ;;
"blastp")
    ;;
"blastx")
    QALPHIN=nucl
    ;;
"tblastn")
    INDEXER=mkindexn
    SALPH=trans
    SALPHIN=nucl
    ;;
"tblastx")
    INDEXER=mkindexn
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

${BINDIR}/bin/lambda2 ${INDEXER} -d db.fasta -i ${DI}.lambda
[ $? -eq 0 ] || errorout "Could not run the indexer"

openssl md5 $(find * -type f) > md5sums
[ $? -eq 0 ] || errorout "Could not run md5 or md5sums"

cp "${SRCDIR}/tests/db_${SALPH}_${DI}.md5sums" md5sums.orig
[ $? -eq 0 ] || errorout "Could not copy md5sums.orig"

[ "$(cat md5sums)" = "$(cat md5sums.orig)" ] || errorout "$(diff -u md5sums md5sums.orig)"

## INDEXER tests end here
if [ "$MODE" = "MKINDEX" ]; then
    rm -r "${MYTMP}"
    exit 0
fi

gunzip < "${SRCDIR}/tests/queries_${QALPHIN}.fasta.gz" > queries.fasta
[ $? -eq 0 ] || errorout "Could not unzip queries.fasta"

${BINDIR}/bin/lambda2 ${SEARCHER} -i ${DI}.lambda -q queries.fasta -t 1 --version-to-outputfile off \
-o output_${PROG}_${DI}.${EXTENSION}
[ $? -eq 0 ] || errorout "Search failed."

[ "$(openssl md5 output_${PROG}_${DI}.${EXTENSION})" = \
"$(zgrep "(output_${PROG}_${DI}.${EXTENSION})" "${SRCDIR}/tests/search_test_outfile.md5sums.gz")" ] || errorout "MD5 mismatch of output file"

rm -r "${MYTMP}"
