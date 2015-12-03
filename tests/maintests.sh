#!/bin/sh

errorout()
{
    echo $1 #> /dev/stderr
    rm -r "${MYTMP}"
    exit 1
}

[ $# -ne 4 ] && exit 1

SRCDIR=$1
BINDIR=$2
PROG=$3
DI=$4

# check existence of commands
which openssl gunzip mktemp diff cat zcat > /dev/null
[ $? -eq 0 ] || errorout "Not all required programs found. Needs: openssl gunzip mktemp diff cat zcat"

ALPH=nucl
if [ "$PROG" = "blastp" ] || [ "$PROG" = "blastx" ]; then
    ALPH=prot
fi
ALPHIN=$ALPH

if [ "$PROG" = "tblastn" ] || [ "$PROG" = "tblastx" ]; then
    ALPH=trans
fi

MYTMP="$(mktemp -q -d -t "$(basename "$0").XXXXXX" 2>/dev/null || mktemp -q -d)"
[ $? -eq 0 ] || errorout "Could not create tmp"

cd "$MYTMP"
[ $? -eq 0 ] || errorout "Could not cd to tmp"

gunzip < "${SRCDIR}/tests/${ALPHIN}_db.fasta.gz" > db.fasta
[ $? -eq 0 ] || errorout "Could not unzip database file"

${BINDIR}/bin/lambda_indexer -d db.fasta -di ${DI} -p ${PROG}
[ $? -eq 0 ] || errorout "Could not run the indexer"

openssl md5 * > md5sums
[ $? -eq 0 ] || errorout "Could not run md5 or md5sums"

gunzip < "${SRCDIR}/tests/${ALPH}_${DI}_db.md5sums.gz" > md5sums.orig
[ $? -eq 0 ] || errorout "Could not unzip md5sums.orig"

[ "$(cat md5sums)" = "$(cat md5sums.orig)" ] || errorout "$(diff -u md5sums md5sums.orig)"

rm -r "${MYTMP}"