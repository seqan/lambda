## Lambda: the Local Aligner for Massive Biological DatA

Lambda is a BLAST compatible local aligner optimized for NGS, metagenomics
and protein space.

For more information, see
 * the [homepage](https://www.seqan.de/projects/lambda/)
 * the [publication](http://bioinformatics.oxfordjournals.org/content/30/17/i349.abstract)
 * or write to [Hannes Hauswedell](mailto:hannes.hauswedell@[molgen.mpg.de|fu-berlin.de])

## Download

The latest release is available from the
[homepage](https://www.seqan.de/projects/lambda/).

Latest git version is available from
[github](https://github.com/h-2/seqan/tree/feature/lambda/extras/apps/lambda).

Lambda is Free Software: you may use it for any purpose and you can
redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
See the LICENSE file for more information.

## Build

    % tar xzf seqan-lambda-v0.4.tar.gz
    % mkdir -p seqan-lambda-build/release
    % cd seqan-lambda-build/release
    % cmake ../../seqan-lambda-v0.4 \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_FLAGS:STRING="-march=native"
    % make -j2 lambda lambda_indexer

GCC >= 4.9.1 is the recommended compiler. Warnings concerning a lack of C++14 can be
ignored. Please be aware that due to excessive use of templating and
compile-time optimizations the build might take well over 10min.


## Run

Optionally mask the database:

    % /path/to/segmasker -infmt fasta -in db.fasta -outfmt interval -out db.seg

Run the indexer:

    % bin/lambda_indexer -d db.fasta [-s db.seg]

Run lambda:

    % bin/lambda -q query.fasta -d db.fasta

For a full list of options, see the help page:

    % bin/lambda --help

## Give feedback

Please report bugs to the SeqAn
[bug-tracker](https://github.com/seqan/seqan).

Any other questions or feedback you can send to
[Hannes Hauswedell](mailto:hannes.hauswedell@[molgen.mpg.de|fu-berlin.de]).

Thanks for using Lambda, we hope that it is useful to you! If yes, please cite
us :)
