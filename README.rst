Lambda: the Local Aligner for Massive Biological DatA
-----------------------------------------------------

Lambda is a BLAST compatible local aligner optimized for NGS,
metagenomics and protein space.

For more information, see 
 * the `homepage <https://www.seqan.de/projects/lambda/>`__
 * the `publication <http://bioinformatics.oxfordjournals.org/content/30/17/i349.abstract>`__
 * or write to `Hannes Hauswedell <mailto:hannes.hauswedell@[molgen.mpg.de|fu-berlin.de]>`__

Download
--------

The latest version is available 
`here <https://github.com/seqan/lambda/releases>`__. Versions prior to 0.9.0 are available 
`here <https://github.com/h-2/seqan/releases>`__.

Lambda is Free Software: you may use it for any purpose and you can
redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version. See the
COPYING.rst file for more information.

Build
-----

::

    % tar xzf lambda-vX.Y.Z.tar.gz
    % mkdir -p lambda-build/release
    % cd lambda-build/release
    % cmake ../../lambda-vX.Y.Z
    % make -j2

Currently only GCC >= 4.9.1 is supported. Please be aware that due to excessive use of
templating and compile-time optimizations the build might take well over
10min.

Run
---

Optionally mask the database:

::

    % /path/to/segmasker -infmt fasta -in db.fasta -outfmt interval -out db.seg

Run the indexer:

::

    % bin/lambda_indexer -d db.fasta [-s db.seg]

Run lambda:

::

    % bin/lambda -q query.fasta -d db.fasta

For a full list of options, see the help page:

::

    % bin/lambda --help

or look at MANUAL.rst

Give feedback
-------------

Please report bugs to the `bug-tracker <https://github.com/seqan/lambda/issues>`__.

Any other questions or feedback you can send to 
`Hannes Hauswedell <mailto:hannes.hauswedell@[molgen.mpg.de|fu-berlin.de]>`__.

Thanks for using Lambda, we hope that it is useful to you! If yes,
please cite us :)