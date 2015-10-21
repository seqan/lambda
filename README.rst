Lambda: the Local Aligner for Massive Biological DatA
-----------------------------------------------------

Lambda is a local aligner optimized for many query sequences and searches in protein space.
It is compatible to BLAST, but much faster than BLAST and many other comparable tools.

cite
----------

Please cite the following if you use Lambda anywhere in your academic work, also as part of pipelines
or comparisons:

*Lambda: the local aligner for massive biological data*;
Hannes Hauswedell, Jochen Singer, Knut Reinert;
`Bioinformatics 2014 30 (17): i349-i355 <http://bioinformatics.oxfordjournals.org/content/30/17/i349.abstract>`__;
doi: 10.1093/bioinformatics/btu439

download
--------

The latest version is available 
`here <https://github.com/seqan/lambda/releases>`__. Versions prior to 0.9.0 are available 
`here <https://github.com/h-2/seqan/releases>`__.

Lambda is Free and open source software, so you can use it for any purpose, free of charge.
However certain conditions apply when you (re-)distribute or modify Lambda, please respect the
`license <./COPYING.rst>`__.

You can also build lambda from source which will result in binaries optimized for your
specific system (and thus faster). For instructions, please see the
`wiki <https://github.com/seqan/lambda/wiki>`__.

run
---

Optionally mask the database:

::

    % /path/to/segmasker -infmt fasta -in db.fasta -outfmt interval -out db.seg

Run the indexer (or check the `wiki <https://github.com/seqan/lambda/wiki>`__ for pre-built indexes!):

::

    % bin/lambda_indexer -d db.fasta [-s db.seg]

Run lambda:

::

    % bin/lambda -q query.fasta -d db.fasta

*Please note that if you downloaded the binaries from the web-site, you might also have* ``bin/lambda-avx`` *which is
measurably faster, so try that first!*

For a list of options, see the help pages:

::

    % bin/lambda --help
    % bin/lambda --full-help

Or visit the Tuning-guide in the `wiki <https://github.com/seqan/lambda/wiki>`__.

give feedback
-------------

Please report bugs to the `bug-tracker <https://github.com/seqan/lambda/issues>`__.

Any other questions or feedback you can send to 
`Hannes Hauswedell <mailto:hannes.hauswedell@[molgen.mpg.de|fu-berlin.de]>`__.

Thank you for using Lambda, we hope that it is useful to you!