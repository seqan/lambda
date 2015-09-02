**CHANGELOG**
=============

Version 0.9.0 [2015/09/04]
--------------------------

bug-fixes and minor changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  error in man-page (#5)
-  crash in TBlastX mode (#10)
-  erroneously report unexpected extensions (#9)
-  parameter for number of matches not working (#13)
-  crash when putative duplicates heuristic is turned off (#22)
-  many small improvements to BlastIO

changes in command line interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  hide most options by default (visible again with ``--full-help``)

new features
~~~~~~~~~~~~

- ported to SeqAn 2.0 bringing in lots of smaller changes (#1)
- gzipped and bzipped input and output files supported (#19)
- some checks are performed on input data to detect wrong alphabets (#21)


Version 0.4.7 [2014/12/12]
--------------------------

bug-fixes and minor changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  fix build on Darwin / MacOS X

availability
~~~~~~~~~~~~

-  `Linux Binaries
   64Bit <http://www.seqan.de/wp-content/plugins/download-monitor/download.php?id=63>`__
-  `Linux Binaries
   Sandybridge <http://www.seqan.de/wp-content/plugins/download-monitor/download.php?id=64>`__
-  `Darwin Binaries
   64Bit <http://www.seqan.de/wp-content/plugins/download-monitor/download.php?id=65>`__
-  `Source
   code <http://www.seqan.de/wp-content/plugins/download-monitor/download.php?id=66>`__

Version 0.4.5 [2014/12/05]
--------------------------

bug-fixes and minor changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  ``lambda_indexer`` now has a different suffix array construction
   algorithm
-  this works on larger files (where the old algorithm sometimes failed)
   and is fully parallelized
-  there is also a rough progress indication when indexing

availability
~~~~~~~~~~~~

-  `Linux Binaries
   64Bit <http://www.seqan.de/wp-content/plugins/download-monitor/download.php?id=61>`__
-  `Linux Binaries
   Sandybridge <http://www.seqan.de/wp-content/plugins/download-monitor/download.php?id=60>`__
-  `Source
   code <https://github.com/h-2/seqan/releases/tag/lambda-v0.4.5>`__

Version 0.4.1 [2014/11/10]
--------------------------

bug-fixes and minor changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  default index type was not set to FM in ``lambda_indexer``

availability
~~~~~~~~~~~~

-  `Linux Binaries
   64Bit <http://www.seqan.de/wp-content/plugins/download-monitor/download.php?id=57>`__
-  `Linux Binaries
   Sandybridge <http://www.seqan.de/wp-content/plugins/download-monitor/download.php?id=58>`__
-  `Source
   code <https://github.com/h-2/seqan/releases/tag/lambda-v0.4.1>`__

Version 0.4.0 [2014/11/07]
--------------------------

performance changes in comparison to published version (0.2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  new default mode with 30-80% speed gains and up to 75% memory
   reduction over published version
-  double-indexing mode with speed gains > 100%
-  sensitivity slightly increased at the same time (1-2%)

changes in command line interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  renamed many parameters and changed some defaults
-  please look at ``lambda --help`` to see all the changes!!
-  better control of verbosity with ``-v`` parameter
-  threads now controlled with ``-t`` instead of environment variable

new features
~~~~~~~~~~~~

-  BlastN mode now usable again and proper parameter-handling added for
   it
-  added percent identity cutoff in addition to e-value cutoff (``-id``)
-  added a limit for maximum number of matches per query sequence
   (``-nm``)
-  added abundancy heuristic (``-pa``) and priorization of hits to not
   look at all hits if number of hits >> chosen limit
-  single-indexing mode which has huge memory advantages (``-qi none``)
   [now default]
-  FM-Index is now also default
-  removed Lambda-Alphabets, since they currently provide little benefit
   over Murphy10

bug-fixes and minor changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  indeces with different settings (index type, alphabet) can now be
   created on the same fasta file without conflicts between them
-  changed pre-scoring heuristic to include region around match (``-ps``
   and ``-pt``)
-  fixed build issues with gcc-4.8.x
-  FastQ support fixed

availability
~~~~~~~~~~~~

-  `Linux Binaries
   64Bit <http://www.seqan.de/wp-content/plugins/download-monitor/download.php?id=55>`__
-  `Linux Binaries
   Sandybridge <http://www.seqan.de/wp-content/plugins/download-monitor/download.php?id=54>`__
-  `Source code
   (tarball) <http://www.seqan.de/wp-content/plugins/download-monitor/download.php?id=56>`__
-  `Source code
   (git) <https://github.com/h-2/seqan/releases/tag/lambda-v0.4.0>`__

Version 0.3 [2014/09/06]
------------------------

performance changes in comparison to published version (0.2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Speed increased by ~20%
-  Suffix-Array index memory consumption reduced from 16x to 6x input
   database size
-  experimental support for FM-index as index (instead of SA) [not
   widely tested, yet]

bug-fixes and minor changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  small bugs in BLAST output formats corrected

availability
~~~~~~~~~~~~

-  `source-code
   .tar.gz <http://www.seqan.de/wp-content/plugins/download-monitor/download.php?id=53>`__
-  `source-code git <https://github.com/h-2/seqan.git>`__ commit
   d41b4b58749282dbca838a7f8506c0b378767b1b)

Version 0.2 [2014/04/07] *published version*
--------------------------------------------

-  multiple optimizations
-  added option to partition the query sequences
-  added overlapping seeds capability

availability
~~~~~~~~~~~~

-  `source-code
   .tar.gz <http://www.seqan.de/wp-content/plugins/download-monitor/download.php?id=48>`__
-  `source-code git <https://github.com/h-2/seqan.git>`__ commit
   b8ca36432d0530dd5d39560f8e2dc2cffb7c5d9d)

Version 0.1 [2014/01/15] *initial release*
------------------------------------------