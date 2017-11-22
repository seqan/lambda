Lambda: the Local Aligner for Massive Biological DatA
-----------------------------------------------------

Lambda is a local aligner optimized for many query sequences and searches in protein space. It...

* is highly compatible to BLAST (bitscore and e-value statistics, tab seperated and verbose output formats)
* is much faster than BLAST and many other comparable tools
* supports many other input and output formats, including standards-conformant ``.sam`` and ``.bam`` and many compression types
* has special features for species annotation and taxonomic analysis
* is well-documented and easy to use (e.g. provides progress-bars and memory usage estimates)

downloads and installation
--------------------------

+------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
|  **Executables**                                                                                                                                                                                                                      |
+------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
| .. image:: https://raw.githubusercontent.com/seqan/lambda/gh-pages/images_readme/appbar.disk.download.png        | Pre-built executables for GNU/Linux, Mac and FreeBSD are available from the                                        |
|    :alt: Download Executables                                                                                    | `releases page <https://github.com/seqan/lambda/releases>`__.                                                      |
|    :target: https://github.com/seqan/lambda/releases                                                             |                                                                                                                    |
|    :width: 76px                                                                                                  |                                                                                                                    |
+------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
|  **Source code**                                                                                                                                                                                                                      |
+------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
| .. image:: https://raw.githubusercontent.com/seqan/lambda/gh-pages/images_readme/appbar.column.three.png         | You can also build lambda from source which will result in binaries optimized for your                             |
|    :alt: Build from source                                                                                       | specific system (and thus faster). For instructions, please see the                                                |
|    :target: https://github.com/seqan/lambda/wiki                                                                 | `wiki <https://github.com/seqan/lambda/wiki>`__.                                                                   |
|    :width: 76px                                                                                                  |                                                                                                                    |
+------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+

usage instructions
------------------


Before you can search, you need to have an index. You can

1. download and unzip a pre-built index from the `wiki <https://github.com/seqan/lambda/wiki>`__; or
2. index one yourself (this can take some time but only has to be done once):

::

    % bin/lambda2 mkindexp -d db.fasta

*(in case you want to create a nucleotide index, instead use ``mkindexn``)*

After that running Lambda is as simple as

::

    % bin/lambda2 searchp -q query.fasta -i db.fasta.lambda

*(in case you want to perform a nucleotide search, instead use ``searchn``)*

For a list of options, see the help pages:

::

    % bin/lambda2 --help
    % bin/lambda2 COMMAND --help

Advanced options are available via ``--full-help`` or the man pages, and more documentation is available
in the `wiki <https://github.com/seqan/lambda/wiki>`__.

authorship and copyright
------------------------

Lambda is being developed by `Hannes Hauswedell <mailto:hannes.hauswedell@fu-berlin.de>`__, but it incorporates a lot of work from other members of the `SeqAn project <http://www.seqan.de>`__.

+------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
|  **Please always cite the publication, also if using Lambda in comparisons and pipelines**                                                                                                                                            |
+------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
| .. image:: https://raw.githubusercontent.com/seqan/lambda/gh-pages/images_readme/appbar.book.hardcover.open.png  | *Lambda: the local aligner for massive biological data*;                                                           |
|    :alt: Please cite                                                                                             | Hannes Hauswedell, Jochen Singer, Knut Reinert;                                                                    |
|    :target: http://bioinformatics.oxfordjournals.org/content/30/17/i349.abstract                                 | `Bioinformatics 2014 30 (17): i349-i355 <http://bioinformatics.oxfordjournals.org/content/30/17/i349.abstract>`__; |
|    :width: 76px                                                                                                  | doi: 10.1093/bioinformatics/btu439                                                                                 |
+------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
| **Please respect the license of the software**                                                                                                                                                                                        |
+------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
| .. image:: https://raw.githubusercontent.com/seqan/lambda/gh-pages/images_readme/copyleft.png                    | Lambda is Free and open source software, so you can use it for any purpose, free of charge.                        |
|    :alt: Respect the license                                                                                     | However certain conditions apply when you (re-)distribute and/or modify Lambda, please respect the                 |
|    :target: https://github.com/seqan/lambda/blob/master/LICENSE.rst                                              | `license <https://github.com/seqan/lambda/blob/master/LICENSE.rst>`__.                                             |
|    :width: 76px                                                                                                  |                                                                                                                    |
+------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+

feedback & updates
------------------

+-------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
| .. image:: https://raw.githubusercontent.com/seqan/lambda/gh-pages/images_readme/appbar.social.github.octocat.png | You can ask questions and report bugs on the `github tracker <https://github.com/seqan/lambda/issues>`__ .         |
|    :alt: GitHub                                                                                                   | Please also `subscribe <https://github.com/seqan/lambda/subscription>`__ and/or star us!                           |
|    :target: https://github.com/seqan/lambda/issues                                                                |                                                                                                                    |
|    :width: 76px                                                                                                   |                                                                                                                    |
+-------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
| .. image:: https://raw.githubusercontent.com/seqan/lambda/gh-pages/images_readme/appbar.email.png                 | To stay up to date via e-mail, please subscribe to the                                                             |
|    :alt: Newsletter                                                                                               | `newsletter <https://lists.fu-berlin.de/listinfo/lambda-users>`__. There is on average less than one e-mail        |
|    :target: https://lists.fu-berlin.de/listinfo/lambda-users                                                      | per month.                                                                                                         |
|    :width: 76px                                                                                                   |                                                                                                                    |
+-------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
| .. image:: https://raw.githubusercontent.com/seqan/lambda/gh-pages/images_readme/appbar.social.twitter.png        | You can also follow SeqAn on `twitter <https://twitter.com/SeqAnLib>`__ to receive updates on Lambda.              |
|    :alt: Newsletter                                                                                               |                                                                                                                    |
|    :target: https://twitter.com/SeqAnLib                                                                          |                                                                                                                    |
|    :width: 76px                                                                                                   |                                                                                                                    |
+-------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+

*icons on this page by Austin Andrews / https://github.com/Templarian/WindowsIcons*
