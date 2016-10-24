"""
------------------------------------------------------------------------------
nwalign: fast `cython`_  - `Needleman-Wunsch`_ alignment
------------------------------------------------------------------------------

.. _`Needleman-Wunsch`: http://en.wikipedia.org/wiki/Needleman-Wunsch_algorithm 
.. _`scoring matrix`: http://en.wikipedia.org/wiki/Substitution_matrix
.. _`cython`: http://cython.org

This module provides a python module and a command-line interface to do global-
sequence alignment using the `Needleman-Wunsch` algorithm. It uses `cython`_ 
and numpy for speed.

Command-Line Usage 
==================
the nwalign executable is installed to the PATH by setuptools
::

    $ nwalign alphabet alpet
    alphabet
    alp---et

specify an alignment `scoring matrix`_ 
::

    $ nwalign --matrix /usr/share/ncbi/data/BLOSUM62 EEAEE EEEEG
    EEAEE-
    EE-EEG

with specified penalties
::

    $ nwalign --gap_open -10 --gap_extend -4 --match 12 ASDFF ASFF
    ASDFF
    AS-FF


Python Usage
============
Alignment
---------
::

    >>> import nwalign as nw
    >>> nw.global_align("CEELECANTH", "PELICAN", matrix='PAM250')
    ('CEELE-CANTH', '-PEL-ICAN--')

    # with a specified penalty for open and extend.
    >>> nw.global_align("CEELECANTH", "PELICAN", gap_open=-10, gap_extend=-4, matrix='PAM250')
    ('CEELECANTH', '-PELICAN--')


the `matrix` is specified as the full path to an `scoring matrix`_ as
is distributed with the NCBI toolset.

Scoring
-------
get the score of an alignment. (the first 2 args are from an alignment
and must have the same length.
::

    >>> nw.score_alignment('CEELECANTH', '-PELICAN--', gap_open=-5,
    ...                     gap_extend=-2, matrix='PAM250')
    11

    >>> nw.score_alignment('CEELE-CANTH', '-PEL-ICAN--', gap_open=-5,
    ...                     gap_extend=-2, matrix='PAM250')
    6


"""
from cnwalign import global_align, global_align_no_matrix, score_alignment


def main():
    import sys
    import optparse
    parser = optparse.OptionParser(usage="""
    %prog [options] seq1 seq2 
    """)
    parser.add_option("--gap_extend", dest="gap_extend", help="gap extend penalty (must be integer <= 0)", type="int", default=-1)
    parser.add_option("--gap_open", dest="gap_open", help="gap open penalty (must be integer <= 0)", type="int", default=-1)
    parser.add_option("--match", dest="match", help="match score (must be integer > 0)", type="int", default=1)
    parser.add_option("--matrix", dest="matrix", help="scoring matrix in ncbi/data/ format,\
                                      if not specificied, match/mismatch are used", default=None)
    parser.add_option("--server", dest="server", default=0, type='int',
                      help="if non-zero integer, a server is started")

    try:
        options, args = parser.parse_args()
    except:
        sys.exit(parser.print_help())
    if options.server != 0:
        import nwserver
        nwserver.main(sys.argv[2])
        
    elif len(args) != 2:
        sys.exit(parser.print_help())
    else:
        print "\n".join(global_align(args[0], args[1], options.match,
                                 options.gap_open, options.gap_extend, options.matrix))

if __name__ == "__main__":
    main()
