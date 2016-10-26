#
# Aarti Jajoo, Evan Biederstedt, Oct 2016
#

"""
This is still very pseudo-code-ish, as we need to definitively define how to work with the field in the "field_line"

See how we work with array `score[i]`, as well as defining `line==1`
"""

import numpy as np

def calculateScore(field_line, match_score, mismatch_score, gap_score):
    # parameters: 
    # field_line --- some file format which inputs fields we need for the JGA; currently thinking individual text files
    # Penalty scores----naturally users can override if they need to, etc.
    # match_score = 1 (default)
    # mismatch_score = -1 (default)
    # gap_score = -1 (default)
    
    for line in field_line:
        # from field line, we get: n_matches, n_mismatches, n_gaps---GET THOSE
        # IDEA: append to each "key" (i.e. each unique 'subgroup' identifier)
        score[i] = n_matches * match_score + n_mismatches * mismatch_score + n_gaps * gap_score
    

    if line == 1: # that is, there's only one diagonal
        #
        # -(-ContigStart + firstPos_read) * gap_score  // this is a gap, start
        #  -(ContigEnd - lastPos_read) * gap_score     // this is a gap, end
        #  score                                       // this is 1-D element
        #
        GlobalScore = (-(-ContigStart + firstPos_read) * gap_score -
                       (ContigEnd - lastPos_read) * gap_score + score)
    else:
        JGA()  # JGA(startContig, lengthContig, start, end, score




