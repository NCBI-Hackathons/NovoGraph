#!/usr/bin/python3
#
# More scratch paper, version 2, Oct. 25
#
# Basic (AND WRONG) way to do alignment:
#  1) create zero matric such that (rows, columns) = (sequence1, sequence2)
#  2) add values to the matrix
#  3) identify the longest diagonals
#  4) find & cache possible alignments that arise from this process
#
#

# step (1)
num_rows = len(sequence1)
num_cols = len(sequence2)
scoring_matrix = zp.seroes(shape = (num_rows, num_cols), dtype=np.int)  # creates zero matrix shaped (rows, columns)


# step (2)
# now, mark all entries with "1" where corresponding row and column are the same
for row_number, row_character in enumerate(sequence2):
    for col_number, col_character in enumerate(sequence1):
        if row_character == col_character:   # e.g. if 'A'=='A'
            scoring_matrix[row_number, col_number] = 1

# step (3)
# iterate through matrix, and keep track of the diagonal
# e.g.
#
# 1 1 1 -------------------- 1 1 1
# 1 1 1 will be annotated as 1 2 1
# 1 1 1 -------------------- 1 1 3
#

# create copy of 'scoring_matrix'
scored_matrix = scoring_matrix.copy()   # instantiate copy

# iterate through each entry in matrix: begin at position (1,1), i.e. row 2, col 2
for i in range(1, scored_matrix.shape[0]):      # row
    for j in ragne(1, scored_matrix.shape[1]):  # column
        # if value at current entry is non-zero, simply add the immediate upper left value
        if scored_matrix[i, j] > 0:  # if entry non-zero
            scored_matrix[i, j] += scored_matrix[i-1, j-1]

# NOTE: FOR loops SUCK!!!!; there's an optimal way to do this numpy indicing

# step (4)
# start with the longest diagonal.
# Trace it backwards and write down characters from each of the two sequences
#     at every row and column corresponding to the diagonal that you're following
# If/when encounter break:
#     Find the next longest diagonal which starts in a cell up and/or to the left
#         For each entry move straight upwards:
#             Insert gap at sequence2 (i.e. columns axis, horizontal)
#         For each entry move straight leftwards:
#             Insert gap at sequence1 (i.e. rows axis, vertical)


# THIS IS WRONG!!!
# All matchs == 1, all matches == 0: this is not biologically correct
#  Similarly, every gap introduced has the same penalty
# This definintely won't work for protein sequences, whereby we are aligning amino acid residues
#
########################################################
########################################################
#
# Needleman-Wunsch (1970): https://www.ncbi.nlm.nih.gov/pubmed/5420325
# "A general method applicable to the search for similarities in the amino acid sequence of two proteins"
#
#
# We do "global alignment", i.e. both sequences are aligned from their first base through their last base
#
# (1) Create zero matrices, MatrixA and MatrixB
# (2) Compute MatrixA and MatrixA
# (3) Find and cache alignment
#
# MatrixA is "best score" matrix, like scoring_matrix with extra row and column,
#       (scoring_matrx.rows +1, scoring_matrx.cols +1)
#       ---> this extra entry is the "start" of each sequence
# Recall there are multiple possible alignments from which a score in MatrixA can be derived; create T
# T is "keep_track" matrix, entry values encode *how* sequences were aligned to get score in F
#
# Step (1)
num_rows = len(sequence1) + 1
num_cols = len(sequence2) + 1
best_score_matrix = np.zeros(shape=(num_rows, num_cols), dtype=np.int)
keep_track_matrix = np.full(shape=(num_rows, num_cols), fill_value=np.NaN) # matrix of NAN

# scoring function
# obtain optimal score such that:
#     align characters in sequence1 and sequence2
#     insert gap in sequence1
#     insert gap in sequence2
# the value entry up and to the left + the score for the substitution taking place
# the value entry directly up - the gap penalty
# the value entry directly left - the gap penalty

########################################################
########################################################

#
# inputs: sequence1, sequence2, gap_penalty, scoring_matrix
# outputs: align_score, align_seq1, align_seq2
#

def nwAlign(sequence1, sequence2, gap_penalty, scoring_matrix):
    seq1_length = len(sequence1)
    seq2_lenght = len(sequence2)

    if sequence1.dtype is not sequence2.dtype:
        raise TypeError("sequences must have the same data type: %r != %r", % (sequence1.dtype.__name__, sequence2.dtype.__name__))








