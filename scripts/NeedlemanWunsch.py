#!/usr/bin/python3                                                                                                                     

# First version, Needleman-Wunsch                                                               

import numpy as np

penalty = {"match":1, "mismatch":-1, "gap": -1} # create dict of penalty scores                                                        

def nwAlign(string1, string2):
    if string1 == string2:
        return penalty["match"]    # if string1 matches string2, match                                                                 
    elif string1 == "-" or string2 == "-":
        return penalty["gap"]
    else:
        return penalty["mismatch"]

def needlemanWunsch(seq1, seq2):  # (M, N) = (rows, columns)                                                                           
    M = len(seq1)  # rows                                                                                                              
    N = len(seq2)  # columns                                                                                                           
    score = np.zeros((M+1, N+1))

    for i in range(M+1):
        score[i][0] = penalty["gap"] * i
    for j in range(N+1):
        score[0][j] = penalty["gap"] * j

    for i in range(1, M+1):
        for j in range(1, N+1):
            diagonal = score[i-1][j-1] + nwAlign(seq1[i-1], seq2[j-1])
            delete =  score[i-1][j] + penalty["gap"]
            insert = score[i][j-1] + penalty["gap"]
            score[i][j] = max(diagonal, delete, insert)

    print("score matrix = \n%s\n" % score)
    align1, align2 = "", ""
    i, j = m, n

    while i > 0 and j > 0:
        score_current = score[i][j]
        score_diag = score[i-1][j-1]
        score_left = score[i][j-1]
        score_up = score[i-1][j]

        if score_current == score_diag + nwAlign(seq1[i-1], seq2[j-1]):
            print("diag")
            align1, align2 = seq1[i-1], seq2[j-1]
            i, j = i-1, j-1
        elif score_current == score_up + penalty["gap"]:
            print("left")
            align1, align2 = "-", seq2[j-1]
            j -= 1
        total_align1 += align1
        total_align2 += align2

    while i > 0:
        align1, align2 = seq1[i-1], "-"
        total_align1 += align1
        total_align2 += align2
        i -= 1

    while j > 0:
        align1, align2 = "-", seq2[j-1]
        total_align1 += align1
        total_align2 += align2
        j -=1

    total_align1 = total_align1[::-1]
    total_align2 = total_align2[::-1]
    seqN = len(align1)
    sym = ""
    seq_score = 0
    ident = 0
    for i in range(seqN):
        align1 = total_align1[i]
        align2 = total_align2[i]
        if align1 == align2:
            sym += align1
            ident += 1
            seq_score += nwAlign(align1, align2)

        else:
            seq_score += nwAlign(align1, align2)
            sym += ""
    ident = ident/seqN * 100
