#!/usr/bin/python3

###########################
#
# JGA
# Jajaoo Global Alignment algorithm
#
# Needleman-Wunsch applied to match contigs to the reference via coordinate
# In order to find the global alignment, we use the following scoring system:
# Map contigs by coordinates against reference by coordinates:
#    (start_seq1, start_ref1) to (end_seq2, start_ref2)
# For each subgroup (i.e. unique chromsome # and +/-):
# ---if diagonal, add +len(diagonal) # units bases
# ---if vertical/horizontl, -len(distance) # units bases
#
# Aarti Jajoo (Baylor College of Medicine)
# Evan Biederstedt (New York Genome Center, Weill Cornell)
#
###########################

import numpy as np


def globalPath(startContig, lenghtContig, start, end):
    # parameter explanation
    # start
    # end
    # stop
    # stop
    start = np.random.randint(10000, size=(10,2))
    """
        print(start)
        
        
        array([[ 413, 5574],
        [1800, 8123],
        [ 331, 8414],
        [8910, 9772],
        [ 837, 1202],
        [9660, 6792],
        [2167,  278],
        [5997, 4494],
        [3104, 3755],
        [7259, 6283]])
        
    """
    trash = start[:, 0]                             # takes the first "column" of `start`
    # array([ 413, 1800,  331, 8910,  837, 9660, 2167, 5997, 3104, 7259])
    end   = np.random.randint(10000, size=(10,2))
    score = np.random.randint(10000, size=(10,1))
    stop = 0
    
    # reindex
    trash = trash.ravel().argsort()  # trash MUST be a numpy array
    start = start[trash, :]
    end   = end[trash, :]
    score = score[trash, :]

    ndiag = len(start)     # based on definition above `size=(10,2)`, len(start) = 10 always
    possiblePath = np.array([[x] for x in range(0, ndiag-1)] # possiblePath = array([[0],[1],[2],[3],[4],[5],[6],[7],[8]])
                            # list partition --- create an array of arrays #### in MATLAB, num2cell()
    pathExtension = np.ones((1, ndiag), dtype=int)  # array([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1]])
    while stop == True:
        extensionPathIndx = np.nonzero(pathExtension[0]==1)   # (array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),)
        nPossiblePath = len(possiblePath)  # len() = 9
        if extensionPathIndx:
            presentPath = possiblePath[extensionPathIndx[0][0]] # array([0])
            presentDiagEnd = end[presentPath[-1],:]             # takes first pair from end, s.t. end = np.random.randint(10000, size=(10,2))
            potentialDiag = np.logical_and(start[:,1] > presentDiagEnd[1] ,  start[:,0] > presentDiagEnd[0]) # e.g. array([False,  True, False, False, False, False, False, False, False, False], dtype=bool)
            potentialDiagIndx = [i for i, x in enumerate(potentialDiag) if x] # [1]
            lenPotential = potentialDiag.sum() # 1
        if lenPotential == 0:
            pathExtension[0][extensionPathIndx[0][0]] =  0
        else
            possiblePath[range(extensionPathIndx[0][0]+lenPotential, nPossiblePath+lenPotential - 1)] = possiblePath[extensionPathIndx[0][0] + 1:nPossiblePath]
                            #array([[0],[1],[2],[3],[4],[5],[6],[7],[8]])
            pathToAppend = possiblePath[extensionPathIndx[0][0]]  # array([0])
            for i in range(0, lenPotential-1):  # lenPotential=1; for i in #output== 'range(0, 0)':
                possiblePath[extensionPathIndx[0][0] + i] = np.append(pathToAppend, potentialDiagIndx[i])  #pathToAppend
        else
            stop = False


    nPossiblePath = len(possiblePath)  # always 9
    startPositionContig = 1
    endPositionContig = 10000
    for i in range(0, nPossiblePath-1): # 0 1 2 3 4 5 6 7
        path = possiblePath[i]          #possiblePath[0]
        pathStarts = start[path,:]      # grab position in `start`
        pathEnds = end[path,:]
                            # newScore = sum
        newScore[i] = (np.sum(score[path])-(pathStarts[0,0]-startPositionContig) - np.sum(pathStarts[1:,0]-pathEnds[0:-1,0])
                           - np.sum(pathStarts[1:,1]-pathEnds[0:-1,1])+ endPositionContig-pathEnds[-1,0]))


"""
    def main():
    
    
    if __name__ == "__main__":
    main()
 
"""

