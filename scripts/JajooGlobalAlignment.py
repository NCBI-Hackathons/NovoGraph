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
# Aarti Jajoo (Baylor)
# Evan Biederstedt (NYGC)
#
###########################

import numpy as np


def globalPath(startContig,lenghtContig,start,end):
    # parameter explanation
    # start
    # end
    # stop
    # stop
    start = np.random.randint(10000, size=(10,2))
    end   = np.random.randint(10000, size=(10,2))
    score = np.random.randint(10000, size=(10,1))
    stop = 0
    
    # reindex
    trash = start[:, 0]
    trash = trash.ravel().argsort()
    start = start[trash, :]
    end   = end[trash, :]
    score = score[trash, :]

    ndiag = len(start)     # based on definition above `size=(10,2)`, len(start) = 10 always
    possiblePath = np.array([[x] for x in range(0, ndiag-1)])
        # possiblePath = array([[0],[1],[2],[3],[4],[5],[6],[7],[8]])
    pathExtension = np.ones((1, ndiag), dtype=int)  # array([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1]])
    while stop == 0:
        extensionPathIndx = np.nonzero(pathExtension[0]==1)   # (array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),)
        nPossiblePath = len(possiblePath)  # len() = 9
        if extensionPathIndx:
            presentPath = possiblePath[extensionPathIndx[0][0]] # array([0])
            presentDiagEnd = end[presentPath[-1],:]             # takes first pair from end
            potentialDiag = np.logical_and(start[:,1] > presentDiagEnd[1] ,  start[:,0] > presentDiagEnd[0])
            potentialDiagIndx = [i for i, x in enumeratepPotentialDiag) if x]
            lenPotential = potentialDiag.sum()
        if lenPotential == 0:
            pathExtension[0][extensionPathIndx[0][0]] =  0
        else
            possiblePath[range(extensionPathIndx[0][0]+lenPotential, nPossiblePath+lenPotential - 1)] = PossiblePath[extensionPathIndx[0][0] + 1:nPossiblePath]
            pathToAppend = PossiblePath[ExtensionPathIndx[0][0]]
            for i in range(0,lenPotential-1):
                PossiblePath[ExtensionPathIndx[0][0]+i] = np.append(PathToAppend,PotentialDiagIndx[i])
        else
            Stop = 1;
    clear trash

    nPossiblePath = len(PossiblePath)
    StartPositionContig = 1
    EndPositionContig = 10000
    for i in range(0,nPossiblePath-1):
        Path = PossiblePath[i];
        PathStarts = Start[Path,:]
        PathEnds = End[Path,:]
        Score[i] = (PathStarts[0,0]-StartPositionContig)+sum(PathStarts[1:,0]-PathEnds[0:-1,0])+sum(PathStarts[1:,1]-PathEnds[0:-1,1])+ EndPositionContig-PathEnds[-1,0];
        Score[i] = -Score[i]


"""
    def main():
    
    
    if __name__ == "__main__":
    main()
 
"""

