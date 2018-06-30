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
    # start --- start coordinate of contig X
    # end --- length of contig X
    # start ---
    # end ---
    
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
    start = np.random.randint(10000, size=(10,2))
    trash = start[:, 0]                             # takes the first "column" of `start`
    # array([ 413, 1800,  331, 8910,  837, 9660, 2167, 5997, 3104, 7259])
    end   = np.random.randint(10000, size=(10,2))
    score = np.random.randint(10000, size=(10,1))
  
    
    # reindex
    trash = trash.ravel().argsort()  # trash MUST be a numpy array
    start = start[trash, :]
    end   = end[trash, :]
    score = score[trash, :]
    

    ndiag = len(start)     # based on definition above `size=(10,2)`, len(start) = 10 always
    possiblePath = list(np.array([[x] for x in range(0, ndiag)])) # possiblePath = array([[0],[1],[2],[3],[4],[5],[6],[7],[8],[9])
                            # list partition --- create an array of arrays #### in MATLAB, num2cell()
    pathExtension = np.ones((1, ndiag), dtype=int)  # array([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1]])
    while stop == False:
        extensionPathIndx = np.nonzero(pathExtension[0]==1)   # (array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),)
        nPossiblePath = len(possiblePath)  # len() = 10
        if extensionPathIndx:
            presentPath = possiblePath[extensionPathIndx[0][0]] # array([0]) --- THIS IS THE START ELEMENT OF extensionPathIndx
            presentDiagEnd = end[presentPath[-1],:]             # takes first pair from end, s.t. end = np.random.randint(10000, size=(10,2)); array([7133, 5983])
            potentialDiag = np.logical_and(start[:,1] > presentDiagEnd[1] ,  start[:,0] > presentDiagEnd[0])
            # e.g. array([False,  True, False, False, False, False, False, False, False, False], dtype=bool)
            potentialDiagIndx = [i for i, x in enumerate(potentialDiag) if x] # [1]
            lenPotential = np.sum(potentialDiag) # 1
        if lenPotential == 0:
            pathExtension[0][extensionPathIndx[0][0]] =  0
        else
            possiblePath.append(possiblePath[extensionPathIndx[0][0] + 1:nPossiblePath])
                            #array([[0],[1],[2],[3],[4],[5],[6],[7],[8]])
            pathToAppend = possiblePath[extensionPathIndx[0][0]]  # array([0])
            for i in range(0, lenPotential-1):  # lenPotential=1; for i in #output== 'range(0, 0)':
                possiblePath[extensionPathIndx[0][0] + i] = np.append(pathToAppend, potentialDiagIndx[i])  #pathToAppend

             """
                 PossiblePath(ExtensionPathIndx(1)+lenPotential:nPossiblePath+lenPotential-1) = PossiblePath(ExtensionPathIndx(1)+1:nPossiblePath) ; # EXTEND! # [ [1], [2], [2], [3]]
                 trash = repmat(PossiblePath{ExtensionPathIndx(1)},lenPotential,1);  # repmat() --> [[1], [1]]
                 trash = [trash find(PotentialDiag==1)];  # find [[1 3], [1 4]]
                 [row,col] = size(trash);  # replaces here, outputs [[1 3], [1 4], [2], [3]]
                 PossiblePath(ExtensionPathIndx(1):ExtensionPathIndx(1)+lenPotential-1) = mat2cell(trash,repmat(1,1,lenPotential),col);
             """


        else
            break # stop = True


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

