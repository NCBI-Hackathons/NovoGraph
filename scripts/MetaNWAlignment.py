#!/usr/bin/python3

###########################
#
# MNW
# Meta Needleman-Wunsch Algorithm 
#
# Needleman-Wunsch applied to match contigs to the reference via coordinate
# In order to find the global alignment, we use the following scoring system:
# Correct description to be added 
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

# simulated data
start = np.random.randint(10000, size=(10,2))
# array([ 413, 1800,  331, 8910,  837, 9660, 2167, 5997, 3104, 7259])
end   = start + np.random.randint(1000, size=(10,2))
score = np.random.randint(1000, size=(10,1))
startContig = 1;
lengthContig = 10000
  
trash = start[:, 0]                             # takes the first "column" of `start`   
# Arrange local allignments by the start posisiton 
trash = trash.ravel().argsort()  
start = start[trash, :]
end   = end[trash, :]
score = score[trash, :]
    


def globalPath(startContig, lenghtContig, start, end,score):
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
    # based on definition above `size=(10,2)`, len(start) = 10 always
    ndiag = len(start)     
    # list partition --- create an array of arrays #### in MATLAB, num2cell()
    # possiblePath = array([[0],[1],[2],[3],[4],[5],[6],[7],[8],[9])
    possiblePath = list(np.array([[x] for x in range(0, ndiag)])) 
    # array([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1]])                        
    pathExtension = np.ones((ndiag,), dtype=int)  
    while True:
        extensionPathIndx = np.nonzero(pathExtension==1)   # (array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),)
        nPossiblePath = len(possiblePath)  # len() = 10
        if len(extensionPathIndx)==0:
            # array([0]) --- THIS IS THE START ELEMENT OF extensionPathIndx]
            presentPath = possiblePath[extensionPathIndx[0][0]] 
            presentDiagEnd = end[presentPath[-1],:]             # takes first pair from end
            potentialDiag = np.logical_and(start[:,1] > presentDiagEnd[1] ,  start[:,0] > presentDiagEnd[0])
            # e.g. array([False,  True, False, False, False, False, False, False, False, False], dtype=bool)
            potentialDiagIndx = [i for i, x in enumerate(potentialDiag) if x] # [1]
            lenPotential = np.sum(potentialDiag) # 1
            if lenPotential == 0:
                pathExtension[extensionPathIndx[0][0]] =  0
            else:
                pathToAppend = possiblePath[extensionPathIndx[0][0]]  # array([0])
                newPathToAdd = [];
                for i in range(0, lenPotential):  # lenPotential=1; for i in #output== 'range(0, 0)':
                    newPathToAdd.append(np.append(pathToAppend, potentialDiagIndx[i]) )   #pathToAppend
                del possiblePath[extensionPathIndx[0][0]]
                possiblePath = possiblePath[0:extensionPathIndx[0] [0]]+newPathToAdd+possiblePath[extensionPathIndx[0][0]:]
                pathExtension = np.append(pathExtension,np.ones((lenPotential-1,), dtype=int))
        else:
            break      


    nPossiblePath = len(possiblePath)  # always 
    endContig = startContig+lenghtContig-1
    newScore = [];
    for i in range(0, nPossiblePath-1): # 0 1 2 3 4 5 6 7
        path = possiblePath[i]          #possiblePath[0]
        pathStarts = start[path,:]      # grab position in `start`
        pathEnds = end[path,:]
                            # newScore = sum
        newScore.append(np.sum(score[path])-(pathStarts[0,0]-startContig) - np.sum(pathStarts[1:,0]-pathEnds[0:-1,0])- np.sum(pathStarts[1:,1]-pathEnds[0:-1,1])+ endContig-pathEnds[-1,0])
    
     
    
if __name__ == "__main__":
    globalPath(startContig, lengthContig, start, end,score)
 


