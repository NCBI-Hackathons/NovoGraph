# implementation of MetaNWAlignment.py, globalPath(startContig, lenghtContig, start, end, score)

import numpy as np

score = np.random.randint(1000, size=(10,1))
end   = start + np.random.randint(1000, size=(10,2))
lengthContig = 10000
startContig = 1
start = np.random.randint(10000, size=(10,2))

def globalPath(startContig, lenghtContig, start, end, score):
    """
    parameters:
    # simulated data
    start = np.random.randint(10000, size=(10,2))
    # array([ 413, 1800,  331, 8910,  837, 9660, 2167, 5997, 3104, 7259])
    end   = start + np.random.randint(1000, size=(10,2))
    score = np.random.randint(1000, size=(10,1))
    startContig = 1;
    lengthContig = 10000
    
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

    return newScore


if __name__ == "__main__":
    globalPath(startContig, lengthContig, start, end, score)
