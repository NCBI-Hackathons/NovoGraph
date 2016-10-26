# 
# Needlman Wunsch Algo
import numpy as np
import random
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# 
def GlobalPath(StartContig,LenghtContig,Start,End):
  Start = numpy.random.randint(10000, size=(10,2))
  End   = numpy.random.randint(10000, size=(10,2))
  Stop = 0
  ndiag = len(Start)
  PossiblePath = numpy.array([[x] for x in range(0,ndiag-1)])
  PathExtension = numpy.ones((1,ndiag))
  while Stop == 0:
    ExtensionPathIndx = numpy.nonzero(PathExtension[0]==1)
    nPossiblePath = len(PossiblePath);
    if ExtensionPathIndx:
      PresentPath = PossiblePath[ExtensionPathIndx[0][0]]
      PresentDiagEnd = End[PresentPath[-1],:]
      PotentialDiag = numpy.logical_and(Start[:,1]>PresentDiagEnd[1] ,  Start[:,0] > PresentDiagEnd[0])
      lenPotential = PotentialDiag.sum()
      if lenPotential==0:



        if lenPotential ==0
            PathExtension(ExtensionPathIndx(1)) = 0;
    else
        PossiblePath(ExtensionPathIndx(1)+lenPotential:nPossiblePath+lenPotential-1) = PossiblePath(ExtensionPathIndx(1)+1:nPossiblePath) ;
            trash = repmat(PossiblePath{ExtensionPathIndx(1)},lenPotential,1);
            trash = [trash find(PotentialDiag==1)];
            [row,col] = size(trash);
            PossiblePath(ExtensionPathIndx(1):ExtensionPathIndx(1)+lenPotential-1) = mat2cell(trash,repmat(1,1,lenPotential),col);
        end









