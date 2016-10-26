#!/usr/bin/python3

#
# JGA
# Jajaoo Global Alignment algorithm
#

import numpy as np
import random

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
      PotentialDiagIndx = [i for i, x in enumerate(PotentialDiag) if x]
      lenPotential = PotentialDiag.sum()
      if lenPotential==0:
         PathExtension[0][ExtensionPathIndx[0][0]] =  0
      else
         PossiblePath[range(ExtensionPathIndx[0][0]+lenPotential,nPossiblePath+lenPotential-1)] = PossiblePath[ExtensionPathIndx[0][0]+1:nPossiblePath]
         PathToAppend = PossiblePath[ExtensionPathIndx[0][0]]
         for i in range(0,lenPotential-1):
             PossiblePath[ExtensionPathIndx[0][0]+i] = np.append(PathToAppend,PotentialDiagIndx[i])
    else
        Stop = 1;
    clear trash

  nPossiblePath = len(PossiblePath)
  StartPositionContig = 1
  EndPositionContig = 10000
  for i in range(0,nPossiblePath-1)
    Path = PossiblePath[i];
    PathStarts = Start[Path,:]
    PathEnds = End[Path,:]
    Score[i] = (PathStarts[0,0]-StartPositionContig)+sum(PathStarts[1:,0]-PathEnds[0:-1,0])+sum(PathStarts[1:,1]-PathEnds[0:-1,1])+ EndPositionContig-PathEnds[-1,0];
    Score[i] = -Score[i]
            

