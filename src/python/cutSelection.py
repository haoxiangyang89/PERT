#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 16:15:31 2019

@author: haoxiangyang
"""
from part_tight import revPar

# cut selection codes
def examineCuts(pData,disData,Omega,cutSet,divSet,that,xhat,thetahat,yhat):
    # check for each cut whether it is tight, if not update the counts
    cutSel = []
    for nc in range(len(cutSet)):
        # how many rounds have been through
        for l in range(len(cutSet[nc][1])):
            omega = cutSet[nc][1][l][0]
            cutV = cutSet[nc][1][l][1]
            # add t
            for i in pData.II:
                cutV += cutSet[nc][1][l][2][i]*(that[i] - cutSet[nc][0][0][i])
                # add x
                for j in pData.Ji[i]:
                    cutV += cutSet[nc][1][l][3][i,j]*(xhat[i,j] - cutSet[nc][0][1][i,j])
                for par in range(len(cutSet[nc][0][3][i])):
                    cutV += cutSet[nc][1][l][4][i,par]*(sum(yhat[i,parNew] for parNew in range(len(divSet[i])) \
                                  if revPar(cutSet[nc][0][3][i],divSet[i][parNew]) == par) - cutSet[nc][0][2][i,par])
            if abs(thetahat[omega] - cutV)/thetahat[omega] <= 1e-4:
                # tight
                cutSel.append((nc,l))
    return cutSel
    
def selectCuts(cutSet,cutSetSel):
    cutSetNew = []
    for nc in range(len(cutSet)):
        cutSetNew.append([cutSet[nc][0],[cutSet[nc][1][l] for l in range(len(cutSet[nc][1])) if (nc,l) in cutSetSel]])
    return cutSetNew
        