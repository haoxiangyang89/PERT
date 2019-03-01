#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 16:47:18 2019

@author: haoxiangyang
"""

from tighten import *
from part_tight import *
from ubCalFunc import *
import copy
import numpy as np
from gurobipy import *
import time
from cutSelection import *
from sub import *
import pdb
pdb.set_trace()

def partBenders(model, where):
    if where == GRB.Callback.MIPSOL:
        # the callback function
        that = {}
        xhat = {}
        thetahat = {}
        yhat = {}
        # obtain the solution at the current node
        for i in pData.II:
            that[i] = model.cbGetSolution(tVar[i])
            for j in pData.Ji[i]:
                xhat[i,j] = model.cbGetSolution(xVar[i,j])
            for par in range(len(divSet[i])):
                yhat[i,par] = round(model.cbGetSolution(yVar[i,par]))
        for omega in Omega:
            thetahat[omega] = model.cbGetSolution(thetaVar[omega])
        tcoreList.append(that)
        xcoreList.append(xhat)
        ycoreList.append(yhat)
    
        # generate cuts
        pidict = {}
        lambdadict = {}
        gammadict = {}
        vk = {}
        thetaInt = {}
        ubCost = min(ubCostList)
        ubTemp,thetaInt = ubCal(pData,disData,Omega,xhat,that,Tmax1,1)
        if ubCost > ubTemp:
            for i in pData.II:
                tbest[i] = that[i]
                for j in pData.Ji[i]:
                    xbest[i,j] = xhat[i,j]
            for omega in Omega:
                thetabest[omega] = thetaInt[omega]
        ubCostList.append(ubTemp)
    
        # obtain the cores
        tcore,xcore,ycore = avgCore(pData,divSet,tcoreList,xcoreList,ycoreList)
        dataList = {}
        GCurrent = {}
        for omega in Omega:
            dataItem = sub_divDual(pData,disData[omega],that,xhat,yhat,divSet,H,lDict,tcore,xcore,ycore)
            pidict[omega] = dataItem[0]
            lambdadict[omega] = dataItem[1]
            gammadict[omega] = dataItem[2]
            vk[omega] = dataItem[3]
            GCurrent[omega] = dataItem[4]
        tV = {}
        xV = {}
        yV = {}
        thetaV = {}
        for i in pData.II:
            tV[i] = model.getVarByName("t[{}]".format(i))
            for j in pData.Ji[i]:
                xV[i,j] = model.getVarByName("x[{},{}]".format(i,j))
            for par in range(len(divSet[i])):
                yV[i,par] = model.getVarByName("y[{},{}]".format(i,par))
        for omega in Omega:
            thetaV[omega] = model.getVarByName("theta[{}]".format(omega))
        cutDual = []
        for omega in Omega:
            if vk[omega] - thetahat[omega] > 1e-4*thetahat[omega]:
                cutDual.append([omega,vk[omega],pidict[omega],lambdadict[omega],gammadict[omega]])
                model.cbLazy(thetaV[omega] >= vk[omega] + sum(pidict[omega][i]*(tV[i] - that[i]) for i in pData.II) +\
                    sum(sum(lambdadict[omega][i,j]*(xV[i,j] - xhat[i,j]) for j in pData.Ji[i]) for i in pData.II) +\
                    sum(sum(gammadict[omega][i,par]*(yV[i,par] - yhat[i,par]) for par in range(len(divSet[i]))) for i in pData.II))
                
        cutSet.append([[that,xhat,yhat,divSet],cutDual])
        GList.append(GCurrent)
        print("lazy",model.cbGet(GRB.Callback.MIPSOL_OBJBND))

# partion solution process
Tmax = disData[len(Omega) - 1].H + longestPath(pData)[0]
pdData = copy.deepcopy(pData)
for i in pData.II:
    if i != 0:
        pdData.D[i] = pData.D[i] + max([disData[omega].d[i] for omega in Omega])
    else:
        pdData.D[i] = pData.D[i]
        
lDict = longestPath(pdData)
Tmax1 = lDict[0]

allSucc = findSuccAll(pData)
distanceDict = {}
for i in pData.II:
    for j in allSucc[i]:
        distanceDict[i,j] = detCal(pData,i,j)

H = {}
H[0] = 0
counter = 1
for omega in Omega:
    if not(disData[omega].H in H.values()):
        H[counter] = disData[omega].H
        counter += 1
H[counter] = Tmax
HOmega = range(1,counter)
HRev = {}
for hIter in H.keys():
    HRev[H[hIter]] = hIter

# start with an upper bound based on the deterministic solution
tdet,xdet,fdet = detBuild(pData)
ubdet = ubCal(pData,disData,Omega,xdet,tdet,Tmax1)
brInfo = precludeRel(pData,H,ubdet)

# initialize the cutSet and divSet
cutSet = []
divSet = {}
divDet = {}
for i in pData.II:
    set1 = [h for h in HOmega if brInfo[pData.II.index(i),h - 1] == 1]
    setn1 = [h for h in HOmega if brInfo[pData.II.index(i),h - 1] == -1]
    
    if set1 != []:
        set1t = partType(0,max(set1))
        if setn1 != []:
            setn1t = partType(min(setn1),len(HOmega) + 1)
            set0t = partType(max(set1),min(setn1))
            divSet[i] = [set1t,set0t,setn1t]
            divDet[i] = [1,0,-1]
        else:
            set0t = partType(max(set1),len(HOmega) + 1)
            divSet[i] = [set1t,set0t]
            divDet[i] = [1,0]
    else:
        if setn1 != []:
            setn1t = partType(min(setn1),len(HOmega) + 1)
            set0t = partType(0,min(setn1))
            divSet[i] = [set0t,setn1t]
            divDet[i] = [0,-1]
        else:
            set0t = partType(0,len(HOmega) + 1)
            divSet[i] = [set0t]
            divDet[i] = [0]
            
# initiate the process
lbCost = -np.inf
lbCostList = [lbCost]
GList = []
cutSel = {}
cutThreshold = 10
tcoreList = []
xcoreList = []
ycoreList = []

# keep the record of the progress
keepIter = True
lbHist = []
ubHist = []
timeHist = []
cutHist = []
intSolHist = []
yhistList = []
noThreads = 4
epsilon = 1e-2
sN = 10
MM = 5

# ubList is the list of upper bounds, tHList is the list of each activity's range of starting time
ubList,tHList,ubbest,tbest,xbest,thetabest,textList,xextList = iniPart(pData,disData,Omega,sN,MM,1)
# pre-separate the partition
divSet,divDet = splitPar_CI(divSet,divDet,tHList)
for ii in range(len(textList)):
    tcoreList.append(textList[ii])
    xcoreList.append(xextList[ii])
    ycoreItem = {}
    for i in pData.II:
        for par in range(len(divSet[i])):
            if (textList[ii][i] >= H[divSet[i][par].startH])and(textList[ii][i] < H[divSet[i][par].endH]):
                ycoreItem[i,par] = 1
            else:
                ycoreItem[i,par] = 0
    ycoreList.append(ycoreItem)

ubCost = min(ubbest,ubdet)
ubCostList = [ubCost]

mp = Model("partitionSolve")
mp.Params.IntFeasTol = 1e-8
mp.Params.FeasibilityTol = 1e-8
mp.Params.Threads = noThreads
mp.Params.LazyConstraints = 1

# add master variables
thetaVar = mp.addVars(Omega, lb = 0, name = "theta")
JJ = [(i,j) for i in pData.II for j in pData.Ji[i]]
tVar = mp.addVars(pData.II, lb = 0, name = "t")
xVar = mp.addVars(JJ, lb = 0, ub = 1, name = "x")
parI = [(i,par) for i in pData.II for par in range(len(divSet[i]))]
yVar = mp.addVars(parI, vtype = GRB.BINARY, name = "y")

# set the starting point of the solution from iniPart
#ybest = {}
#for i in pData.II:
#    tVar[i].Start = tbest[i]
#    for j in pData.Ji[i]:
#        xVar[i,j].Start = xbest[i,j]
#    for par in range(len(divSet[i])):
#        if (tbest[i] >= H[divSet[i][par].startH])and(tbest[i] < H[divSet[i][par].endH]):
#            yVar[i,par].Start = 1
#            ybest[i,par] = 1
#        else:
#            yVar[i,par].Start = 0
#            ybest[i,par] = 0
#for omega in Omega:
#    #thetaVar[omega].Start = thetabest[omega]
#    thetaVar[omega].Start = sub_div(pData,disData[omega],tbest,xbest,ybest,divSet,H,lDict)

# add constraints
budgetConstr = mp.addConstr((sum(sum(pData.b[i][j]*xVar[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B), name = "budgetConstr")
durationConstr = mp.addConstrs((tVar[k[1]] - tVar[k[0]] >= pData.D[k[0]]*(1 - sum(pData.eff[k[0]][j]*xVar[k[0],j] for j in pData.Ji[k[0]])) \
                                for k in pData.K), name = "durationConnstr")
xConstr = mp.addConstrs((sum(xVar[i,j] for j in pData.Ji[i]) <= 1 for i in pData.II), name = "xConstr")
tubConstr = mp.addConstrs((tVar[i] <= sum(H[divSet[i][par].endH]*yVar[i,par] for par in range(len(divSet[i]))) for i in pData.II), name = "tub")
tlbConstr = mp.addConstrs((tVar[i] >= sum(H[divSet[i][par].startH]*yVar[i,par] for par in range(len(divSet[i]))) for i in pData.II), name = "tlb")
yConstr = mp.addConstrs((sum(yVar[i,par] for par in range(len(divSet[i]))) == 1 for i in pData.II), name = "yConstr")
yLimit = mp.addConstrs((yVar[i,par] == 0 for i in pData.II for par in range(len(divSet[i])) if divDet[i][par] != 0), name = "yLimit")
mp.setObjective(pData.p0*tVar[0] + sum(disData[omega].prDis*thetaVar[omega] for omega in Omega), GRB.MINIMIZE)

while keepIter:
    tCurrent = {}
    xCurrent = {}
    thetaCurrent = {}
    yCurrent = {}
    
    # add the cut
    # cutInfo = 2 dimensional vector, first dimention record the primal solution,
    # second dimension record the dual solution for every scenario
    mp.addConstrs((thetaVar[cutSet[nc][1][l][0]] >= cutSet[nc][1][l][1] + sum(cutSet[nc][1][l][2][i]*(tVar[i] - cutSet[nc][0][0][i]) +\
                   sum(cutSet[nc][1][l][3][i,j]*(xVar[i,j] - cutSet[nc][0][1][i,j]) for j in pData.Ji[i]) +\
                   sum(cutSet[nc][1][l][4][i,par]*(sum(yVar[i,parNew] for parNew in range(len(divSet[i])) if revPar(cutSet[nc][0][3][i],divSet[i][parNew]) == par) - cutSet[nc][0][2][i,par])\
                       for par in range(len(cutSet[nc][0][3][i]))) for i in pData.II)\
                   for nc in range(len(cutSet)) for l in range(len(cutSet[nc][1]))), name = "cuts")
    
    # add the constraints between y
    # obtain the set of y's all predecessors
    yLogic = {}
    for i in pData.II:
        # for each precedence relationship
        for j in allSucc[i]:
            k = (i,j)
            for par1 in range(len(divSet[k[0]])):
                for par2 in range(len(divSet[k[1]])):
                    if H[divSet[k[1]][par2].endH] < H[divSet[k[0]][par1].startH] + distanceDict[i,j]:
                        yLogic[k,par1,par2] = mp.addConstr(yVar[k[0],par1] + yVar[k[1],par2] <= 1)
    
    # set the starting point to the best among all core points
    lbmin = np.inf
    coreMinind = -1
    thetabest = {}
    for ii in range(len(tcoreList)):
        thetacore = {}
        for omega in Omega:
            thetacore[omega] = sub_div(pData,disData[omega],tcoreList[ii],xcoreList[ii],ycoreList[ii],divSet,H,lDict)
        lbcore = pData.p0*(tcoreList[ii][0]) + sum(thetacore[omega]*disData[omega].prDis for omega in Omega)
        if lbcore < lbmin:
            lbmin = lbcore
            coreMinind = ii
            thetabest = thetacore
    for i in pData.II:
        tVar[i].Start = tcoreList[coreMinind][i]
        for j in pData.Ji[i]:
            xVar[i,j].Start = xcoreList[coreMinind][i,j]
        for par in range(len(divSet[i])):
            yVar[i,par].Start = ycoreList[coreMinind][i,par]
    for omega in Omega:
        thetaVar[omega].Start = thetabest[omega]
    
    mp.update()
    tic = time.time()
    mp.optimize(partBenders)
    tIter = time.time() - tic
    timeHist.append(tIter)
    lbCost = mp.objVal
    lbHist.append(lbCost)
    for i in pData.II:
        tCurrent[i] = tVar[i].X
        for j in pData.Ji[i]:
            xCurrent[i,j] = xVar[i,j].X
        for par in range(len((divSet[i]))):
            yCurrent[i,par] = yVar[i,par].X
    for omega in Omega:
        thetaCurrent[omega] = thetaVar[omega].X
    
    # refine the partition
    ubCost = min(ubCostList)
    ubHist.append(ubCost)
    intSolHist.append(len(ubCostList))
    if (ubCost - lbCost)/ubCost < epsilon:
        keepIter = False
    else:
        cutSel = examineCuts(pData,disData,Omega,cutSet,divSet,tCurrent,xCurrent,thetaCurrent,yCurrent)
        cutSetNew = selectCuts(cutSet,cutSel)
        GCurrent = GList[-1]
        GFrac = {}
        for i in pData.II:
            GFraciList = [disData[omega].H for omega in Omega if (GCurrent[omega][i] < 1 - 1e-6)and(GCurrent[omega][i] > 1e-6)]
            if GFraciList != []:
                GFrac[i] = [HRev[min(GFraciList)],HRev[max(GFraciList)]]
            else:
                GFrac[i] = []
        newPartition = []
        for i in pData.II:
            if GFrac[i] != []:
                newItem = (i, GFrac[i][0], GFrac[i][1])
                newPartition.append(newItem)
        divSet,divDet = splitPar3(divSet,divDet,newPartition)
        cutHist.append(sum(len(cutSet[l][1]) for l in range(len(cutSet))))
        cutSet = copy.deepcopy(cutSetNew)
    
    # correct all ycoreList
    for ll in range(len(ycoreList)):
        for i in pData.II:
            for par in range(len(divSet[i])):
                if (tcoreList[ll][i] >= H[divSet[i][par].startH])and(tcoreList[ll][i] < H[divSet[i][par].endH]):
                    ycoreList[ll][i,par] = 1
                else:
                    if (divSet[i][par].endH == len(Omega) + 1)and(abs(tcoreList[ll][i] - H[len(Omega) + 1]) <= 1e-4):
                        ycoreList[ll][i,par] = 1
                    else:
                        ycoreList[ll][i,par] = 0
    
    # redefine the master with the updated partition
    mp = Model("partitionSolve")
    mp.Params.IntFeasTol = 1e-8
    mp.Params.FeasibilityTol = 1e-8
    mp.Params.Threads = noThreads
    mp.Params.LazyConstraints = 1
    
    # add master variables
    thetaVar = mp.addVars(Omega, lb = 0, name = "theta")
    JJ = [(i,j) for i in pData.II for j in pData.Ji[i]]
    tVar = mp.addVars(pData.II, lb = 0, name = "t")
    xVar = mp.addVars(JJ, lb = 0, ub = 1, name = "x")
    parI = [(i,par) for i in pData.II for par in range(len(divSet[i]))]
    yVar = mp.addVars(parI, vtype = GRB.BINARY, name = "y")
    
    # add constraints
    budgetConstr = mp.addConstr((sum(sum(pData.b[i][j]*xVar[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B), name = "budgetConstr")
    durationConstr = mp.addConstrs((tVar[k[1]] - tVar[k[0]] >= pData.D[k[0]]*(1 - sum(pData.eff[k[0]][j]*xVar[k[0],j] for j in pData.Ji[k[0]])) \
                                    for k in pData.K), name = "durationConnstr")
    xConstr = mp.addConstrs((sum(xVar[i,j] for j in pData.Ji[i]) <= 1 for i in pData.II), name = "xConstr")
    tubConstr = mp.addConstrs((tVar[i] <= sum(H[divSet[i][par].endH]*yVar[i,par] for par in range(len(divSet[i]))) for i in pData.II), name = "tub")
    tlbConstr = mp.addConstrs((tVar[i] >= sum(H[divSet[i][par].startH]*yVar[i,par] for par in range(len(divSet[i]))) for i in pData.II), name = "tlb")
    yConstr = mp.addConstrs((sum(yVar[i,par] for par in range(len(divSet[i]))) == 1 for i in pData.II), name = "yConstr")
    yLimit = mp.addConstrs((yVar[i,par] == 0 for i in pData.II for par in range(len(divSet[i])) if divDet[i][par] != 0), name = "yLimit")
    mp.setObjective(pData.p0*tVar[0] + sum(disData[omega].prDis*thetaVar[omega] for omega in Omega), GRB.MINIMIZE)