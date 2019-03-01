#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 17:52:06 2019

@author: haoxiangyang
"""
from gurobipy import *

# deterministic formulations

def detBuild(pData):
    # construct the deterministic crashing problem
    mp = Model('deterministic_model')
    mp.setParam('OutputFlag',0)
    
    # set up variables
    t = mp.addVars(pData.II, lb = 0, name = "t")
    JJ = [(i,j) for i in pData.II for j in pData.Ji[i]]
    x = mp.addVars(JJ, lb = 0, ub = 1, name = "x")
    
    # set up constraints
    durationConstr = mp.addConstrs((t[k[1]] - t[k[0]] >= pData.D[k[0]]*(1 - sum(pData.eff[k[0]][j]*x[k[0],j] for j in pData.Ji[k[0]])) for k in pData.K), name = "durationConstr")
    budgetConstr = mp.addConstr((sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B), name = "budgetConstr")
    xConstr = mp.addConstrs((sum(x[i,j] for j in pData.Ji[i]) <= 1 for i in pData.II), name = "xConstr")
    mp.setObjective(t[0],GRB.MINIMIZE)
    mp.update()
    mp.optimize()
    
    tdet = {}
    xdet = {}
    fdet = mp.objVal
    for i in pData.II:
        tdet[i] = t[i].X
    for j in JJ:
        xdet[j] = x[j].X
    return tdet,xdet,fdet

def expModel(pData,eH,ed,M = 9999999):
    # script for two-stage expected time expected magnitude disruption program
    mp = Model('expected_model')
    mp.Params.OutputFlag = 0
    
    # set up variables
    t0 = mp.addVars(pData.II, lb = 0, name = "t0")
    JJ = [(i,j) for i in pData.II for j in pData.Ji[i]]
    x0 = mp.addVars(JJ, lb = 0, ub = 1, name = "x0")
    t = mp.addVars(pData.II, lb = 0, name = "t")
    x = mp.addVars(JJ, lb = 0, ub = 1, name = "x")
    G = mp.addVars(pData.II, vtype=GRB.BINARY, name = "G")
    s = mp.addVars(JJ, lb = 0, ub = 1, name = "s")

    # set up constraints
    durationConstr0 = mp.addConstrs((t0[k[1]] - t0[k[0]] >= pData.D[k[0]]*(1 - sum(pData.eff[k[0]][j]*x0[k[0],j] for j in pData.Ji[k[0]])) for k in pData.K), name = "durationConstr0")
    budgetConstr0 = mp.addConstr((sum(sum(pData.b[i][j]*x0[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B), name = "budgetConstr0")
    xConstr0 = mp.addConstrs((sum(x0[i,j] for j in pData.Ji[i]) <= 1 for i in pData.II), name = "xConstr0")
    
    durationConstr = mp.addConstrs((t[k[1]] - t[k[0]] >= pData.D[k[0]] + ed[k[0]]*G[k[0]] - \
                               sum(pData.D[k[0]]*pData.eff[k[0]][j]*x[k[0],j] + ed[k[0]]*pData.eff[k[0]][j]*s[k[0],j] for j in pData.Ji[k[0]]) \
                               for k in pData.K), name = "durationConstr")
    budgetConstr = mp.addConstr((sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B), name = "budgetConstr")
    xConstr = mp.addConstrs((sum(x[i,j] for j in pData.Ji[i]) <= 1 for i in pData.II), name = "xConstr")
    
    FConstr = mp.addConstrs((eH - (1 - G[i])*M <= t0[i] for i in pData.II), name = "FConstr")
    GConstr = mp.addConstrs((eH + G[i]*M >= t0[i] for i in pData.II), name = "GConstr")
    
    tConstr1 = mp.addConstrs((t[i] + G[i]*M >= t0[i] for i in pData.II), name = "tConstr1")
    tConstr2 = mp.addConstrs((t[i] - G[i]*M <= t0[i] for i in pData.II), name = "tConstr2")
    tConstr3 = mp.addConstrs((t[i] >= G[i]*eH for i in pData.II), name = "tConstr3")
    
    xConstr1 = mp.addConstrs((x[i,j] + G[i] >= x0[i,j] for i in pData.II for j in pData.Ji[i]), name = "xConstr1")
    xConstr2 = mp.addConstrs((x[i,j] - G[i] <= x0[i,j] for i in pData.II for j in pData.Ji[i]), name = "xConstr2")

    sLinear1 = mp.addConstrs((s[i,j] <= G[i] for i in pData.II for j in pData.Ji[i]), name = "sLinear1")
    sLinear2 = mp.addConstrs((s[i,j] <= x[i,j] for i in pData.II for j in pData.Ji[i]), name = "sLinear2")
    sLinear3 = mp.addConstrs((s[i,j] >= G[i] + x[i,j] - 1 for i in pData.II for j in pData.Ji[i]), name = "sLinear3")
    
    mp.setObjective(pData.p0*t0[0] + (1 - pData.p0)*t[0],GRB.MINIMIZE)
    mp.update()
    mp.optimize()
    
    texp = {}
    xexp = {}
    gexp = {}
    fexp = mp.objVal
    for i in pData.II:
        texp[i] = t0[i].X
        for j in pData.Ji[i]:
            xexp[i,j] = x0[i,j].X
        gexp[i] = G[i].X
    return texp,xexp,fexp,gexp,mp

    
def extForm(pData,disData,Omega,prec = 1e-4,TL = 36000):
    # extensive formulation, time limit 10 hours
    M = {}
    for omega in Omega:
        M[omega] = sum(max(pData.D[i],pData.D[i]+disData[omega].d[i]) for i in pData.II if i != 0)
        
    mp = Model('extensive_model')
    mp.Params.IntFeasTol = 1e-8
    mp.Params.TimeLimit = TL
    mp.Params.MIPGap = prec
    mp.Params.OutputFlag = 0
    
    HDiff = {}
    for omega in Omega:
        if disData[omega].H in HDiff.keys():
            HDiff[disData[omega].H].append(omega)
        else:
            HDiff[disData[omega].H] = [omega]
    HDict = {}
    HList = []
    for hIter in HDiff.keys():
        for item in range(len(HDiff[hIter])):
            HDict[HDiff[hIter][item]] = HDiff[hIter][0]
        HList.append(HDiff[hIter][0])
    HList = sorted(HList)
    
    t0 = mp.addVars(pData.II, lb = 0, name = "t0")
    JJ = [(i,j) for i in pData.II for j in pData.Ji[i]]
    x0 = mp.addVars(JJ, lb = 0, ub = 1, name = "x0")
    t = mp.addVars(pData.II, Omega, lb = 0, name = "t")
    x = mp.addVars(JJ, Omega, lb = 0, ub = 1, name = "x")
    G = mp.addVars(pData.II, HList, vtype=GRB.BINARY, name = "G")
    s = mp.addVars(JJ, Omega, lb = 0, ub = 1, name = "s")
    
    durationConstr0 = mp.addConstrs((t0[k[1]] - t0[k[0]] >= pData.D[k[0]]*(1 - sum(pData.eff[k[0]][j]*x0[k[0],j] for j in pData.Ji[k[0]])) for k in pData.K), name = "durationConstr0")
    budgetConstr0 = mp.addConstr((sum(sum(pData.b[i][j]*x0[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B), name = "budgetConstr0")
    xConstr0 = mp.addConstrs((sum(x0[i,j] for j in pData.Ji[i]) <= 1 for i in pData.II), name = "xConstr0")
    
    durationConstr = mp.addConstrs((t[k[1],omega] - t[k[0],omega] >= pData.D[k[0]] + disData[omega].d[k[0]]*G[k[0],HDict[omega]] - \
                                   sum(pData.D[k[0]]*pData.eff[k[0]][j]*x[k[0],j,omega] + disData[omega].d[k[0]]*pData.eff[k[0]][j]*s[k[0],j,omega] for j in pData.Ji[k[0]]) \
                                   for k in pData.K \
                                   for omega in Omega), name = "durationConstr")
    budgetConstr = mp.addConstrs((sum(sum(pData.b[i][j]*x[i,j,omega] for j in pData.Ji[i]) for i in pData.II) <= pData.B for omega in Omega), name = "budgetConstr")
    xConstr = mp.addConstrs((sum(x[i,j,omega] for j in pData.Ji[i]) <= 1 for i in pData.II for omega in Omega), name = "xConstr")
    
    FConstr = mp.addConstrs((disData[omega].H - (1 - G[i,HDict[omega]])*M[omega] <= t0[i] for i in pData.II for omega in Omega), name = "FConstr")
    GConstr = mp.addConstrs((disData[omega].H + G[i,HDict[omega]]*M[omega] >= t0[i] for i in pData.II for omega in Omega), name = "GConstr")
    
    tConstr1 = mp.addConstrs((t[i,omega] + G[i,HDict[omega]]*M[omega] >= t0[i] for i in pData.II for omega in Omega), name = "tConstr1")
    tConstr2 = mp.addConstrs((t[i,omega] - G[i,HDict[omega]]*M[omega] <= t0[i] for i in pData.II for omega in Omega), name = "tConstr2")
    tConstr3 = mp.addConstrs((t[i,omega] >= G[i,HDict[omega]]*disData[omega].H for i in pData.II for omega in Omega), name = "tConstr3")
    
    xConstr1 = mp.addConstrs((x[i,j,omega] + G[i,HDict[omega]] >= x0[i,j] for i in pData.II for j in pData.Ji[i] for omega in Omega), name = "xConstr1")
    xConstr2 = mp.addConstrs((x[i,j,omega] - G[i,HDict[omega]] <= x0[i,j] for i in pData.II for j in pData.Ji[i] for omega in Omega), name = "xConstr2")
    
    GGConstr1 = mp.addConstrs((G[i,HList[omega]] >= G[i,HList[omega+1]] for i in pData.II for omega in range(len(HList) - 1)), name = "GGConstr1")
    GGConstr2 = mp.addConstrs((G[k[1],omega] >= G[k[0],omega] for k in pData.K for omega in HList), name = "GGConstr2")
    
    sLinear1 = mp.addConstrs((s[i,j,omega] <= G[i,HDict[omega]] for i in pData.II for j in pData.Ji[i] for omega in Omega), name = "sLinear1")
    sLinear2 = mp.addConstrs((s[i,j,omega] <= x[i,j,omega] for i in pData.II for j in pData.Ji[i] for omega in Omega), name = "sLinear2")
    sLinear3 = mp.addConstrs((s[i,j,omega] >= G[i,HDict[omega]] + x[i,j,omega] - 1 for i in pData.II for j in pData.Ji[i] for omega in Omega), name = "sLinear3")
    
    mp.setObjective(pData.p0*t0[0] + sum(disData[omega].prDis*t[0,omega] for omega in Omega),GRB.MINIMIZE)
    mp.update()
    mp.optimize()
    
    text = {}
    xext = {}
    gext = {}
    fext = mp.objVal
    for i in pData.II:
        text[i] = t0[i].X
        for j in pData.Ji[i]:
            xext[i,j] = x0[i,j].X
        for omega in HList:
            gext[i,omega] = G[i,omega].X
    return text,xext,fext,gext,mp
    