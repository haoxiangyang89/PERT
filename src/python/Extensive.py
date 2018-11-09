#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 22:52:07 2018

@author: haoxiangyang
"""
from gurobipy import *

# solving the extensive formulation
def extForm_cheat(pData,disData,Omega,prec = 1e-4,TL = 99999999):
    M = {}
    for omega in Omega:
        M[omega] = sum(max(pData.D[i],pData.D[i]+disData[omega].d[i]) for i in pData.II if i != 0)

    mp = Model('extensive')
    mp.Params.TimeLimit = TL
    # mp = Model(solver = CbcSolver());
    t0 = mp.addVars(pData.II)
    xInd = []
    for i in pData.II:
        for j in pData.Ji[i]:
            xInd.append((i,j))
    x0 = mp.addVars(xInd, ub = 1)
    t = mp.addVars(pData.II, Omega)
    x = mp.addVars(xInd, Omega, ub = 1)
    G = mp.addVars(pData.II, Omega, vtype = GRB.BINARY)
    s = mp.addVars(xInd, Omega, ub = 1)

    FConstr = mp.addConstrs(disData[omega].H - (1 - G[i,omega])*M[omega] <= t0[i] for i in pData.II
                                                                                  for j in pData.Ji[i]);
    GConstr = mp.addConstrs(disData[omega].H + G[i,omega]*M[omega] >= t0[i] for i in pData.II for omega in Omega)
    tConstr1 = mp.addConstrs(t[i,omega] + G[i,omega]*M[omega] >= t0[i] for i in pData.II for omega in Omega)
    tConstr2 = mp.addConstrs(t[i,omega] - G[i,omega]*M[omega] <= t0[i] for i in pData.II for omega in Omega)
    tConstr3 = mp.addConstrs(t[i,omega] >= disData[omega].H * G[i,omega] for i in pData.II for omega in Omega)
    xConstr1 = mp.addConstrs(x[i,j,omega] + G[i,omega] >= x0[i,j] for i in pData.II for j in pData.Ji[i] for omega in Omega)
    xConstr2 = mp.addConstrs(x[i,j,omega] - G[i,omega] <= x0[i,j] for i in pData.II for j in pData.Ji[i] for omega in Omega)
    durationConstr1 = mp.addConstrs(t[k[1],omega] - t[k[0],omega] >= pData.D[k[0]] + disData[omega].d[k[0]]*G[k[0],omega]
                      - sum(pData.D[k[0]]*pData.eff[k[0]][j]*x[k[0],j,omega] + disData[omega].d[k[0]]*pData.eff[k[0]][j]*s[k[0],j,omega] for j in pData.Ji[k[0]])
                      for k in pData.K for omega in Omega)
    durationConstr2 = mp.addConstrs(t0[k[1]] - t0[k[0]] >= pData.D[k[0]]*(1 - sum(pData.eff[k[0]][j]*x0[k[0],j] for j in pData.Ji[k[0]])) for k in pData.K)
    xConstr = mp.addConstrs(sum(x[i,j,omega] for j in pData.Ji[i]) <= 1 for i in pData.II for omega in Omega)
    GGcons1 = mp.addConstrs(G[i,omega] >= G[i,omega + 1] for i in pData.II for omega in range(len(Omega) - 1))
    GGcons2 = mp.addConstrs(G[k[1],omega] >= G[k[0],omega] for k in pData.K for omega in Omega)
    budgetConstr = mp.addConstrs(sum(sum(pData.b[i][j]*x[i,j,omega] for j in pData.Ji[i]) for i in pData.II) <= pData.B for omega in Omega)
    xConstr0 = mp.addConstrs(sum(x0[i,j] for j in pData.Ji[i]) <= 1 for i in pData.II)
    budgetConstr0 = mp.addConstr(sum(sum(pData.b[i][j]*x0[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B)
    Slinear1 = mp.addConstrs(s[i,j,omega] <= G[i,omega] for i in pData.II for j in pData.Ji[i] for omega in Omega)
    Slinear2 = mp.addConstrs(s[i,j,omega] <= x[i,j,omega] for i in pData.II for j in pData.Ji[i] for omega in Omega)
    Slinear3 = mp.addConstrs(s[i,j,omega] >= x[i,j,omega] - 1 + G[i,omega] for i in pData.II for j in pData.Ji[i] for omega in Omega)

    mp.setObjective(pData.p0*t0[0] + sum(disData[omega].prDis*t[0,omega] for omega in Omega), GRB.MINIMIZE)

    mp.optimize()

    text = {}
    xext = {}
    gext = {}
    for i in pData.II:
        text[i] = t0[i].X;
        for j in pData.Ji[i]:
            xext[i,j] = x0[i,j].X
        for omega in Omega:
            gext[i,omega] = G[i,omega].X
    fext = mp.ObjVal
    return text,xext,fext,gext,mp