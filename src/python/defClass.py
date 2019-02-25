#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 15:16:10 2019

@author: haoxiangyang
"""

# definition of variables/structures
# activity network information
class pInfo:
    def __init__(self,II,Ji,D,b,eff,B,p0,K,Pre,Succ):
        self.II = II
        self.Ji = Ji
        self.D = D
        self.b = b
        self.eff = eff
        self.B = B
        self.p0 = p0
        self.K = K
        self.Pre = Pre
        self.Succ = Succ
        
# disruption information
class disInfo:
    def __init__(self,H,d,prDis):
        self.H = H
        self.d = d
        self.prDis = prDis
        