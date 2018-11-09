#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 22:56:46 2018

@author: haoxiangyang
"""

# define the variable type
class pInfo:
    def __init__(self, II, Ji, D, b, eff, B, p0, K, Pre, Succ):
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
        
class disInfo:
    def __init__(self, H, d, prDis):
        self.H = H
        self.d = d
        self.prDis = prDis
        
class partType:
    def __init__(self, startH, endH):
        self.startH = startH
        self.endH = endH