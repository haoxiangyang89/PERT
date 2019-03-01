#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 14:18:42 2019

@author: haoxiangyang
"""

def mycallback(model, where):
    if where == GRB.Callback.SIMPLEX:
      print(model.cbGet(GRB.Callback.SPX_OBJVAL))