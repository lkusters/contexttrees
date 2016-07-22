# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 17:28:45 2016

@author: Lieneke Kusters
"""

import pickle

filenamein = 'result.pi'
with open(filenamein, 'rb') as file:
    result = pickle.load(file)
    
tree = result['tree']