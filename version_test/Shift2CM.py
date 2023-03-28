# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: Shift2CM.py
@date: 3/26/23 23:48
@desc: 
"""

import numpy as np

def Shift2CM(coor,mass):
    OG = 0
    p_coor = np.empty([len(mass), 6], dtype=float)
    for i in range(len(mass)):
        OG += mass[i]*coor[i,:]
    OG = OG/np.sum(mass)

    for i in range(len(mass)):
        p_coor[i,:] = coor[i,:] - OG
    return p_coor
