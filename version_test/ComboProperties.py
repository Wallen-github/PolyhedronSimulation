# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: ComboProperties.py
@date: 3/28/23 14:21
@desc: 
"""
import numpy as np


def Get_CenterMass(coor,mass):
    CM = np.zeros([3])
    coor_cm = np.zeros([len(mass),3])
    for i in range(len(mass)):
        CM += mass[i]*coor[i,:]
    CM = CM/np.sum(mass)

    for i in range(len(mass)):
        coor_cm[i,:] = coor[i,:] - CM

    return CM,coor_cm



if __name__ == '__main__':

    rho = 2000
    coor = np.array([[1,0,0],[-1,0,0]])
    mass = np.array([8*rho,16*rho])
    CM,coor_cm = Get_CenterMass(coor,mass)
    print('CM = ', CM)
    print('coor_cm = ', coor_cm)