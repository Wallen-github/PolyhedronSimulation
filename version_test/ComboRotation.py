# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: ComboRotationY.py
@date: 3/28/23 13:42
@desc: 
"""
import numpy as np


def ComboRotstaionY(vertices,theta):

    rotVer = np.zeros([len(vertices),3])
    rotY = np.array([[np.cos(theta), 0, np.sin(theta)],
                     [0, 1, 0],
                     [-np.sin(theta), 0, np.cos(theta)]])
    for i in range(len(vertices)):
        rotVer[i,:] = np.dot(rotY,vertices[i,:])

    return rotVer



if __name__ == '__main__':
    vertices1 = np.array([[2, 1, -1], [2, -1, -1], [2, -1, 1], [2, 1, 1],
                          [0, 1, -1], [0, -1, -1], [0, -1, 1], [0, 1, 1]])

    theta = -np.pi/2
    rotVer = ComboRotstaionY(vertices1,theta)
    print(rotVer)
