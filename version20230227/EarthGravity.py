# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: EarthGravity.py
@date: 3/8/23 18:35
@desc: 
"""
import numpy as np


def EarthAccel(time,muA):
    '''
    This function will compute the Earth acceleration based on the given Geocentric hyperbolic orbit
    :param time: epoch
    :param muA: Gravitational parameter of asteroid Assumes density of 2 g∕cm^3
    :return: Earth acceleration
    '''
    Accel = np.zeros(1,3)
    q = 3.72E7 # m, Periapsis radius
    ecc = 4.229 # eccentricity
    axi = q/(ecc-1) # axis
    inc = 47.8 # deg, Inclination
    ome = -143.9 # deg, Argument of periapsis
    Ome = 1.8 # deg, Longitude of asc. node
    f = 90 # deg, true anomaly
    ele = np.array([axi,ecc,inc,ome,Ome,f])
    PosVec = FlybyOrbit(ele)
    M0 = 5.972E24 # kg Earth mass
    Accel = - M0*muA/np.norm(PosVec)*PosVec

    return Accel

def FlybyOrbit(ele):
    PosVec = np.zeros(1,3)

    return PosVec


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    muA = 2.650
    time = 10
    Accel = EarthAccel(time,muA)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
