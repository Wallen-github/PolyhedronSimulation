# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: EarthGravity.py
@date: 3/8/23 18:35
@desc: 
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import *

def EarthAccel(Posi,PosE,GG):
    ME = 5.972E24  # kg, Earth mass
    r0i = PosE - Posi
    Accel = - GG*ME/np.linalg.norm(r0i) * r0i - GG*ME/np.linalg.norm(PosE) * PosE
    return Accel

def EarthPos(time,PosVec0,Unit):
    '''
    This function will compute the Earth acceleration based on the given Geocentric hyperbolic orbit and initial position/velocity
    :param time: sec, epoch
    :param MA: kg, Mass of asteroid Assumes density of 2 g∕cm^3
    :param GG: N*m^2/kg^2, gravitational constant, 6.6742e-11
    :return PosVec: position (m) and velocity (m/sec) of Earth
    '''
    Lunit, Munit, Tunit = Unit
    PosVec0[0:3] = PosVec0[0:3]/Lunit
    PosVec0[3:6] = PosVec0[3:6] * Tunit / Lunit
    timespan = [0, time/Tunit]
    mu = 1
    PosVecSol = solve_ivp(fun=FlybyOrbit,t_span=timespan,y0=PosVec0,args=(mu,))
    PosVec = PosVecSol.y[:,-1]
    PosVec[0:3] = PosVec0[0:3] * Lunit
    PosVec[3:6] = PosVec0[3:6] * Lunit / Tunit

    return PosVec,PosVecSol.y

def InitialEarthPV(MA,GG=1.):
    M0 = 5.972E24  # kg, Earth mass
    R0 = 6378.1370E3  # m, Earth radius
    q = 3.72E7  # m, Periapsis radius
    ecc = 4.229  # eccentricity
    axi = q / (ecc - 1)  # axis
    inc = 47.8  # deg, Inclination
    ome = -143.9  # deg, Argument of periapsis
    Ome = 1.8  # deg, Longitude of asc. node
    f = -10.  # deg, true anomaly
    ele_deg = np.array([axi, ecc, inc, ome, Ome, f])
    ele_rad = ele_deg
    ele_rad[2:6] = ele_deg[2:6] * np.pi / 180  # degree to radian
    mu = GG * (M0 + MA)
    PosVec0 = Keplerian2hyperbola(mu, ele_rad)
    PosVec0 = PosVec0.flatten()

    # Unit normalization
    Lunit = R0
    Munit = MA + M0
    Tunit = np.sqrt(Lunit ** 3 / (GG * Munit))
    Unit = np.array([Lunit, Munit, Tunit])

    return PosVec0, Unit

def FlybyOrbit(t,PosVec,mu):
    # [GG,M0,Mi] = parameters
    Accel = - mu / np.linalg.norm(PosVec) * PosVec
    return Accel

def Keplerian2hyperbola(mu,ele):
    '''
    This function transport the hyperbola orbital elements to position and velocity in inertial frame
    :param mu: G*(m1+m2) gravitational parameter
    :param ele: orbital elements
    :return: position and velocity
    '''
    [axi,ecc,inc,ome,Ome,f] = ele

    fboundary = np.arccos(-1./ecc)
    if f>fboundary or f< -fboundary:
        print('Error in "Keplerian2hyperbola": The true anomaly exceed the boundary f_boundary = ',
              fboundary*180/np.pi, 'deg')

    p = axi * (ecc ** 2 - 1)

    r = abs(p / (1 + ecc * np.cos(f)))

    Phat = np.array([[np.cos(Ome)*np.cos(ome) - np.sin(Ome)*np.sin(ome)*np.cos(inc)],
                     [np.sin(Ome)*np.cos(ome) + np.cos(Ome)*np.sin(ome)*np.cos(inc)],
                     [np.sin(ome)*np.sin(inc)]])
    Qhat = np.array([[-np.cos(Ome)*np.sin(ome) - np.sin(Ome)*np.cos(ome)*np.cos(inc)],
                     [-np.sin(Ome)*np.sin(ome) + np.cos(Ome)*np.cos(ome)*np.cos(inc)],
                     [np.cos(ome)*np.sin(inc)]])

    pos = r*np.cos(f)*Phat + r*np.sin(f)*Qhat

    vel = np.sqrt(mu/p)*(-np.sin(f)*Phat + (np.cos(f)+ecc)*Qhat)

    PV = np.append(pos,vel,axis=0)

    return PV

# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    # Test: Keplerian2hyperbola
    mu = 1
    ele = np.array([1.,2.,0.,0.,0.,110.])
    print('ele = ', ele,' deg')
    ele[2:6] = ele[2:6] * np.pi / 180  # degree to radian
    PV = Keplerian2hyperbola(mu,ele)
    print('PV = ',PV.flatten())

    muA = 2.650
    time = 100000
    GG = 6.6742e-11
    MA = muA/GG
    PosVec0, Unit = InitialEarthPV(MA,GG=GG)
    print('PosVec0 = ', PosVec0)
    PosVec,Sol = EarthPos(time,PosVec0,Unit)
    print('PosVec = ', PosVec)

    # 定义坐标轴
    fig = plt.figure()
    ax1 = plt.axes(projection='3d')
    ax1.plot3D(Sol[0, :], Sol[1, :], Sol[2, :], 'blue')  # 绘制空间曲线
    ax1.scatter3D(Sol[0, :], Sol[1, :], Sol[2, :], 'red')
    plt.show()


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
