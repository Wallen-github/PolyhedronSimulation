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

def EarthAccel(GG,Mi,Posi,PosE):
    ME = 5.972E24  # kg, Earth mass
    Accel = np.zeros(6)
    Accel[0:3] = - GG*ME*Mi*(Posi - PosE)/np.linalg.norm(Posi - PosE)**3
    return Accel

def EarthPos(timespan,PosVec0,Unit,Gravorder):
    '''
    This function will compute the Earth acceleration based on the given Geocentric hyperbolic orbit and initial position/velocity
    :param timespan: sec, epoch span
    :param PosVec0: Initial position (m) and velocity (m/sec) of Earth
    :param Unit: unit length (m), mass (kg), and time (sec)
    :return PosVec: Final position (m) and velocity (m/sec) of Earth
    :return PosVecSol.y: Normalized positions and velocities during whole integration
    '''
    Lunit, Munit, Tunit = Unit
    PosVec0_norm = PosVec0.copy()
    PosVec0_norm[0:3] = PosVec0[0:3] / Lunit
    PosVec0_norm[3:6] = PosVec0[3:6] * Tunit / Lunit
    timespan = timespan / Tunit
    mu = 1
    tol = 1E-13
    PosVecSol = solve_ivp(fun=FlybyOrbit, t_span=timespan, y0=PosVec0_norm, args=(mu,Gravorder), method='RK45', rtol=tol, atol=tol)
    PosVec = PosVecSol.y[:,-1].copy()
    PosVec[0:3] = PosVec[0:3] * Lunit
    PosVec[3:6] = PosVec[3:6] * Lunit / Tunit

    return PosVec,PosVecSol.y

def InitialEarthPV(MA,GG=1.):
    '''
    This function initialized the Earth orbit around Apophis.
    :param MA: kg, Asteroid mass
    :param GG: N*m^2/kg^2, Gravitaional parameter, default value is 1.
    :return:
    PosVec0: Initial position (m) and velocity (m/sec) of Earth
    Unit: unit length (m), mass (kg), and time (sec)
    '''
    M0 = 5.972E24  # kg, Earth mass
    R0 = 6378.1370E3  # m, Earth radius
    q = 3.72E7  # m, Periapsis radius
    ecc = 4.229  # eccentricity
    axi = q / (ecc - 1)  # axis
    inc = 47.8  # deg, Inclination
    ome = -143.9  # deg, Argument of periapsis
    Ome = 1.8  # deg, Longitude of asc. node
    f = -100.  # deg, true anomaly
    ele_deg = np.array([axi, ecc, inc, ome, Ome, f])
    ele_rad = ele_deg.copy()
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

def FlybyOrbit(t,PosVec,mu,order):
    '''
    This function provide a EOM of two-body problem
    Args:
        t: time
        PosVec: position and velocity
        mu: mass, parameter, G * M

    Returns:
        Accel: acceleration
    '''

    Vterm = np.empty([2,3])
    Accel = np.zeros(6)

    # Zero order term
    Vterm[0,:] = - mu * PosVec[0:3] / np.linalg.norm(PosVec[0:3])**3

    # 2nd order term
    C20 = -1.082635E-3
    C21 = -2.58198E-10
    C22 = 1.57457E-6
    x,y,z = PosVec[0:3]
    Vterm[1,0] = (3*x*(C20*(x**2+y**2-4*z**2)+2*C22*(-3*x**2+7*y**2+2*z**2)))/(2*(x**2+y**2+z**2)**(7/2))
    Vterm[1,1] = (3*y*(C20*(x**2+y**2-4*z**2)-2*C22*(7*x**2-3*y**2+2*z**2)))/(2*(x**2+y**2+z**2)**(7/2))
    Vterm[1,2] = (3*z*(10*C22*(-x**2+y**2)+C20*(3*x**2+3*y**2-2*z**2)))/(2*(x**2+y**2+z**2)**(7/2))

    Accel[0:3] = PosVec[3:6]
    for i in range(order):
        Accel[3:6] = Accel[3:6] + Vterm[i,:]
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
    ele = np.array([1.,2.,0.,0.,0.,-80.])
    print('ele = ', ele,' deg')
    ele[2:6] = ele[2:6] * np.pi / 180  # degree to radian
    PV = Keplerian2hyperbola(mu,ele)
    print('PV = ',PV.flatten())

    Unit = np.array([1.,1.,1.])
    timespan = [0,10.]
    PV0 = PV.flatten()
    # PV0= np.array([1,0,0,0,0.5,0])
    Gravorder = 1
    PosVec, Sol = EarthPos(timespan, PV0, Unit, Gravorder)
    print('PosVec = ', PosVec)

    # Test integration with unit
    muA = 2.650 # m^3/sec^2
    timespan = [0,3*24*3600] # sec
    GG = 6.6742e-11 # G, N*m^2/kg^2
    MA = muA/GG # kg
    PosVec0, Unit = InitialEarthPV(MA,GG=GG)
    print('PosVec0 = ', PosVec0)
    Gravorder = 1
    PosVec,Sol = EarthPos(timespan,PosVec0,Unit, Gravorder)
    print('PosVec = ', PosVec)

    # Plot
    fig = plt.figure()
    ax1 = plt.axes(projection='3d')
    ax1.plot3D(Sol[0, :], Sol[1, :], Sol[2, :],c='black',label='geocentric hyperbolic orbit')  # 绘制空间曲线
    ax1.scatter3D(Sol[0, 0], Sol[1, 0], Sol[2, 0], c='red',s=20,label='initial position')
    ax1.scatter3D(0,0,0, c='blue',s=20,label='Earth')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('z')
    plt.legend()
    plt.show()

    # Compare different gravity order
    muA = 2.650  # m^3/sec^2
    timespan = [0, 3 * 24 * 3600]  # sec
    GG = 6.6742e-11  # G, N*m^2/kg^2
    MA = muA / GG  # kg
    PosVec0, Unit = InitialEarthPV(MA, GG=GG)
    print('PosVec0 = ', PosVec0)
    Gravorder = 1
    PosVec1, Sol1 = EarthPos(timespan, PosVec0, Unit, Gravorder)
    Gravorder = 2
    PosVec2, Sol2 = EarthPos(timespan, PosVec0, Unit, Gravorder)

    # Plot
    fig = plt.figure()
    ax1 = fig.add_subplot(3,2,1)
    ax1.plot(abs(Sol1[0,:]-Sol2[0, :]))
    ax1.set_ylabel('delta x')
    ax1.grid()
    ax1 = fig.add_subplot(3, 2, 3)
    ax1.plot(abs(Sol1[1, :] - Sol2[1, :]))
    ax1.set_ylabel('delta y')
    ax1.grid()
    ax1 = fig.add_subplot(3, 2, 5)
    ax1.plot(abs(Sol1[2, :] - Sol2[2, :]))
    ax1.set_ylabel('delta z')
    ax1.grid()
    ax1 = fig.add_subplot(3, 2, 2)
    ax1.plot(abs(Sol1[3, :] - Sol2[3, :]))
    ax1.set_ylabel('delta x_dot')
    ax1.grid()
    ax1 = fig.add_subplot(3, 2, 4)
    ax1.plot(abs(Sol1[4, :] - Sol2[4, :]))
    ax1.set_ylabel('delta y_dot')
    ax1.grid()
    ax1 = fig.add_subplot(3, 2, 6)
    ax1.plot(abs(Sol1[5, :] - Sol2[5, :]))
    ax1.set_ylabel('delta z_dot')
    ax1.grid()
    plt.show()



# See PyCharm help at https://www.jetbrains.com/help/pycharm/
