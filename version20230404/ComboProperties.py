# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: ComboProperties.py
@date: 3/28/23 14:21
@desc: 
"""
import numpy as np
from pylmgc90 import chipy


def Get_CenterMass(coor,mass):
    '''
    This function computes the center of mass
    Args:
        coor (numpy array): positions in the inertial frame
        mass (numpy array): mass

    Returns:
        CM center of mass
        coor_cm coordinate origin at center of mass
    '''
    CM = np.zeros([3])
    coor_cm = np.zeros([len(mass),3])
    for i in range(len(mass)):
        CM += mass[i]*coor[i,:]
    CM = CM/np.sum(mass)

    for i in range(len(mass)):
        coor_cm[i,:] = coor[i,:] - CM

    return CM,coor_cm

def Get_InitialVelocity(coor,mass,omega):
    '''
    This function computes initial velocity w.r.t center of mass
    Args:
        coor (numpy array): coordinate in the inertial frame
        mass (numpy array): mass
        omega (numpy array): combo's spin rate

    Returns:
        vel_cm: velocity w.r.t center of mass
        vel: velocity inertial origin

    '''
    vel_cm = np.zeros([len(mass), 3])
    vel = np.zeros([len(mass), 3])

    CM, coor_cm = Get_CenterMass(coor,mass)

    for i in range(len(mass)):
        vel_cm[i,:] = np.cross(omega, coor_cm[i,:])

    II = np.eye(len(mass))
    massmatrix = np.tile(mass[:1], (len(mass), 1))/np.sum(mass)
    massmatrix_inv = np.linalg.inv(II-massmatrix)
    vel = np.dot(massmatrix_inv,vel_cm)

    return vel_cm, vel

def Get_EnergyMomentum(nbR3,GG=1):

    momentum = np.empty([1,4], dtype=float)
    mass = np.empty([nbR3], dtype=float)
    coor = np.empty([nbR3, 6], dtype=float)
    vel = np.empty([nbR3, 6], dtype=float)
    vel_cm = np.zeros([3], dtype=float)
    Energy = 0
    Knetic = 0
    potent = 0

    for i in range(nbR3):
        coor[i,:] = chipy.RBDY3_GetBodyVector('Coorb', i + 1)
        mass[i] = chipy.RBDY3_GetMass(i + 1)
        vel[i,:] = chipy.RBDY3_GetBodyVector('Vfree', i + 1)

    Msum = np.sum(mass)
    for i in range(nbR3 - 1):
        for j in range(i + 1, nbR3):
            vij = vel[j, 0:3] - vel[i, 0:3]
            rij = coor[j, 0:3] - coor[i, 0:3]
            Knetic += mass[i] * mass[j] * np.dot(vij, vij) / (2. * Msum)
            momentum[0,0:3] += mass[i] * mass[j] * (np.cross(rij, vij)) / Msum
            potent += - GG * mass[i] * mass[j] / np.linalg.norm(rij) #gravpot[i]  #

    for i in range(nbR3):
        vel_cm += mass[i] * vel[i, 0:3] / Msum
        inertia = chipy.RBDY3_GetGlobInertia(i + 1)
        omega_i = vel[i, 3:6]
        Knetic += 0.5 * np.dot(omega_i,np.dot(inertia,omega_i))
        momentum[0,0:3] += np.dot(inertia,omega_i)
    Energy = Knetic + potent
    momentum[0,3] = np.dot(momentum[0,0:3],momentum[0,0:3])
    return Energy, Knetic, momentum

def Get_TotalInertia(nbR3,mass,coor_cm):
    total_inertia = np.zeros([3, 3], dtype=float)
    for i in range(nbR3):
        Rc = np.array([[0,-coor_cm[i,2],coor_cm[i,1]],
                       [coor_cm[i,2],0,-coor_cm[i,0]],
                       [-coor_cm[i,1],coor_cm[i,0],0]])
        inertia = chipy.RBDY3_GetGlobInertia(i + 1)
        inertia1 = chipy.RBDY3_GetBodyInertia(i + 1)
        total_inertia += inertia - mass[i]*np.dot(Rc,Rc)

    return total_inertia


def Get_EffectiveQuantity(kentic, momentum2, total_inertia):
    Inertia_t = np.nan
    omega_e = 2*kentic/np.sqrt(momentum2)
    inertia_d = momentum2/(2*kentic)
    diagInertia = np.linalg.eig(total_inertia)
    Il = min(diagInertia[0])
    Is = max(diagInertia[0])
    Ii = np.median(diagInertia[0])
    if Ii<inertia_d<Is:
        Inertia_t = (inertia_d - Ii)/(Is-Ii)
    elif Il < inertia_d < Ii:
        Inertia_t = (inertia_d - Ii)/(Ii - Il)

    return omega_e, inertia_d, Inertia_t

def Shift2CM(vertices,CM):
    '''
    This function convert polyhedron combo from inertial coordinate to center of mass coordinate
    Args:
        vertices (numpy array n by 3): vertices of one polyhedron
        CM (numpy array 1 by 3): center of mass of combo

    Returns:
        vertices_cm: vertices coordinate w.r.t center of mass

    '''
    vertices_cm = np.empty([len(vertices), 3], dtype=float)
    for i in range(len(vertices)):
        vertices_cm[i,:] = vertices[i,:] - CM
    return vertices_cm





if __name__ == '__main__':

    rho = 2000
    coor = np.array([[1,0,0],[-2,0,0]])
    mass = np.array([8*rho,16*rho])
    CM,coor_cm = Get_CenterMass(coor,mass)
    print('CM = ', CM)
    print('coor_cm = ', coor_cm)

    omega = np.array([0,0,1])
    vel_cm, vel = Get_InitialVelocity(coor, mass, omega)
    print('vel_cm = ', vel_cm)
    print('vel = ', vel)


    vel_test = np.zeros([len(mass), 3])
    for i in range(len(mass)):
        vel_test[i,:] = np.cross(omega, coor[i,:])

    print('vel_test = ', vel_test)