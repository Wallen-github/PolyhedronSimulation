# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: ShapeModel.py
@date: 3/8/23 15:44
@desc: 
"""

import numpy as np
from pylmgc90.pre import *
from scipy.optimize import fsolve

def ShapeModel_customize1():

    rho = 2000  # kg/m^3

    mat = material(name='STONE', materialType='RIGID', density=rho)
    mod = model(name='rigid', physics='MECAx', element='Rxx3D', dimension=3)


    return mat,mod,poly1,poly2


def ShapeModel_regular1():
    '''
    This Shape model is made by two regular polygons generated from spheres with radii 'R' and 'r'. The two radiis of
    sphere can be solved from known parameter 'volume' and 'largest extent'. For Apophis, the 'volume = 1.986E7 m^3',
    the 'largest extent = 410 m'. The 'largest extent' also can be subtituted by the 'mean radius = 168 m'.
    :return: Polygon avatar
    '''

    volume = 1.986E7 # m^3
    meanradii = 410. # m
    radii = fsolve(RadiiFunc,[50,50],args = [volume,meanradii])
    print(radii)
    print(radii[radii != max(radii)])

    rho = 2000 # kg/m^3
    radii1 = max(radii)
    radii2 = min(radii)
    pos = np.array([[radii1,0,0],[-radii2,0,0]])
    center = CenterFunc(pos, [radii1,radii2], 2000)
    print(center)
    center1 = center[0, :]
    center2 = center[1, :]
    nbvertices1 = 20
    nbvertices2 = 10

    mat = material(name='STONE', materialType='RIGID', density=rho)
    mod = model(name='rigid', physics='MECAx', element='Rxx3D', dimension=3)

    poly1 = rigidPolyhedron(model=mod, material=mat, center=center1, color='BLEUx',
                           generation_type='regular', nb_vertices=nbvertices1, vertices=None,
                           faces=None, radius=radii1, tol=0., number=None, seed=None,
                        xr=1., yr=1., zr=1.)
    poly2 = rigidPolyhedron(model=mod, material=mat, center=center2, color='BLEUx',
                           generation_type='regular', nb_vertices=nbvertices2, vertices=None,
                           faces=None, radius=radii2, tol=0., number=None, seed=None,
                        xr=1., yr=1., zr=1.)

    # Compute the initial spin rotation vector
    SpinPeriod = 30.6*3600 # sec
    omegaScalar = 2*np.pi/SpinPeriod # 1/sec
    Lati = -59.3*np.pi/180 # rad
    Long = 246.8*np.pi/180 # rad
    omegaVec = np.array([[omegaScalar*np.cos(Lati)*np.cos(Long)],
                         [omegaScalar*np.cos(Lati)*np.sin(Long)],
                         [omegaScalar*np.sin(Long)]])
    # poly2.imposeInitValue(component=[4,5,6], value=[0,0,2*np.pi/10000])
    # poly1.imposeDrivenDof(component=[1, 2, 3], dofty='vlocy')
    # poly2.imposeDrivenDof(component=[1, 2, 3], dofty='vlocy')
    poly1.imposeInitValue(component=[4,5,6], value=[omegaVec[0],omegaVec[1],omegaVec[2]])
    poly2.imposeInitValue(component=[4, 5, 6], value=[omegaVec[0], omegaVec[1], omegaVec[2]])

    volume = 0.
    # pour chaque contacteur
    for tact in poly1.contactors:
        # on ajoute la contribution du contacteur courant
        volume += tact.volume
    for tact in poly2.contactors:
        # on ajoute la contribution du contacteur courant
        volume += tact.volume
    print(volume)
    print(np.pi * 4 / 3 * radii[0] ** 3 + np.pi * 4 / 3 * radii[1] ** 3)

    return mat,mod,poly1,poly2

def RadiiFunc(radii,par):
    '''
    This function provide a equation set for solving radius based on known volume and mean radius.
    :param radii: radius of two particles
    :param par: volume and largest mean radius
    :return:
    '''
    volume = par[0]
    meanradii = par[1]
    Eq1 = -volume
    Eq2 = -meanradii
    for i in range(len(radii)):
        Eq1 = Eq1 + np.pi*4/3*radii[i]**3
        Eq2 = Eq2 + radii[i]
    Eq3 = radii[0] - 2/3 * radii[1]
    if len(radii)==3:
        return [Eq1, Eq2, Eq3]
    elif len(radii)==2:
        return [Eq2, Eq3]
    else:
        return 'error in RadiiFunc'
def CenterFunc(pos,radii,rho):
    '''
    This function will correct the position of particles into the center mass frame
    :param pos: position of particles; array([n,3])
    :param radii: radius of particles; array(n)
    :param rho: density of particles, all particles have same density
    :return: corrected positions
    '''
    mass = np.zeros([len(radii)])
    masscenter = np.zeros(3)
    corpos = np.zeros([len(radii),3])
    for i in range(0,len(radii)):
        mass[i] = 4.0*np.pi/3.0 * rho * radii[i]**3
        masscenter = masscenter + mass[i]*pos[i,:]

    MassTotal = np.sum(mass)
    masscenter = masscenter/MassTotal
    for i in range(0,len(radii)):
        corpos[i,:] = pos[i,:] - masscenter

    return corpos


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    [mat,mod,poly1,poly2] = ShapeModel_regular1()
