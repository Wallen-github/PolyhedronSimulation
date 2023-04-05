# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: ShapeModel.py
@date: 3/29/23 13:07
@desc: 
"""
import numpy as np
from pylmgc90.pre import *
from ComboProperties import *
from FlybyOrbit import InitialEarthPV


def UnitPoly(mat,dim,vertices1,faces1,vertices2,faces2,rho,volume_cos,shiftpos_cos):
    stone = material(name='STONE', materialType='RIGID', density=rho)
    mat.addMaterial(stone)
    mod = model(name='rigid', physics='MECAx', element='Rxx3D', dimension=dim)

    poly1 = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
                            generation_type='full', vertices=vertices1,
                            faces=faces1, tol=0., number=None, seed=None,
                            xr=1., yr=1., zr=1.)
    poly2 = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
                            generation_type='full', vertices=vertices2,
                            faces=faces2, tol=0., number=None, seed=None,
                            xr=1., yr=1., zr=1.)
    coor = np.vstack((poly1.nodes[1].coor.copy(), poly2.nodes[1].coor.copy()))
    mass = np.array([poly1.contactors[0].volume * rho, poly2.contactors[0].volume * rho])
    CM, coor_cm = Get_CenterMass(coor, mass)

    Vsum = poly1.contactors[0].volume + poly2.contactors[0].volume
    Lunit = (Vsum) ** (1 / 3)
    vertices3 = Shift2CM(vertices1, CM) / Lunit * (volume_cos) ** (1 / 3) + shiftpos_cos
    vertices4 = Shift2CM(vertices2, CM) / Lunit * (volume_cos) ** (1 / 3) + shiftpos_cos

    poly3 = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
                            generation_type='full', vertices=vertices3,
                            faces=faces1, tol=0., number=None, seed=None,
                            xr=1., yr=1., zr=1.)
    poly4 = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
                            generation_type='full', vertices=vertices4,
                            faces=faces2, tol=0., number=None, seed=None,
                            xr=1., yr=1., zr=1.)
    return poly3,poly4

def InitializeVelocity(poly1,poly2,omega_comb,rho,shiftvel_cos):
    coor = np.vstack((poly1.nodes[1].coor.copy(), poly2.nodes[1].coor.copy()))
    mass = np.array([poly1.contactors[0].volume * rho, poly2.contactors[0].volume * rho])
    vel_cm, vel_in = Get_InitialVelocity(coor, mass, omega_comb)

    vel = np.empty([len(vel_cm),3], dtype=float)
    for i in range(len(vel_cm)):
        vel[i,:] = vel_cm[i,:] + shiftvel_cos
    # vel = vel_cm
    velocity1 = list(np.concatenate((vel[0, :], omega_comb)))
    velocity2 = list(np.concatenate((vel[1, :], omega_comb)))
    # poly1.imposeDrivenDof(component=[1, 2, 3], dofty='vlocy')
    # poly2.imposeDrivenDof(component=[1, 2, 3], dofty='vlocy')
    poly1.imposeInitValue(component=[1, 2, 3, 4, 5, 6], value=velocity1)
    poly2.imposeInitValue(component=[1, 2, 3, 4, 5, 6], value=velocity2)

def SetAttitude(vertices,lam,bet):

    RM1 = RotationMatrix(-23*np.pi/180, 'Z')
    RM2 = RotationMatrix(149*np.pi/180, 'X')
    RM3 = RotationMatrix(0*np.pi/180, 'Z')
    DCM = np.dot(RM3,np.dot(RM2,RM1))

    if len(vertices.shape)==2:
        vertices_new = np.zeros([vertices.shape[0],3])
        for i in range(vertices.shape[0]):
            vertices_new[i,:] = np.dot(DCM,vertices[i,:])
    elif len(vertices.shape)==1:
        vertices_new = np.dot(DCM, vertices)
    return vertices_new

def SetSpinRate():
    Lunit, Munit, Tunit = SetUnit(GG=6.6742e-11, volume_cos=1.986E7, rho=2000)
    Pl = 30.6 * 3600 / Tunit
    PPhi = 265.7 * 3600 / Tunit
    omegal = np.pi * 2 / Pl
    omegaPhi = np.pi * 2 / PPhi
    omega_comb = np.array([omegal, 0, omegaPhi])
    return omega_comb


def SetUnit(GG,volume_cos,rho):
    Munit = volume_cos * rho
    Lunit = (volume_cos) ** (1 / 3)
    Tunit = np.sqrt(Lunit ** 3 / GG / Munit)
    return Lunit, Munit, Tunit


def SetUnit_Apophis():
    GG = 6.6742e-11  # G, N*m^2/kg^2
    Volume_apophis = 1.986E7  # m^3
    Lunit = (Volume_apophis) ** (1 / 3)
    Munit = 2000 * Volume_apophis
    Tunit = np.sqrt(Lunit ** 3 / GG / Munit)
    return Lunit,Munit,Tunit

def RotationMatrix(theta,axis):
    if axis == 'X':
        RotationX = np.array([[1, 0, 0],
                              [0, np.cos(theta), np.sin(theta)],
                              [0, -np.sin(theta), np.cos(theta)]])
        return RotationX
    elif axis == 'Y':
        RotationY = np.array([[np.cos(theta), 0, -np.sin(theta)],
                              [0, 1, 0],
                              [np.sin(theta), 0, np.cos(theta)]])
        return RotationY
    elif axis == 'Z':
        RotationZ = np.array([[np.cos(theta), np.sin(theta), 0],
                              [-np.sin(theta), np.cos(theta), 0],
                              [0, 0, 1]])
        return RotationZ
    else:
        print('error in func RotationMatrix: the "axis" should be one of X,Y,Z.')


if __name__ == '__main__':
    bodies = avatars()
    mat = materials()
    svs = see_tables()
    tacts = tact_behavs()

    dim = 3
    rho = 1
    volume_cos = 1
    GG = 1
    Lunit, Munit, Tunit = SetUnit(GG=6.6742e-11, volume_cos=1.986E7, rho=2000)
    PosVec0, UnitFlyby = InitialEarthPV(Munit, GG=6.6742e-11)

    vertices1 = np.array([[2, 1, -1], [2, -1, -1], [2, -1, 1], [2, 1, 1],
                          [0, 1, -1], [0, -1, -1], [0, -1, 1], [0, 1, 1]])
    faces1 = np.array([[1, 2, 3], [1, 4, 3], [1, 4, 8], [1, 5, 8], [1, 2, 6], [1, 5, 6],
                       [7, 8, 4], [7, 3, 4], [7, 3, 2], [7, 6, 2], [7, 8, 5], [7, 6, 5]])
    vertices2 = np.array([[0, 2, -2], [-3, 2, -2], [-3, 2, 2], [0, 2, 2],
                          [0, -2, -2], [-3, -2, -2], [-3, -2, 2], [0, -2, 2]])
    faces2 = np.array([[1, 2, 3], [1, 4, 3], [1, 4, 8], [1, 5, 8], [1, 2, 6], [1, 5, 6],
                       [7, 8, 4], [7, 3, 4], [7, 3, 2], [7, 2, 6], [7, 8, 5], [7, 6, 5]])

    lam = np.pi/2
    bet = np.pi/4
    vertices_new1 = SetAttitude(vertices1, lam, bet)
    vertices_new2 = SetAttitude(vertices2, lam, bet)

    shiftpos_cos = np.array([0,0,0]) #PosVec0[0:3] / Lunit #
    poly1,poly2 = UnitPoly(mat, dim, vertices_new1, faces1, vertices_new2, faces2, rho, volume_cos,shiftpos_cos)

    shiftvel_cos = np.array([0,0,0]) #PosVec0[3:6] * Tunit / Lunit #
    SpinPeriod = 30.6*3600/Tunit
    omega_comb = SetSpinRate() #np.array([0, 0, 2*np.pi/SpinPeriod])
    omega_comb = SetAttitude(omega_comb, lam, bet)
    InitializeVelocity(poly1, poly2, omega_comb, rho,shiftvel_cos)

    bodies.addAvatar(poly1)
    bodies.addAvatar(poly2)

    LawSPSPx = tact_behav(name='rst01', law='RST_CLB', rstn=1.0, rstt=0.0, fric=0)

    tacts += LawSPSPx

    svSPSPx = see_table(CorpsCandidat='RBDY3', candidat='POLYR', colorCandidat='BLEUx', behav=LawSPSPx,
                        CorpsAntagoniste='RBDY3', antagoniste='POLYR', colorAntagoniste='BLEUx', alert=0.)
    svs += svSPSPx

    writeBodies(bodies, chemin='DATBOX/')
    writeBulkBehav(mat, chemin='DATBOX/', dim=dim, gravy=[0., 0., 0.])
    writeTactBehav(tacts, svs, chemin='DATBOX/')
    writeDrvDof(bodies, chemin='DATBOX/')
    writeDofIni(bodies, chemin='DATBOX/')
    writeVlocRlocIni(chemin='DATBOX/')

    print('volume #1 = ',poly1.contactors[0].volume)
    print('volume #2 = ',poly2.contactors[0].volume)
    print('Flyby initial pos = ', PosVec0[0:3],' m')
    print('Flyby initial vel = ', PosVec0[3:6],' m/sec')
    print('Flyby initial pos = ', shiftpos_cos)
    print('Flyby initial vel = ', shiftvel_cos)

    # try:
    #   visuAvatars(bodies)
    # except:
    #   pass

