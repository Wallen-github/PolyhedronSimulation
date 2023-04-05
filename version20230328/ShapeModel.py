# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: ShapeModel.py
@date: 3/29/23 13:07
@desc: 
"""

from pylmgc90.pre import *

from version20230328.ComboProperties import *
from version20230328.FlybyOrbit import InitialEarthPV


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

def SetUnit():
    GG = 6.6742e-11  # G, N*m^2/kg^2
    Volume_apophis = 1.986E7  # m^3
    Lunit = (Volume_apophis) ** (1 / 3)
    Munit = 2000 * Volume_apophis
    Tunit = np.sqrt(Lunit ** 3 / GG / Munit)
    return Lunit,Munit,Tunit

if __name__ == '__main__':
    bodies = avatars()
    mat = materials()
    svs = see_tables()
    tacts = tact_behavs()

    GG = 6.6742e-11  # G, N*m^2/kg^2
    Lunit, Munit, Tunit = SetUnit()
    PosVecE0, UnitFly = InitialEarthPV(MA=Munit, GG=GG)  # output with unit
    FlybyInitialState = PosVecE0.copy()

    dim = 3
    rho = 1
    volume_cos = 1
    vertices1 = np.array([[2, 1, -1], [2, -1, -1], [2, -1, 1], [2, 1, 1],
                          [0, 1, -1], [0, -1, -1], [0, -1, 1], [0, 1, 1]])
    faces1 = np.array([[1, 2, 3], [1, 4, 3], [1, 4, 8], [1, 5, 8], [1, 2, 6], [1, 5, 6],
                       [7, 8, 4], [7, 3, 4], [7, 3, 2], [7, 6, 2], [7, 8, 5], [7, 6, 5]])
    vertices2 = np.array([[0, 2, -2], [-3, 2, -2], [-3, 2, 1], [0, 2, 1],
                          [0, -2, -2], [-3, -2, -2], [-3, -2, 1], [0, -2, 1]])
    faces2 = np.array([[1, 2, 3], [1, 4, 3], [1, 4, 8], [1, 5, 8], [1, 2, 6], [1, 5, 6],
                       [7, 8, 4], [7, 3, 4], [7, 3, 2], [7, 2, 6], [7, 8, 5], [7, 6, 5]])
    shiftpos_cos = np.zeros([3]) #FlybyInitialState[0:3] / Lunit
    poly1,poly2 = UnitPoly(mat, dim, vertices1, faces1, vertices2, faces2, rho, volume_cos,shiftpos_cos)

    shiftvel_cos = np.zeros([3]) #FlybyInitialState[3:6] * Tunit / Lunit
    omega_comb = np.array([0, 0, 2*np.pi*Tunit/(30.6*60*60)])
    InitializeVelocity(poly1, poly2, omega_comb, rho,shiftvel_cos)

    bodies.addAvatar(poly1)
    bodies.addAvatar(poly2)

    LawSPSPx = tact_behav(name='rst01', law='RST_CLB', rstn=1.0, rstt=0.0, fric=1.)

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

    try:
      visuAvatars(bodies)
    except:
      pass

