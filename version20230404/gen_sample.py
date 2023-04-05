# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: gen_sample.py
@date: 1/10/23 13:38
@desc: 
"""
import numpy as np
from pylmgc90 import chipy
from pylmgc90.pre import *

from version20230328.ComboProperties import *

if not os.path.isdir('./DATBOX'):
    os.mkdir('./DATBOX')

if 'norand' in sys.argv:
    seed = list(range(13))
else:
    seed = None

dim = 3

bodies = avatars()
mat    = materials()
svs    = see_tables()
tacts  = tact_behavs()

rho = 2000
stone = material(name='STONE', materialType='RIGID', density=rho)
mat.addMaterial(stone)
mod = model(name='rigid', physics='MECAx', element='Rxx3D', dimension=dim)

displacement = np.array([[0,0,0]])

vertices1 = np.array([[2,1,-1],[2,-1,-1],[2,-1,1],[2,1,1],
                      [0,1,-1],[0,-1,-1],[0,-1,1],[0,1,1]]) + displacement
faces1 = np.array([[1,2,3],[1,4,3],[1,4,8],[1,5,8],[1,2,6],[1,5,6],
                   [7,8,4],[7,3,4],[7,3,2],[7,6,2],[7,8,5],[7,6,5]])
vertices2 = np.array([[0,2,-2],[-3,2,-2],[-3,2,2],[0,2,2],
                      [0,-2,-2],[-3,-2,-2],[-3,-2,2],[0,-2,2]]) - displacement
faces2 = np.array([[1,2,3],[1,4,3],[1,4,8],[1,5,8],[1,2,6],[1,5,6],
                   [7,8,4],[7,3,4],[7,3,2],[7,2,6],[7,8,5],[7,6,5]])

poly1 = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
                       generation_type='full', nb_vertices=4, vertices=vertices1,
                       faces=faces1, radius=10., tol=0., number=None, seed=None,
                    xr=1., yr=1., zr=1.)
poly2 = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
                       generation_type='full', nb_vertices=4, vertices=vertices2,
                       faces=faces2, radius=10., tol=0., number=None, seed=None,
                    xr=1., yr=1., zr=1.)



omega_comb = np.array([0,0,2*np.pi])

coor = np.vstack((poly1.nodes[1].coor.copy(), poly2.nodes[1].coor.copy()))
mass = np.array([poly1.contactors[0].volume * rho, poly2.contactors[0].volume * rho])
CM,coor_cm = Get_CenterMass(coor,mass)
print('CM = ', CM)
vertices3 = Shift2CM(vertices1,CM)#/Lunit*(volume_apophis)**(1/3)
vertices4 = Shift2CM(vertices2,CM)#/Lunit*(volume_apophis)**(1/3)
vel_cm, vel_in = Get_InitialVelocity(coor, mass, omega_comb)

poly3 = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
                       generation_type='full', vertices=vertices3,
                       faces=faces1, radius=10., tol=0., number=None, seed=None,
                    xr=1., yr=1., zr=1.)
poly4 = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
                       generation_type='full', vertices=vertices4,
                       faces=faces2, radius=10., tol=0., number=None, seed=None,
                    xr=1., yr=1., zr=1.)

print(poly3.contactors[0].volume)
print(poly4.contactors[0].volume)

r1 = poly3.nodes[1].coor.copy()
r2 = poly4.nodes[1].coor.copy()
vel1 = np.cross(omega_comb,r1)
vel2 = np.cross(omega_comb,r2)
vel_con = np.vstack((vel1, vel2))

vel = vel_cm
velocity3 = list(np.concatenate((vel[0,:], omega_comb)))
velocity4 = list(np.concatenate((vel[1,:], omega_comb)))
# poly1.imposeDrivenDof(component=[1, 2, 3], dofty='vlocy')
# poly2.imposeDrivenDof(component=[1, 2, 3], dofty='vlocy')
poly3.imposeInitValue(component=[1,2,3,4,5,6], value=velocity3)
poly4.imposeInitValue(component=[1,2,3,4,5,6], value=velocity4)

bodies.addAvatar(poly3)
bodies.addAvatar(poly4)


# MATHIEU TO HAI-SHUO
# Be careful. Here you define a normal and a tangential restitution. If there is a sense for sphere restitution
# is more ambiguous for polyedra.
# You can used a normal one but i recommand to put to zero the tangential one.
LawSPSPx = tact_behav(name='rst01', law='RST_CLB', rstn=1.0, rstt=0.0, fric=1.)

tacts   += LawSPSPx

svSPSPx = see_table(CorpsCandidat='RBDY3', candidat='POLYR',colorCandidat='BLEUx', behav=LawSPSPx,
                     CorpsAntagoniste='RBDY3', antagoniste='POLYR', colorAntagoniste='BLEUx', alert=0.)
svs += svSPSPx

writeBodies(bodies, chemin='DATBOX/')
writeBulkBehav(mat, chemin='DATBOX/', dim=dim , gravy=[0.,0.,0.])
writeTactBehav(tacts, svs, chemin='DATBOX/')
writeDrvDof(bodies, chemin='DATBOX/')
writeDofIni(bodies, chemin='DATBOX/')
writeVlocRlocIni(chemin='DATBOX/')

# try:
#   visuAvatars(bodies)
# except:
#   pass

