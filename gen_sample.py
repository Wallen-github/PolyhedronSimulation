# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: gen_sample.py
@date: 1/10/23 13:38
@desc: 
"""
import numpy as np
from pylmgc90.pre import *

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

rho = 5164

Omega = 3e-5

stone = material(name='STONE', materialType='RIGID', density=rho)
mat.addMaterial(stone)
mod = model(name='rigid', physics='MECAx', element='Rxx3D', dimension=dim)

center1 = np.array([30,0,0])
center2 = np.array([-30,0,0])

poly1 = rigidPolyhedron(model=mod, material=stone, center=center1, color='BLEUx',
                       generation_type='random', nb_vertices=10, vertices=None,
                       faces=None, radius=10., tol=0., number=None, seed=None,
                    xr=1., yr=1., zr=1.)
poly2 = rigidPolyhedron(model=mod, material=stone, center=center2, color='BLEUx',
                       generation_type='random', nb_vertices=15, vertices=None,
                       faces=None, radius=10., tol=0., number=None, seed=None,
                    xr=1., yr=1., zr=1.)

bodies.addAvatar(poly1)
bodies.addAvatar(poly2)


# MATHIEU TO HAI-SHUO
# Be careful. Here you define a normal and a tangential restitution. If there is a sense for sphere restitution
# is more ambiguous for polyedra.
# You can used a normal one but i recommand to put to zero the tangential one.
LawSPSPx = tact_behav(name='rst01', law='RST_CLB', rstn=1.0, rstt=0.0, fric=0.)

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

try:
  visuAvatars(bodies)
except:
  pass

