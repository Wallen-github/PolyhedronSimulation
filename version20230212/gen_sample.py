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

rho = 2000

Omega = 3e-5

stone = material(name='STONE', materialType='RIGID', density=rho)
mat.addMaterial(stone)
mod = model(name='rigid', physics='MECAx', element='Rxx3D', dimension=dim)

# center1 = np.array([0,0,0])
# center2 = np.array([0,5,0])

vertices1 = np.array([[0,50,10],[10,50,0],[0,60,0],[0,50,0]])
faces1 = np.array([[1,2,3],[1,2,4],[1,3,4],[2,3,4]])
vertices2 = np.array([[0,0,10],[10,10,0],[-10,10,0],[0,-10,0],[0,0,0]]) # Convex
# vertices2 = np.array([[0,0,10],[10,10,0],[-10,10,0],[0,-10,0],[0,15,0]]) # Concave
faces2 = np.array([[1,2,5],[1,2,4],[1,3,5],[1,3,4],[2,4,5],[3,4,5]])
poly1 = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
                       generation_type='full', nb_vertices=4, vertices=vertices1,
                       faces=faces1, tol=0., number=None, seed=None)
poly2 = rigidPolyhedron(model=mod, material=stone, color='BLEUx',
                       generation_type='full', nb_vertices=4, vertices=vertices2,
                       faces=faces2, tol=0., number=None, seed=None,)

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

