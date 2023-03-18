# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: Computation.py
@date: 1/10/23 13:50
@desc: 
"""

from __future__ import print_function
import os, sys

from matplotlib import pyplot as plt
from pylmgc90 import chipy
from pykdgrav import *
import EarthGravity

sys.path.append('..')

chipy.Initialize()

chipy.checkDirectories()

chipy.utilities_DisableLogMes()
# timer gravity
timer_id = chipy.timer_GetNewTimer('gravity computation')

chipy.SetDimension(3)

# dt = 5.e-5
dt = 1. # sec
theta = 0.5
nb_steps = 3*24*3600 # sec

echo = 0

freq_display = 400
ref_radius = 5.
freq_write = 400

# freq_update = 600

chipy.nlgs_3D_DiagonalResolution()
chipy.RBDY3_NewRotationScheme()

chipy.PRPRx_UseCpCundallDetection(300)
chipy.PRPRx_LowSizeArrayPolyr(20)

# type = 'Stored_Delassus_Loops         '
stype = 'Exchange_Local_Global         '
norm = 'QM/16'
tol = 0.1e-3
relax = 1.0
gs_it1 = 10
gs_it2 = 200

chipy.utilities_logMes('INIT TIME STEPPING')
chipy.TimeEvolution_SetTimeStep(dt)
chipy.Integrator_InitTheta(theta)

chipy.utilities_logMes('READ BODIES')
chipy.ReadBodies()
chipy.LoadTactors()

chipy.utilities_logMes('READ INI DOF')
chipy.ReadIniDof()

chipy.utilities_logMes('READ BEHAVIOURS')
chipy.ReadBehaviours()
chipy.LoadBehaviours()

chipy.utilities_logMes('READ INI Vloc Rloc')
chipy.ReadIniVlocRloc()

chipy.utilities_logMes('READ DRIVEN DOF')
chipy.ReadDrivenDof()

chipy.WriteBodies()
#
chipy.WriteBehaviours()

chipy.WriteDrivenDof()

chipy.ComputeMass()

# chipy.RBDY3_FatalDamping()

nbR3 = chipy.RBDY3_GetNbRBDY3()

mass = np.zeros(nbR3)
for i in range(nbR3):
    mass[i] = chipy.RBDY3_GetMass(i + 1)

coor = np.empty([nbR3, 6], dtype=float)
p_coor = np.empty([nbR3, 3], dtype=float)
vbeg = np.empty([nbR3, 6], dtype=float)
fext = np.empty([nbR3, 6], dtype=float)
cij = np.zeros(3)
PosVecE = np.empty([nb_steps+1, 6], dtype=float)
Record_posvel = np.zeros([nb_steps+1, 6], dtype=float)

chipy.OpenDisplayFiles()
chipy.WriteDisplayFiles(1)

# Compute the Earth's initial state
muA = 2.650 # Gravitational parameter of asteroid Assumes density of 2 g∕cm^3
GG = 6.6742e-11 # G \sim N*m^2/kg^2
MA = muA/GG
PosVecE0, Unit = EarthGravity.InitialEarthPV(MA,GG=GG)
PosVecE[0,:] = PosVecE0.copy()

for k in range(1, nb_steps + 1):
    print(k, '/', (nb_steps + 1))
    #
    chipy.IncrementStep()

    chipy.ComputeFext()

    fext = np.empty([nbR3, 6], dtype=float)

    for i in range(0, nbR3, 1):
        coor[i, :] = chipy.RBDY3_GetBodyVector('Coorb', i + 1)
        p_coor[i, 0:3] = coor[i, 0:3]
        fext[i, :] = 0.

    chipy.timer_StartTimer(timer_id)
    fext[:, 0:3] = Accel(p_coor, mass, G=GG)
    # Compute the Earth position and velocity
    timespan = [(k-1)*dt,k*dt]
    Gravorder = 2
    PosVecE[k,:], PosVecSol = EarthGravity.EarthPos(timespan, PosVecE[k-1,:], Unit, Gravorder)

    # Record the position and velocity
    for i in range(0, nbR3, 1):
        Record_posvel[k, :] = Record_posvel[k, :] + mass[i] * coor[i, :]
    Record_posvel[k, :] = Record_posvel[k, :] / np.sum(mass)

    fextCM = EarthGravity.CenterMassAccel(GG, np.sum(mass), PosVecE[i, 0:3] + Record_posvel[k, 0:3],
                                          Record_posvel[k, 0:3])

    for i in range(0, nbR3, 1):
        fextE = EarthGravity.EarthAccel(GG,mass[i],p_coor[i, 0:3],PosVecE[i, 0:3]+Record_posvel[k, 0:3])
        fext[i, :] = fext[i, :] * mass[i] + fextE #+ fextCM
    chipy.timer_StopTimer(timer_id)

    for i in range(0, nbR3, 1):
        chipy.RBDY3_PutBodyVector('Fext_', i + 1, fext[i, :])

    chipy.ComputeBulk()
    chipy.ComputeFreeVelocity()

    chipy.SelectProxTactors()
    chipy.RecupRloc()

    chipy.ExSolver(stype, norm, tol, relax, gs_it1, gs_it2)

    chipy.StockRloc()

    chipy.ComputeDof()
    chipy.UpdateStep()

    chipy.WriteDisplayFiles(freq_display, ref_radius)

    chipy.WriteOutDof(freq_write)
    chipy.WriteOutVlocRloc(freq_write)

    chipy.overall_CleanWriteOutFlags()

chipy.WriteLastDof()
chipy.CloseDisplayFiles()

chipy.Finalize()

#定义坐标轴
fig = plt.figure()
ax1 = plt.axes(projection='3d')
ax1.plot3D(PosVecE[:,0],PosVecE[:,1],PosVecE[:,2],label='geocentric hyperbolic orbit',color='black')
# ax1.plot3D(Record_posvel[:,0],Record_posvel[1,:],Record_posvel[2,:],label='real orbit',color='blue')
ax1.scatter3D(PosVecE[0,0],PosVecE[0,1],PosVecE[0,2],label='initial position',c='red',s=20)
ax1.scatter3D(0,0,0, c='blue',s=20,label='Earth')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel('z')
plt.legend()
plt.show()

