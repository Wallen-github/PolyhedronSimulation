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

import pandas as pd
from pylmgc90 import chipy
from pykdgrav import *
from version_test.Shift2CM import Shift2CM

sys.path.append('..')

chipy.Initialize()

chipy.checkDirectories()

chipy.utilities_DisableLogMes()
# timer gravity
timer_id = chipy.timer_GetNewTimer('gravity computation')

chipy.SetDimension(3)

dt = 1E-4
theta = 0.5
nb_steps = 10000

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
tol = 0.1e-8
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

inertia = chipy.RBDY3_GetGlobInertia(1)

coor = np.empty([nbR3, 6], dtype=float)
coor_post = np.empty([nbR3, 6], dtype=float)
p_coor = np.empty([nbR3, 6], dtype=float)
vbeg = np.empty([nbR3, 6], dtype=float)
fext = np.empty([nbR3, 6], dtype=float)
gravpot = np.empty([nbR3, 1], dtype=float)

record_pvw = np.empty([nb_steps+1,nbR3*9],dtype=float)
record_pvw_post = np.empty([nb_steps+1,nbR3*9],dtype=float)
record_momentum = np.empty([nb_steps+1,(nbR3+1)*3],dtype=float)
record_energy = np.zeros([nb_steps+1,nbR3+1],dtype=float)
record_effquanti = np.empty([nb_steps+1,2],dtype=float)

chipy.OpenDisplayFiles()
chipy.WriteDisplayFiles(1)

for k in range(1, nb_steps + 1):
    print(k, '/', (nb_steps + 1))
    #
    chipy.IncrementStep()

    chipy.ComputeFext()

    fext = np.empty([nbR3, 6], dtype=float)

    for i in range(0, nbR3, 1):
        coor[i, :] = chipy.RBDY3_GetBodyVector('Coorb', i + 1)
        fext[i, :] = 0.

    # convert the position and velocity to the CM frame
    p_coor = Shift2CM(coor, mass)

    # record data
    gravpot = Potential(coor[:, 0:3], mass, G=1)
    for i in range(0, nbR3, 1):
        record_pvw[k - 1, i * 9:i * 9 + 3] = coor[i, 0:3]
        record_pvw[k - 1, i * 9 + 3:(i + 1) * 9] = chipy.RBDY3_GetBodyVector('Vfree', i + 1)
        inertia = chipy.RBDY3_GetGlobInertia(i+1)
        omega_i = record_pvw[k - 1, i * 9 + 6:(i + 1) * 9]
        record_momentum[k - 1,i*3:i*3+3] = mass[i]*np.cross(coor[i, 0:3],coor[i, 3:6]) + np.dot(inertia,omega_i)
        record_momentum[k - 1, nbR3 * 3:nbR3 * 3 + 3] += record_momentum[k - 1,i*3:i*3+3]
        record_energy[k-1, i] = 0.5*mass[i]*np.dot(coor[i, 3:6],coor[i, 3:6]) + 0.5 * np.dot(omega_i,np.dot(inertia,omega_i)) #- gravpot[i]
        record_energy[k - 1, nbR3] += record_energy[k-1, i]
    print(record_energy[k - 1, :])

    # effective spin rate
    # momentum = np.linalg.norm(record_momentum[k - 1, nbR3 * 3:nbR3 * 3 + 3])
    # record_effquanti[k-1,0] = 2* record_energy[k - 1, nbR3] / momentum
    # record_effquanti[k - 1, 1] = momentum**2 / 2* record_energy[k - 1, nbR3]


    chipy.timer_StartTimer(timer_id)
    # fext[:, 0:3] = Accel(p_coor, mass, G=6.6742e-11)
    fext[:, 0:3] = Accel(coor[:,0:3], mass, G=1)
    for i in range(0, nbR3, 1):
        fext[i, :] = fext[i, :] * mass[i] #np.array([0.,2E6,0.,0.,0.,0.])
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




# Record and output data
fid = open('./DATA_Record/pos_vel_spin.dat', 'w+')
for k in range(1,len(record_pvw)-1,10):
    line = ' '
    for j in range(len(record_pvw[k,:])):
        line += '%14.7E \t' % (record_pvw[k,j])
    line += '\n'
    fid.write(line)
fid.close()

fid = open('./DATA_Record/angular_momentum.dat', 'w+')
for k in range(1,len(record_momentum)-1,10):
    line = ' '
    for j in range(len(record_momentum[k,:])):
        line += '%14.7E \t' % (record_momentum[k,j])
    line += '\n'
    fid.write(line)
fid.close()

fid = open('./DATA_Record/energy.dat', 'w+')
for k in range(1,len(record_energy)-1,10):
    line = ' '
    for j in range(len(record_energy[k,:])):
        line += '%14.7E \t' % (record_energy[k,j])
    line += '\n'
    fid.write(line)
fid.close()

fid = open('./DATA_Record/effctive_quatity.dat', 'w+')
for k in range(1,len(record_effquanti)-1,10):
    line = ' '
    for j in range(len(record_effquanti[k,:])):
        line += '%14.7E \t' % (record_effquanti[k,j])
    line += '\n'
    fid.write(line)
fid.close()






