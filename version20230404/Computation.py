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

import numpy as np
import pandas as pd
from pylmgc90 import chipy
from pykdgrav import *
from ComboProperties import *
from FlybyOrbit import *
from ShapeModel import SetUnit
from version_test.Shift2CM import Shift2CM

sys.path.append('..')

chipy.Initialize()

chipy.checkDirectories()

chipy.utilities_DisableLogMes()
# timer gravity
timer_id = chipy.timer_GetNewTimer('gravity computation')

chipy.SetDimension(3)


GG = 1
Lunit, Munit, Tunit = SetUnit(GG=6.6742e-11,volume_cos=1.986E7,rho=2000)
# time = int(30.6*3600/Tunit) # one spin period
time = int(3 * 24 * 3600 / Tunit) # 3 days

dt = 1E-3
theta = 0.5
nb_steps = 1000 * time

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
velo = np.empty([nbR3, 6], dtype=float)
p_coor = np.empty([nbR3, 6], dtype=float)
vbeg = np.empty([nbR3, 6], dtype=float)
fext = np.empty([nbR3, 6], dtype=float)
gravpot = np.empty([nbR3, 1], dtype=float)

record_pvw = np.empty([nb_steps+1,nbR3*9],dtype=float)
record_pvw_cm = np.empty([nb_steps+1,nbR3*9],dtype=float)
record_momentum = np.empty([nb_steps+1,4],dtype=float)
record_energy = np.zeros([nb_steps+1,2],dtype=float)
record_effquanti = np.empty([nb_steps+1,3],dtype=float)

chipy.OpenDisplayFiles()
chipy.WriteDisplayFiles(1)

PosVecFly, UnitFlyby = InitialEarthPV(Munit, GG=6.6742e-11)

for k in range(1, nb_steps + 1):
    print(k, '/', (nb_steps + 1))
    #
    chipy.IncrementStep()

    chipy.ComputeFext()

    for i in range(0, nbR3, 1):
        coor[i, :] = chipy.RBDY3_GetBodyVector('Coorb', i + 1)
        velo[i, :] = chipy.RBDY3_GetBodyVector('Vfree', i + 1)
        fext[i, :] = 0.

    timespan = [(k - 1) * dt * Tunit, k * dt * Tunit]  # sec
    PosVecFly, Sol = EarthPos(timespan, PosVecFly, UnitFlyby, Gravorder=1)

    # record data
    CMP, coor_cm = Get_CenterMass(coor[:,0:3], mass)
    CMV, velo_cm = Get_CenterMass(velo[:,0:3], mass)
    for i in range(0, nbR3, 1):
        record_pvw[k - 1, i * 9:i * 9 + 3] = coor[i, 0:3]*Lunit + PosVecFly[0:3]
        record_pvw[k - 1, i * 9 + 3:i * 9 + 6] = velo[i,0:3] * Lunit/Tunit + PosVecFly[3:6]
        record_pvw[k - 1, i * 9 + 6:(i + 1) * 9] = velo[i,3:6] * 3600/Tunit
        record_pvw_cm[k - 1, i * 9:i * 9 + 3] = coor_cm[i, 0:3] * Lunit
        record_pvw_cm[k - 1, i * 9 + 3:i * 9 + 6] = velo_cm[i, 0:3] * Lunit / Tunit
        record_pvw_cm[k - 1, i * 9 + 6:(i + 1) * 9] = velo[i, 3:6] * 3600 / Tunit

    if k>=2:
        record_energy[k-1,0], record_energy[k-1,1], record_momentum[k-1,:] = Get_EnergyMomentum(nbR3,GG=GG)
        total_inertia = Get_TotalInertia(nbR3,mass,coor)
        omega_e, inertia_d, Inertia_t = Get_EffectiveQuantity(record_energy[k-1,1], record_momentum[k-1,3],total_inertia)
        record_effquanti[k - 1, :] = np.array([omega_e * 180/np.pi * 3600 / Tunit, inertia_d, Inertia_t])


    r0i = np.zeros([6])
    qcm = np.zeros([6])
    qcm[0:3] = PosVecFly[0:3] / Lunit
    chipy.timer_StartTimer(timer_id)
    # fext[:, 0:3] = Accel(p_coor, mass, G=6.6742e-11)
    fext[:, 0:3] = Accel(coor[:, 0:3], mass, G=GG)
    for i in range(0, nbR3, 1):
        mE = 5.972E24 / Munit
        r0i[0:3] = coor[i, 0:3] + qcm[0:3]
        # rE = np.linalg.norm(coor[i, 0:3])
        fext[i, :] = fext[i, :] * mass[i] - GG*mE*mass[i]/np.linalg.norm(r0i[0:3])**3 * r0i + GG*mE*mass[i]/np.linalg.norm(qcm[0:3])**3 * qcm #- GG*mE*mass[i]/rE**3 * coor[i, 0:3]

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

fid = open('./DATA_Record/pos_vel_spin_cm.dat', 'w+')
for k in range(1,len(record_pvw_cm)-1,10):
    line = ' '
    for j in range(len(record_pvw_cm[k,:])):
        line += '%14.7E \t' % (record_pvw_cm[k,j])
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

fid = open('./DATA_Record/effectquantity.dat', 'w+')
for k in range(1,len(record_effquanti)-1,10):
    line = ' '
    for j in range(len(record_effquanti[k,:])):
        line += '%14.7E \t' % (record_effquanti[k,j])
    line += '\n'
    fid.write(line)
fid.close()







