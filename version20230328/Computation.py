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
from ComboProperties import *
from version20230328.FlybyOrbit import *
from version20230328.ShapeModel import SetUnit
from version_test.Shift2CM import Shift2CM

sys.path.append('..')

chipy.Initialize()

chipy.checkDirectories()

chipy.utilities_DisableLogMes()
# timer gravity
timer_id = chipy.timer_GetNewTimer('gravity computation')

chipy.SetDimension(3)

Lunit, Munit, Tunit = SetUnit()
# time = (30.6*60*60)/Tunit
time = 3*24*3600/Tunit

dt = 1E-3
theta = 0.5
nb_steps = 1000*int(time)

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
record_pvw_post = np.empty([nb_steps+1,nbR3*9],dtype=float)
record_momentum = np.empty([nb_steps+1,4],dtype=float)
record_energy = np.zeros([nb_steps+1,2],dtype=float)
record_effquanti = np.empty([nb_steps+1,3],dtype=float)
record_flyby = np.empty([nb_steps+1,7],dtype=float)

GG = 1#6.6742e-11
record_flyby[0, 0] = 0.
PosVecE0, Unit = InitialEarthPV(MA=Munit, GG=6.6742e-11) # output with unit
record_flyby[0, 1:7] = PosVecE0.copy()

chipy.OpenDisplayFiles()
chipy.WriteDisplayFiles(1)

for k in range(1, nb_steps + 1):
    print(k, '/', (nb_steps + 1))
    #
    chipy.IncrementStep()

    chipy.ComputeFext()

    for i in range(0, nbR3, 1):
        coor[i, :] = chipy.RBDY3_GetBodyVector('Coorb', i + 1)
        velo[i, :] = chipy.RBDY3_GetBodyVector('Vfree', i + 1)
        fext[i, :] = 0.

    # record data
    CM, coor_cm = Get_CenterMass(coor[:,0:3], mass)
    CM, velo_cm = Get_CenterMass(velo[:,0:3], mass)
    for i in range(0, nbR3, 1):
        record_pvw[k - 1, i * 9:i * 9 + 3] = coor_cm[i, 0:3]*Lunit
        record_pvw[k - 1, i * 9 + 3:i * 9 + 6] = velo_cm[i,0:3] * Lunit/Tunit
        record_pvw[k - 1, i * 9 + 6:(i + 1) * 9] = velo[i,3:6] * 3600/Tunit

    chipy.timer_StartTimer(timer_id)
    # fext[:, 0:3] = Accel(p_coor, mass, G=6.6742e-11)
    fext[:, 0:3] = Accel(coor[:,0:3], mass, G=GG)

    timespan = [(k-1)*dt * Tunit, k*dt * Tunit]  # sec
    record_flyby[k, 0] = float(k * dt)
    record_flyby[k, 1:7], Sol = EarthPos(timespan, record_flyby[k - 1, 1:7], Unit, Gravorder=2)
    rE = record_flyby[k, 1:7] / Lunit
    mE = 5.972E24 / Munit

    for i in range(0, nbR3, 1):
        fext[i, :] = fext[i, :] * mass[i] - mass[i] * mE / (np.linalg.norm(rE))**3 * rE

    chipy.timer_StopTimer(timer_id)

    if k>=2:
        record_energy[k-1,0], record_energy[k-1,1], record_momentum[k-1,:] = Get_EnergyMomentum(nbR3,GG=GG)
        total_inertia = Get_TotalInertia(nbR3)
        record_effquanti[k-1,:] = Get_EffectiveQuantity(record_energy[k-1,1], record_momentum[k-1,3],total_inertia)

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

fid = open('./DATA_Record/effectquantity.dat', 'w+')
for k in range(1,len(record_effquanti)-1,10):
    line = ' '
    for j in range(len(record_effquanti[k,:])):
        line += '%14.7E \t' % (record_effquanti[k,j])
    line += '\n'
    fid.write(line)
fid.close()

fid = open('./DATA_Record/flybyOrbit.dat', 'w+')
for k in range(1,len(record_flyby)-1,10):
    line = ' '
    for j in range(len(record_flyby[k,:])):
        line += '%14.7E \t' % (record_flyby[k,j])
    line += '\n'
    fid.write(line)
fid.close()






