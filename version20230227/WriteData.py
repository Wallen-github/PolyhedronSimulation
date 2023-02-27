# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: ReadWriteData.py
@date: 1/10/23 19:31
@desc: 
"""
import numpy as np
import os
import ReadData


def Write_bdyty(fid):
    bodyatype = 'RBDY3'
    bodynumber = 1
    fid.write('$bdyty\n')
    fid.write(' %5s  %5d\n' % (bodyatype, bodynumber))
    fid.write('$blmty\n')


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    chemin = 'DATBOX/'
    fid = open(os.path.join(chemin, 'BODIES.DAT'), 'w+')
    # w:    creat and write
    # w+:   creat, write and read
    # r:    read, error in no file
    # r+:   write and read, error in no file
    # a:    append write
    # a+:   append write and read
    Write_bdyty(fid)

    #%%
    # on incremente le numero du prochain rigide
    bulk_number = 1
    bulkmaterialnom = 'STONE'
    # bulkavrd = 7.8159264E+00
    # inertia = np.array([6.8426699E+03, 6.8426699E+03, 6.8426699E+03])
    bulkavrd = 1.0742178E-01
    inertia = np.array([1.2083981E-05,1.6518386E-05,2.5501754E-05]) # Apophis
    # to generalize
    line = ' %5s  %5d  behav  %5s' % ('PLAIN', bulk_number, bulkmaterialnom)
    line += ' %5s=%14.7E' % ('avrd', bulkavrd)
    line += '\n'
    line += '                           '
    line += ' %5s=%14.7E' % ('I1  ', inertia[0])
    line += ' %5s=%14.7E' % ('I2  ', inertia[1])
    line += ' %5s=%14.7E' % ('I3  ', inertia[2])
    line += '\n'
    fid.write(line)

    #%%
    # on initialise l'ecriture des noeuds
    fid.write('$nodty\n')
    ntype = 'NO6xx'
    coorsize = 6
    coor = np.array([60,0,0,0,0,0])
    # on ecrit les noeuds comme precedemment (cf. writeElementNodesInBodies)
    ligne = ' %5s  %5d                ' % (ntype, 1)
    for j in range(coorsize):
        if j % 3 == 2 or j == coorsize - 1:
            ligne += 'coo%d=%14.7E  \n' % (j + 1, coor[j])
        else:
            ligne += 'coo%d=%14.7E  ' % (j + 1, coor[j])
        if j % 3 == 2 and not j == coorsize - 1:
            ligne += '                             '
    fid.write(ligne)

    #%%
    f_apo = open("apophis_v233s7_vert2_new.mod.wf", "r")
    vertices, faces = ReadData.Read_vertices(f_apo)
    f_apo.close()
    shape = 'POLYR'
    number=1
    color = 'BLEUx'
    nb_vertices = 4
    nb_faces = 4
    vertices = np.array([[0,0,10],
                        [9.4280904E+00,0,-3.3333333E+00],
                        [-4.7140452E+00, 8.1649658E+00, -3.3333333E+00],
                        [-4.7140452E+00, -8.1649658E+00, -3.3333333E+00]])
    connectivity = np.array([[4,2,3],
                             [1,2,3],
                             [1,4,3],
                             [1,4,2]])

    # f_apo = open("apophis_v233s7_vert2_new.mod.wf", "r")
    # vertices, faces = ReadData.Read_vertices(f_apo)
    # connectivity = faces
    # f_apo.close()
    # nb_vertices = vertices.shape[0]
    # nb_faces = faces.shape[0]

    # nb_vertices = 2000
    # nb_faces = 3996
    # f_con = open("Connectivity.txt", "r")
    # import ReadData
    # connectivity = ReadData.Read_connectivity(f_con)
    # f_con.close()
    # f_apo = open("apophis_v233s7_vert2_new.mod.wf", "r")
    # vertices, faces = ReadData.Read_vertices(f_apo)
    # f_apo.close()

    # tacts writing
    fid.write('$tacty                                                                  \n')
    # on ecrit la premiere ligne decrivant le contacteur
    line = ' %5s  %5d  color  %5s  nb_vertex=%7d    nb_faces=%7d\n' % (shape, number, color, \
                                                                       nb_vertices, nb_faces)
    # on ajoute les lignes donnant les coordonnees des sommets
    for i in range(nb_vertices):
        line += '                             '
        line += 'coo1=%14.7E  coo2=%14.7E  coo3=%14.7E\n' % (
        vertices[i, 0], vertices[i, 1], vertices[i, 2])
    # on ajoute les lignes donnant les connectivites des faces
    for i in range(nb_faces):
        line += '                             '
        line += 'ver1=%7d         ver2=%7d         ver3=%7d\n' % \
                (connectivity[i, 0], connectivity[i, 1], connectivity[i, 2])
    fid.write(line)
    fid.write('$$$$$$\n')
    fid.close()
