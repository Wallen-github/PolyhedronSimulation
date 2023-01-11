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


def Write_bdyty(bodyatype, bodynumber, chemin, filename):
    fid = open(os.path.join(chemin, filename), 'a+')
    fid.write('$bdyty\n')
    fid.write(' %5s  %5d\n' % (bodyatype, bodynumber))
    fid.write('$blmty\n')


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    chemin = 'DATBOX/'
    bodyatype = 'RBDY3'
    bodynumber = 1
    fid = open(os.path.join(chemin, 'BODIES.DAT'), 'a+')
    fid.write('$bdyty\n')
    fid.write(' %5s  %5d\n' % (bodyatype, bodynumber))
    fid.write('$blmty\n')

    #%%
    # on incremente le numero du prochain rigide
    bulk_number = bulk_number + 1
    bulkmaterialnom = 'STONE'
    bulkavrd = 7.8159264E+00
    # to generalize
    line = ' %5s  %5d  behav  %5s' % ('PLAIN', bulk_number, bulkmaterialnom)
    line += ' %5s=%14.7E' % ('avrd', bulkavrd)
    line += '                           '
    line += ' %5s=%14.7E' % ('I1  ', bulk.inertia[0])
    line += ' %5s=%14.7E' % ('I2  ', bulk.inertia[1])
    line += ' %5s=%14.7E' % ('I3  ', bulk.inertia[2])
    line += '\n'
    fid.write(line)

    #%%
    # on initialise l'ecriture des noeuds
    fid.write('$nodty\n')
    ntype = 'NO6xx'
    coorsize = 6
    coor = np.array([30,0,0,0,0,0])
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
    shape = 'POLYR'
    number=1
    color = 'BLEUx'
    nb_vertices = 2000
    nb_faces = 3996
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
