# coding=utf-8

"""
@author: Hai-Shuo Wang
@email: wallen1732@gamil.com
@file: DataCombine.py
@date: 1/11/23 14:04
@desc: 
"""
import re
import numpy as np

def Read_connectivity(f_con):
    # read data from Connectivity
    i = 0
    connectivity = []
    line_con = f_con.readline()
    while line_con:
        # extract number from string
        string_numbers = np.array(re.findall(r"\d+", line_con))

        # transfer valid numbers from string to int
        int_numbers = list(map(int, [string_numbers[1], string_numbers[3], string_numbers[5]]))
        # print(int_numbers)

        # record these numbers
        connectivity.append(int_numbers)

        # read a string line
        line_con = f_con.readline()

        i += 1
        # if i > 5: break

    connectivity = np.array(connectivity)
    return connectivity

def Read_vertices(f_apo):

    # read data from apophis_v233s7_vert2_new
    i = 0
    vertices = []
    faces = []
    line_apo = f_apo.readline()
    while line_apo:

        # extract number from string
        string_line = np.array(re.findall(r"\d+\.?\d*", line_apo))

        # transfer these number from string to float
        double_line = list(map(float, string_line))

        # record these numbers
        if line_apo[0]== 'v':
            vertices.append(double_line)
            # print(double_line)
        if line_apo[0]== 'f':
            faces.append(double_line)

        # read a string line
        line_apo = f_apo.readline()

        i += 1
        # if i > 5: break
    vertices = np.array(vertices)
    faces = np.array(faces)
    return vertices,faces


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    f_con = open("Connectivity.txt","r")
    connectivity = Read_connectivity(f_con)
    f_con.close()
    print(connectivity)

    f_apo = open("apophis_v233s7_vert2_new.mod.wf", "r")
    vertices,faces = Read_vertices(f_apo)
    f_apo.close()
    print(vertices)
    print(faces)


