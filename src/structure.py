# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 14:30:21 2025

@author: User
"""

import numpy as np
import matplotlib.pyplot as plt
import time
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

class LatticeLayer:
    def __init__(self, translational_vectors, sublattice, rotation_angle):
        self.translational_vectors = translational_vectors
        self.sublattice = sublattice
        self.rotation_angle = rotation_angle
        self.lattice_points = None
        

    def generate_lattice_points(self, x_range, y_range):
        size = x_range.shape[0] * y_range.shape[0]
        
        # Initialize the lattice point array
        lattice_points = np.zeros((size*2, 2), dtype=np.float64)
                
        # Use meshgrid to generate a grid coordinate matrix
        X, Y = np.meshgrid(x_range, y_range)

        # Compute unrotated lattice point coordinates
        lattice_points[:, 0] = (self.sublattice[:, 0, np.newaxis] + X.flatten()).flatten()
        lattice_points[:, 1] = (self.sublattice[:, 1, np.newaxis] + Y.flatten()).flatten()

        # Save calculation results
        self.lattice_points = lattice_points


    def apply_rotation(self):
        # rotation matrix 
        rot_mat = np.array([[np.cos(self.rotation_angle/2), -np.sin(self.rotation_angle/2)], 
                            [np.sin(self.rotation_angle/2),  np.cos(self.rotation_angle/2)]])
        self.lattice_points = np.dot(self.lattice_points, rot_mat)
        return self.lattice_points


class TwistedLayer:
    def __init__(self, d, x_range, y_range, name='rot_layer1_sub_a.npy'):
        self.layers = None
        self.c = 3**(1/2)
        self.d = d
        self.translational_vectors = [np.array([self.c*d, 0]), np.array([0, 3*d])] #translational vector of unitcell
        self.npy_file_name = name
        self.x_range = x_range
        self.y_range = y_range

    
    def add_layer(self,  stack_conf, rotation_angle, round_digit = 12, 
                         select_sublattice=None, save_points=False):
        
        # Determine the sublattice structure of each layer based on "stack_conf"
        # for the base layer
        if stack_conf == 'Base':  
            if select_sublattice == 'A':
                layer_sub_a = np.array([[0,      0], [self.c*self.d/2, 3*self.d/2]]) #unit cell base layer sublattice A
            elif select_sublattice == 'B':
                layer_sub_b = np.array([[0, self.d], [self.c*self.d/2, 5*self.d/2]]) #unit cell base layer sublattice B   
        
        # stacking configuration correspond to the base layer
        elif stack_conf == 'AA':  
            if select_sublattice == 'A':
                layer_sub_a = np.array([[0,      0], [self.c*self.d/2, 3*self.d/2]]) #unit cell AA sublattice A
            if select_sublattice == 'B':
                layer_sub_b = np.array([[0, self.d], [self.c*self.d/2, 5*self.d/2]]) #unit cell AA sublattice B
        
        # stacking configuration correspond to the base layer
        elif stack_conf == 'AB':  
            if select_sublattice == 'A':
                layer_sub_a = np.array([[0, -self.d], [self.c*self.d/2, 3*self.d/2-self.d]]) #unit cell AB sublattice A
            if select_sublattice == 'B':
                layer_sub_b = np.array([[0,       0], [self.c*self.d/2, 5*self.d/2-self.d]]) #unit cell AB sublattice B
        
        
        # Create a LatticeLayer object based on the selected sublattice
        if select_sublattice == 'A':
            new_layer = LatticeLayer(self.translational_vectors, layer_sub_a, rotation_angle)
        elif select_sublattice == 'B':
            new_layer = LatticeLayer(self.translational_vectors, layer_sub_b, rotation_angle)
           
        new_layer.generate_lattice_points(self.x_range, self.y_range)
        self.layers = new_layer
        
        points = np.around(self.layers.apply_rotation(),round_digit) # round it!
        
        if save_points == True:
            np.save(self.npy_file_name, points)
        return points
    
    def find_coincident_atoms(self, layer1, layer2):
        nrows, ncols = layer1.shape
        dtype = {'names': ['f{}'.format(i) for i in range(ncols)],
                 'formats': ncols * [layer1.dtype]}

        # Find coincident atoms
        coincident_atoms = np.intersect1d(layer1.view(dtype), layer2.view(dtype))
        coincident_atoms = coincident_atoms.view(layer1.dtype).reshape(-1, ncols)
        return coincident_atoms
    

class SupercellFinder:
    def __init__(self, coincident_atoms):
        self.coincident_atoms = coincident_atoms
        self.vertex = np.array([])
        self.atomO = None
        self.atomA = None
        self.atomB = None
        self.atomC = None
    
    def find_nearest_vector(self, value):
        idx = np.array([np.linalg.norm(x + y) for (x, y) in self.coincident_atoms - value]).argmin()
        return self.coincident_atoms[idx]

    def find_vertices(self):
        # Find the (0,0) or the closest point as the starting point of the supercell
        pt = np.array([0, 0])
        self.vertex = np.append(self.vertex, self.find_nearest_vector(pt))  # find closest O(0,0) point
        print('vertex =', self.vertex)
        self.atomO = self.vertex
        
        # Find candidate vertex atoms
        left_O_max = self.coincident_atoms[np.where(self.coincident_atoms[:, 0] < self.vertex[0])][:, 0].max()
        right_O_min = self.coincident_atoms[np.where(self.coincident_atoms[:, 0] > self.vertex[0])][:, 0].min()
        
        test1 = self.coincident_atoms[np.where(self.coincident_atoms[:, 0] < self.vertex[0])]  # for atom C
        test2 = self.coincident_atoms[np.where(self.coincident_atoms[:, 0] > self.vertex[0])]  # for atom B
        
        test_array1 = test1[np.where(test1[:, 0] == left_O_max)]
        test_array2 = test2[np.where(test2[:, 0] == right_O_min)]
        
        test_y1 = test_array1[np.where(test_array1[:, 1] > self.vertex[1])][:, 1].min()
        test_y2 = test_array2[np.where(test_array2[:, 1] > self.vertex[1])][:, 1].min()
        
        # Define atoms
        Bx = test_array2[np.where(test_array2[:, 1] > self.vertex[1])][:, 0].min()
        By = test_array2[np.where(test_array2[:, 1] > self.vertex[1])][:, 1].min()
        self.atomB = np.array([Bx, By])
        
        test3 = self.coincident_atoms[np.where(self.coincident_atoms[:, 0] == self.vertex[0])]  # for atom A'
        Apy = test3[np.where(test3[:, 1] > self.vertex[1])][:, 1].min()
        Apx = self.vertex[0]
        self.atomA = np.array([Apx, Apy])
        
        Cx = test_array1[np.where(test_array1[:, 1] > self.vertex[1])][:, 0].max()
        Cy = test_array1[np.where(test_array1[:, 1] > self.vertex[1])][:, 1].min()
        self.atomC = np.array([Cx, Cy])
        
        from scipy.spatial import distance
        if By == Cy:
            #pseudo code:
            #if L1 == L2 => this is still tilt diamond 
            #if not => this is hexagon(vertical hexagon) 
            if round(distance.euclidean(self.atomA, self.atomO),10) == round(distance.euclidean(self.atomB, self.atomO),10):
                #the supercell is tilt diamond shape as before
                print("supercell is 'horizontal diamond' ")
            else:
                if round(distance.euclidean(self.atomA, self.atomB),10) == round(distance.euclidean(self.atomO, self.atomB),10):
                    print("supercell is 'vertical diamond' ")
                else:
                    print("supercell is 'vertical hexagon' ")
                    if By > Apy:
                        #type 1 vertical hexagon A' is lower than B
                        Ay = test3[np.where(test3[:,1]>Apy)][:,1].min()
                        Ax = Apx
                        self.atomA = [Ax,Ay]
                    else:
                        #type 2 vertical hexagon A' is higher than B
                        Bpx = test2[np.where(test2[:,0]==Bx)][:,0].min()
                        Bpy = test2[np.where(test2[:,1]>By)][:,1].min()
                        self.atomB = [Bpx,Bpy]
                        Cpx = test1[np.where(test1[:,0]==Cx)][:,0].min()
                        Cpy = test1[np.where(test1[:,1]>Cy)][:,1].min()
                        self.atomC = [Cpx,Cpy]
                        Ay = test3[np.where(test3[:,1]>Apy)][:,1].min()
                        Ax = Apx
                        self.atomA = [Ax,Ay]
        else:
            #the supercell is hexagon shape, this hexagon is horizontal hexagon   
            print("supercell is 'horizontal hexagon' ")
            if Apy > By:
                #type 1 horizontal hexagon A' is higher than B
                Bpx = test2[np.where(test2[:,0]>Bx)][:,0].min()
                Bpy = test2[np.where(test2[:,1]==By)][:,1].max()
                self.atomB = [Bpx,Bpy]
                Cpx = test1[np.where(test1[:,0]<Cx)][:,0].max()
                Cpy = test1[np.where(test1[:,1]<Cy)][:,1].max()
                self.atomC = [Cpx,Cpy]
            else:
                #type 2 horizontal hexagon A' is equal high with B
                Bpx = test2[np.where(test2[:,0]>Bx)][:,0].min()
                Bpy = test2[np.where(test2[:,1]<By)][:,1].max()
                self.atomB = [Bpx,Bpy]
                Cpx = test1[np.where(test1[:,0]<Cx)][:,0].max()
                Cpy = Cy #in this case, Cpy = Cy = Bpy, so they are equal high
                self.atomC = [Cpx,Cpy]

    def get_vertices(self):
        return self.atomO, self.atomA, self.atomB, self.atomC
    


class Count_atom_num:
    def __init__(self, points, atomo, atomb, atoma, atomc):
        self.points = points
        self.o = atomo
        self.a = atoma
        self.b = atomb 
        self.c = atomc 
        
    def cell_info(self):
        """
        This method will give you the basic information about the Moire cell.
        """
        
        ob = self.b - self.o
        oc = self.c - self.o
        
        dot_product = np.dot(ob, oc)
        norm_a1 = np.linalg.norm(ob)
        norm_a2 = np.linalg.norm(oc)
        cos_theta = dot_product / (norm_a1 * norm_a2)
        theta = np.arccos(np.clip(cos_theta, -1.0, 1.0))  # clip 避免數值誤差超出 [-1,1]
        cell_area = norm_a1 * norm_a2 * np.sin(theta)
        return cell_area, ob, oc, norm_a1
        
    def count_atom_num(self, stack_conf, sublattice_type, h=0):

        polygon = Polygon([self.o, self.b, self.a, self.c])

        # Create temporary data
        temp_points = self.points[(self.points[:, 1] > -0.1) & (self.points[:, 0] > self.c[0] - 0.1) & (self.points[:, 0] < self.b[0] + 0.1)]
        in_spcell = np.array([point for point in temp_points if polygon.contains(Point(point[0], point[1]))])
        
        if stack_conf == 'AA' and sublattice_type == 'A':
            in_spcell = np.vstack((self.o, in_spcell))
        elif stack_conf == 'AB' and sublattice_type == 'B':
            in_spcell = np.vstack((self.o, in_spcell))

        in_spcell = np.hstack((in_spcell, np.full((in_spcell.shape[0], 1), h)))
        
        return in_spcell
    
    def add_corrugation(self, atom_coord):
        """
        先算oa和oc兩線段間的夾角，
        如果是120度:
            o,a,b,c四個角是AA的位置
            1/3*oa - 1/3*oc 是AB的位置，2/3*oa - 2/3*oc是BA的位置
        60度:
            o,a,b,c四個角是AA的位置
            1/3*oa + 1/3*oc 是AB的位置，2/3*oa + 2/3*oc是BA的位置
        
        atom_coord就是所有在cell內的原子座標，由can[1:]取得
        """
        pass
    


    
