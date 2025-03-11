# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 15:50:53 2025

@author: User
"""

import numpy as np
import time

class TBModel:
    def __init__(self, atomo, atomb, atoma, atomc, kx, ky, can=None, sp_zm=0, b_mag=0, beta_d=0, strain_m=1):
        # initialized the input parameters
        self.atomo = atomo
        self.atomb = atomb
        self.atoma = atoma
        self.atomc = atomc
        self.kx = kx
        self.ky = ky
        self.sp_zm = sp_zm
        self.b_mag = b_mag
        self.beta_d = beta_d
        self.strain_m = strain_m
        self.cans = can  # 可以是一些外部計算的結果
        
        # Supercell related parameters
        l1 = np.array(atomb) - np.array(atomo)
        l2 = np.array(atomc) - np.array(atomo)
        self.L1 = np.array([l1[0], l1[1], 0]) 
        self.L2 = np.array([l2[0], l2[1], 0])
        self.noa = self.cans[1].shape[0]  # number of atoms in layer1_sublattice A, for constructing Hamiltonian
        self.nol = len(self.cans) - 1     # number of layer or type of sublattice  of the system
        
        # Basic parameters appear in the Hamiltonian
        self.a0 = 1.42         #Å, bond length
        self.d0 = 3.35         #Å, interlayer distance
        self.Vpppi0 = -2.7     #eV, default = -2.7
        self.Vppsg0 = 0.48     #eV
        self.delta  = 0.184*2.46  #Å, decay length of the transfer integral
        self.cutoff_r = 4*self.a0+0.01  #cutoff distance for hopping, is d=4*a0 parameter in Koshino's paper
        self.Ham = np.zeros((self.kx.shape[0], self.noa*self.nol, self.noa*self.nol), dtype=complex)  # Hamiltonian matrix
        
    '''
    hopping energy formula:
        
    -t(Ri,Rj) = Vpppi*[1-((Ri-Rj)•ez/cutoff_r)^2] + Vppsg*[((Ri-Rj)•ez/cutoff_r)^2]
    
    Vpppi = Vpppi0*exp(-(cutoff_r-a0)/delta)
    Vppsg = Vppsg0*exp(-(cutoff_r-d0)/delta)
    
    H = 𝚺ij t(Ri,Rj)*|𝛙i><𝛙j| + h.c.
    
    '''
    
    def hopping_t(self, Ri, Rj, RiRj):
        #RiRj  = np.linalg.norm(Ri-Rj) # is equal to parameter "d" in Koshino's paper
        ez = np.array([0,0,1])
        Vpppi = self.Vpppi0 * np.exp(-(RiRj - self.a0) / self.delta)
        Vppsg = self.Vppsg0 * np.exp(-(RiRj - self.d0) / self.delta)
        if RiRj > self.cutoff_r:
            t = 0
        elif abs(RiRj - 0) < 0.0001:
            t = 0     # if onsite = 0 means no onsite
        else:
            if self.strain_m == 1:
                t = Vpppi*(1-(np.dot((Ri-Rj),ez)/RiRj)**2) + Vppsg*((np.dot((Ri-Rj),ez)/RiRj)**2)
            else:
                #this term is strain term: np.exp(-beta_d*(RiRj_strain/a0-1))
                RiRj_s = np.linalg.norm((Ri-Rj)*np.array([self.strain_m, self.strain_m,1]))
                t = Vpppi*(1-(np.dot((Ri-Rj),ez)/RiRj)**2)*np.exp(-self.beta_d*(RiRj_s/self.a0-1)) + Vppsg*((np.dot((Ri-Rj),ez)/RiRj)**2)
        return t

    
    @staticmethod  #設為靜態方法(static method)
    def onsite_e_field(noa, ls_index_e, onsite_list): 
        """
         通用的 onsite electric field 函數。
        
         Args:
         noa (int):        每層的 sublattice 原子數
         ls_index_e (int): 對應的 sublattice 種類數
         onsite_list (list or array): 每層的 onsite 參數清單。
        
         Returns:
         np.ndarray: 對應的對角矩陣。
        """

        diagnal = np.hstack([np.full(noa * ls_index_e, onsite) for onsite in onsite_list])
        diagonal_matrix = np.diag(diagnal)
        return diagonal_matrix
    
    
    @staticmethod
    def onsite_stagger_potential(noa, onsite_a, onsite_b, num_layers): 
        """
         通用的 staggered potential 函數，適用於任意層數。
        
         Args:
         noa (int): 每層的 sublattice 原子數
         onsite_a (float): sublattice a 的位能
         onsite_b (float): sublattice b 的位能
         num_layers (int): 層數
        
         Returns:
         np.ndarray: 對應的對角矩陣。
        """

        onsite_l_a = np.full(noa, onsite_a)
        onsite_l_b = np.full(noa, onsite_b)
        
        # Combine the onsites of each layer
        diagnal = np.hstack([onsite_l_a if i % 2 == 0 else onsite_l_b for i in range(2 * num_layers)])
        
        # Create a diagonal matrix
        diagonal_matrix = np.diag(diagnal)
        return diagonal_matrix
    

    def construct_ham(self, ls_index_s, ls_index_h, hopping_list=False): # hopping for atoms in layer_i sublattice_alpha
        #ls_index_s stands for index for select_atom cans[i], where i = 1,2,3,4...
        #ls_index_h stands for index for hopping cans[i], where i = 1,2,3,4...
        '''
        new_cans = np.zeros((9,noa,3))
        for x in range(-1,2,1):     # x times L1
            for y in range(-1,2,1): # y times L2
                #for i in range(cans[ls_index_h].shape[0]):   # shift layer1 sublattice A's atom
                new_cans[(x+1)*3+(y+1)] = cans[ls_index_h] + x*L1 + y*L2
        #print('new_i=',new_i)
        '''
        
        '''
        lambda x: None 是一個匿名函數，它接受一個參數 x，但什麼都不做，只返回 None。
        這樣，當我們在迴圈中調用 save_action(result) 時，不會有任何動作，省去了不必要的操作。
        '''
        saved_data = []
        save_action = saved_data.append if hopping_list else lambda x: None
        
        # 
        x_values = np.arange(-1, 2)
        y_values = np.arange(-1, 2)
        
        # 生成所有 x 和 y 值的組合
        X, Y = np.meshgrid(x_values, y_values)
        shifts = np.column_stack((X.ravel(), Y.ravel()))
        
        # 利用廣播機制進行向量化計算
        new_cans = self.cans[ls_index_h] + shifts.dot(np.vstack((self.L1, self.L2)))[:, np.newaxis]
        
        start_time = time.perf_counter()
        
        for k in range(self.noa):    # k is the selected atom
            select_atom = self.cans[ls_index_s][k]   
            
            for ii in range(9): #The 3x3 grid consists of 9 cells, and from (-1,-1), (-1,0) ... to (1,1), objects that can hop are selected.
                #print(k,ii)
                
                for i in range(self.noa):
                    #print(new_cans[ii][i])
                    distance = np.linalg.norm(select_atom - new_cans[ii][i]) #is equal to vector "d" in Koshino's paper

                    if distance < self.cutoff_r:
                        Rij = select_atom - new_cans[ii][i] 
                        hopping_energy = self.hopping_t(select_atom,new_cans[ii][i],distance)
                        for kxy in range(self.kx.shape[0]):
                            self.Ham[kxy][k+self.noa*(ls_index_s-1)][i+self.noa*(ls_index_h-1)] += hopping_energy*np.exp(1j*(Rij[0]*self.kx[kxy]+Rij[1]*self.ky[kxy]))
                        
                        save_action([Rij[0], Rij[1], Rij[2], 
                                     (ls_index_s-1) * self.noa + k, 
                                     (ls_index_h-1) * self.noa + i, 
                                     hopping_energy])

                            
        end_time = time.perf_counter()
        print(f"({ls_index_s})({ls_index_h}) Elapsed time: {end_time - start_time:.6f} s")
        return saved_data if hopping_list else None
        
    '''
    需要測試 spin_pn 不算完全修改完成
    '''
    def spin_pn(self,h,b): # h is Hamiltonian
        pauli_z = np.array([[1, 0], [0,-1]]) # Pauli matrix, sigma_z, diagonal_multiplier
        pauli_i = np.array([[1, 0], [0, 1]]) # Identity matrix, off_diagonal_multiplier

        diagonal = np.diag(h)  # get diagonal term, an 1d array
        diagonal = diagonal*np.identity(self.noa*self.nol) # generate diagonal matrix using diagonal term
        off_diagonal = h - diagonal*np.identity(self.noa*self.nol) # generate off-diagonal matrix

        new_dia = np.kron(h, b*pauli_z)     # all the diagonal term times pauli_sigma_z matrix
        new_off_dia = np.kron(h, pauli_i)   # all the off-diagonal term times identity matrix
        new_matrix = new_dia + new_off_dia  # combine two matrix as unity
        return new_matrix
    
    def open_spin(self): #open spin dependent onsite term
        if self.sp_zm != 0:  
            self.Ham = self.spin_pn(self.Ham,self.sp_zm)
    
    def zeeman(self, h, b_m):    
        pauli_i = np.array([[1, 0], [0, 1]]) # Identity matrix, off_diagonal_multiplier
        new_ham = np.kron(pauli_i, h[:])
        N = new_ham[1].shape[0] // 2
        added_zeeman_matrix = np.zeros((2*N,2*N))
        added_zeeman_matrix[:N, :N] += b_m * np.eye(N)
        added_zeeman_matrix[N:, N:] -= b_m * np.eye(N)
        new_ham = new_ham + added_zeeman_matrix
        return new_ham
    
    def open_zeeman(self): #open magnetic field for Zeeman term
        if self.b_mag != 0:  
            self.Ham = self.zeeman(self.Ham, self.b_mag)
    
    
    def finalize_ham(self, hopping_list=False):  # 新的方法來處理 Ham
        
        hopping_lists = []
    
        for ls_index_s in range(1, 5):
            for ls_index_h in range(1, 5):
                if hopping_list==True:
                    hopping_lists.append(self.construct_ham(ls_index_s, ls_index_h, hopping_list))
                else:
                    self.construct_ham(ls_index_s, ls_index_h, hopping_list)
        
        # 增加onsite或stagger potential改成由實例化後從外部控制是否增加
        
        if self.Ham is not None:
            Ham_t = self.Ham.conj().transpose(0, 2, 1)
            self.Ham += Ham_t
        
        #print(self.Ham)  # check Hamiltonian
        if hopping_list:
            return self.Ham, hopping_lists 
        else:
            return self.Ham
