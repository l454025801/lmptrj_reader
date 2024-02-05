import numpy as np
import pandas as pd
import plotly.graph_objects as go
import os
import numba as nb


# ===============================================================================================================
'''
Functions by norok2 from STACKOVERFLOW
https://stackoverflow.com/questions/66674537/python-numpy-get-difference-between-2-two-dimensional-array
To find different 2d elements in 2 np arrays
'''

#@nb.njit
def mul_xor_hash(arr, init=65537, k=37):
    result = init
    for x in arr.view(np.uint64):
        result = (result * k) ^ x
    return result


#@nb.njit
def setdiff2d_nb(arr1, arr2):
    # : build `delta` set using hashes
    
    delta = {mul_xor_hash(arr2[0])}
    for i in range(1, arr2.shape[0]):
        delta.add(mul_xor_hash(arr2[i]))
    # : compute the size of the result
    n = 0
    for i in range(arr1.shape[0]):
        if mul_xor_hash(arr1[i]) not in delta:
            n += 1
    # : build the result
    result = np.empty((n, arr1.shape[-1]), dtype=arr1.dtype)
    j = 0
    for i in range(arr1.shape[0]):
        if mul_xor_hash(arr1[i]) not in delta:
            result[j] = arr1[i]
            j += 1
    return result

def setdiff2d_bc(arr1, arr2):
    idx = (arr1[:, None] != arr2).any(-1).all(1)
    return arr1[idx]

def setdiff2d_bc(arr1, arr2):
    idx = (arr1[:, None] != arr2).any(-1).all(1)
    return arr1[idx]
# ===============================================================================================================

def get_molID(index, fiber_beads_num=1800, IgG_beads_num=6, IgG_num=10):
    # return the molecule ID of a atom, and the type of the molecule
    # in one unit cell, there are 1 fiber + 10 IgG (1860 beads in run4)
    cell_num = 1+index//(fiber_beads_num + IgG_beads_num*IgG_num)
    index_residual = index%(fiber_beads_num + IgG_beads_num*IgG_num)
    if index_residual <= fiber_beads_num:
        return cell_num*11+1, "fiber"
    else:
        return cell_num*11+1+(index_residual-fiber_beads_num)//IgG_beads_num+1, "IgG"

def cal_bridge_loop(bonds_arr, prev_bonds_arr):
    # calculate how many bridges and loops **formed in a frame
    # take the bonds_arr as input from read_results
    
    single = 0
    bridge = 0
    loop = 0
    
    # only explore new formed bonds
    new_bonds = setdiff2d_bc(bonds_arr, prev_bonds_arr)
    print(new_bonds)
    explored_IgG = []
    for i in range(len(new_bonds)):
        IgG_molID = new_bonds[i][3]
        if IgG_molID in explored_IgG:
            continue
        else:
            explored_IgG.append(IgG_molID)
            tmp_list = np.where(prev_bonds_arr[:,3]==IgG_molID)[0]
            tmp_list2 = np.where(new_bonds[:,3]==IgG_molID)[0]
            if len(tmp_list) == 0 and len(tmp_list2) == 1:
                # if the IgG not in previous list and only appear once in new list
                # it is a single bond
                single += 1
            elif len(tmp_list2) == 2:
                # if it appears twice in new list
                # determine if it is a bridge or loop
                if new_bonds[tmp_list2[0]][1] == new_bonds[tmp_list2[1]][1]:
                    loop += 1
                else:
                    bridge += 1
            elif len(tmp_list) == 1 and len(tmp_list2) == 1:
                # if it appears once in previous frame
                # a new loop/bridge forms and remove one single bonds
                if new_bonds[tmp_list2[0]][1] == prev_bonds_arr[tmp_list2[0]][1]:
                    loop += 1
                    single -= 1
                else:
                    bridge += 1
                    single -= 1
    return np.array([[single, bridge, loop]])
    
def read_results(slurm_out, range1, bond_out):
    '''
    slurm_out: log file containing thermo output
    range1: range of rows contain thermo output
    bond_out: dumped file for bonds
    '''
    log_file=open(slurm_out,'r')
    # create a DF containing all thermo output
    log_df = pd.DataFrame(columns=["Step","PotEng","temp","new bonds","total bonds"])
    for i, line in enumerate(log_file):
        if i > range1[0] and i < range1[1] :
            log_df.loc[len(log_df.index)] = line.split()

    log_df = log_df.astype({'Step': int, 'PotEng': float, 'temp':float, 'new bonds':int, 'total bonds':int})
    log_df2 = log_df[log_df['new bonds'] !=0] # get a DF with steps when new bonds form
    log_file.close()
    bond_file = open(bond_out,'r')
    
    timestep_index=[] # list hold line index for each time step
    bond_file_content = bond_file.readlines()
    for num, line in enumerate(bond_file_content):
            if line == "ITEM: TIMESTEP\n":
                timestep_index.append([num+1, int(bond_file_content[num+1])])
                                      # index of the line,  Timestep in simulation
    
    print(len(timestep_index))
    prev_bonds_arr = np.zeros((1,4))
    cross_link_count = np.zeros((1,3))
    for i,pair in enumerate(timestep_index):
        if pair[1] in log_df2['Step'].values:
            print('prev')
            print(prev_bonds_arr)
            index = pair[0]
            bonds_arr = np.zeros((1,4)) # fiber atom, fiber molID, IgG atom, IgG molID
            for i in range(index+8, timestep_index[i+1][0]-1):
                line_content_list = bond_file_content[i].split()
                if get_molID(int(line_content_list[1]))[1] == 'fiber':
                    bonds_arr = np.vstack([bonds_arr,[int(line_content_list[1]), get_molID(int(line_content_list[1]))[0],\
                                    int(line_content_list[2]), get_molID(int(line_content_list[2]))[0]]])
                else:
                    bonds_arr = np.vstack([bonds_arr,[int(line_content_list[2]), get_molID(int(line_content_list[2]))[0],\
                                    int(line_content_list[1]), get_molID(int(line_content_list[1]))[0]]])
            print('now')
            print(bonds_arr)
            print('new')
            cross_link_count = np.vstack([cross_link_count, cal_bridge_loop(bonds_arr, prev_bonds_arr)])
            print(cross_link_count)
            print('--------------------------')
            prev_bonds_arr = bonds_arr
            #log_df2.loc[]
    bond_file.close()
    
    return log_df, log_df2, cross_link_count

#https://stackoverflow.com/questions/66674537/python-numpy-get-difference-between-2-two-dimensional-array
