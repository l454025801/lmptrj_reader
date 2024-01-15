from io import StringIO
import numpy as np
import pandas as pd

def read_lmptrj(lmptrj, steps, frames, trj_format):
    '''
        lmptrj: lammps dump trajectory
    trj_format: List. the format of the dump file. Ex: ["index", "type", "x", "y", "z"]
    '''
    lmptrj_file = open(lmptrj,'r')
    
    # first store each frame as a key in a dictionary
    tmp_dict={}
    frame = -1
    tmp_txt = ''
    for line in lmptrj_file.readlines():
        if line=="ITEM: TIMESTEP\n":
            print('*', end='')
            if frame != -1:
                frame_array = np.loadtxt(StringIO(tmp_txt),skiprows=9)
                tmp_dict[frame] = frame_array
            frame += 1
            tmp_txt=''
        tmp_txt+=line
    # last frame 
    frame_array = np.loadtxt(StringIO(tmp_txt),skiprows=9)
    tmp_dict[frame] = frame_array
    # release memory
    del(tmp_txt)
    
    
    lmptrj_file.close()
    return tmp_dict
