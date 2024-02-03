import numpy as np
import pandas as pd
import plotly.graph_objects as go
import os

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
    
    log_df2 = log_df[log_df['new bonds'] !="0"] # get a DF with steps when new bonds form
    log_file.close()
    bond_file = open(bond_out,'r')
    
    line_index=[] # list hold index for each time step
    bond_file_content = bond_file.readlines()
    for num, line in enumerate(bond_file_content):
            if line == "ITEM: TIMESTEP\n":
                line_index.append(num)
    for index in line_index:
        if log_df.loc[log_df["Step"]==bond_file_content[index+1].split()[0],'new bonds'].iloc[0] != "0":
            
    bond_file.close()
    
    return log_df#, bond_df
