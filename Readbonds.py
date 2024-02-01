def read_results(slurm_out, range1, bond_out, range2):
    log_file=open(slurm_out,'r')
    log_df = pd.DataFrame(columns=["Step","PotEng","temp","new bonds","total bonds"])
    for i, line in enumerate(log_file):
        if i > range1[0] and i < range1[1] :
            log_df.loc[len(log_df.index)] = line.split()
    bond_file = open(bond_out,'r')
    
    return log_df#, bond_df
