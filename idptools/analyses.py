# Some standard analyses, streamlined for IDP's
from functools import wraps
import numpy as np
import pandas as pd
import mdtraj
from utility import stats
from utility import log


# ===== Helpers =====
def load_chains(func):
    def loaded_func(traj,top=None,stride=1,atom_indices=None,*args,**kwargs):
        if top is None:
            t = mdtraj.load(traj)
        else:
            tmp = mdtraj.load(top,stride=stride)
            chain_selection = tmp.topology.select("protein")    

            if atom_indices is not None:
                t = mdtraj.load(traj,top=top,atom_indices=atom_indices,stride=stride) 
            elif chain_selection is None:
                t = mdtraj.load(traj,top=top,stride=stride)
            else:
                t = mdtraj.load(traj,top=top,atom_indices=chain_selection,stride=stride)
        return func(t,*args,**kwargs)
    return loaded_func

# ===== Analyses =====
# Rg
@load_chains
def rg(t,prefix=None):
    """Calculate Rg for chains in a system
    Args:
        t (mdtraj Trajectory): with decorator, can also take trajfile, topfile, stride, atom_indices
    Returns:
        array: Rg/frame, [nframes X nchains]
        array: statistics of each chain (mean, std, err)
        array: statistics of whole system (mean, std, err)
    Note: 
        also saves files
    """
    """
    # load
    if top is None:
        t = mdtraj.load(traj)
    else:
        t = mdtraj.load(traj,top=top)
    """
    t_chains = t.atom_slice( t.topology.select("protein") )
    # get Rg for each chain
    rgs = []
    s = []
    for c in t_chains.topology.chains:
        inds = [a.index for a in c.atoms]
        rg = mdtraj.compute_rg(t.atom_slice(inds))
        rgs.append(rg)
        stat = stats.stats(rg)
        s.append(stat)

    # save
    rgs = np.array(rgs).transpose()
    np.savetxt("rgs.dat",rgs)

    report ="Rg report:\n"
    report+="==========\n"
    for ii,stat in enumerate(s):
        report+=">>> Chain {}\n".format(ii)
        line = stats.format_stats(*stat)
        report += line
    with open("rgs.txt","w") as f:
        f.write(report)    

    return rgs,s


def ree(traj):
    pass    


@load_chains
def dssp(t,simplified=True):
    """Calculate Rg for chains in a system
    Args:
        t (mdtraj Trajectory): with decorator, can also take trajfile, topfile, stride, atom_indices
    Returns:
        array: Rg/frame, [ statistics X nframes X nchains]
        array: statistics of each chain (mean, std, err)
        array: statistics of whole system (mean, std, err)
    Note: 
        also saves files
        TODO: handle multiple chains
    """
    t_chains = t.atom_slice( t.topology.select("protein") )
    if simplified:
        #codes_simplified = ["H","E","C"]
        codes = ["H","E","C"]
    else:
        codes = ["H","B","E","G","I","T","S",""]

    # collect data
    dssps = []
    dssp_summaries = []
    dssp_summary_stats = []
     
    df_dssps = pd.DataFrame()
    df_dssp_summaries = pd.DataFrame()
    df_dssp_summary_stats = pd.DataFrame()
    for ic,c in enumerate(t_chains.topology.chains):
        if t_chains.n_chains == 1:
            prefix = "dssp"
        else:
            prefix = "dssp_chain{:02g}".format(ic)

        # compute
        inds = [a.index for a in c.atoms]
        dssp = mdtraj.compute_dssp(t.atom_slice(inds),simplified=simplified)
        dssps.append(dssp)
        df = pd.DataFrame(dssp)
        df.to_csv(prefix+".csv")
        df["chain"] = ic
        if ic == 0:
            df_dssps = df.copy()
        else:#collate
            df_dssps = pd.concat([df_dssps,df]) 
        
        #get summary counts 
        summary = [] 
        for code in codes:
            n_count = (dssp==code).sum(1)
            summary.append(n_count)
        summary = np.array(summary).T
        dssp_summaries.append(summary)
        df = pd.DataFrame(summary,columns=codes)
        df.to_csv(prefix+"_summaries.csv")
        df["chain"] = ic
        if ic == 0:
            df_dssp_summaries = df.copy()
        else:
            df_dssp_summaries = pd.concat([df_dssp_summaries,df])

        #get statistics
        statistics = []
        for ii,code in enumerate(codes):
            s = stats.stats(summary,col=ii)
            statistics.append(s)
        statistics = np.array(statistics).T
        dssp_summary_stats.append(statistics)
        df = pd.DataFrame(statistics,columns=codes,index=["mean","std","err","t0","g","Neff"])
        df.to_csv(prefix+"_summary_stats.csv") 
        df["chain"] = ic
        if ic == 0:
            df_dssp_summary_stats = df.copy()
        else:
            df_dssp_summary_stats = pd.concat([df_dssp_summary_stats,df])
            
        
    # save data
    # collate all into one massive dataframe, indexed by chain & frame???
    df_dssps.to_csv("dssp_all.csv")
    df_dssp_summaries.to_csv("dssp_summaries_all.csv")
    df_dssp_summary_stats.to_csv("dssp_summary_stats_all.csv")

    #alternatively: return pandas collated dataframes
    return dssps, dssp_summaries, dssp_summary_stats




