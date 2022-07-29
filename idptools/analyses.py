# Some standard analyses, streamlined for IDP's
from functools import wraps
import itertools

import numpy as np
import pandas as pd
import mdtraj
from utility import stats
from utility import log


# ===== Helpers =====
def load_chains(func):
    @wraps(func)
    def loaded_func(traj,top=None,stride=1,atom_indices=None,*args,**kwargs):
        if isinstance(traj,mdtraj.Trajectory):
            t = traj

        elif top is None:
            t = mdtraj.load(traj,stride=stride)
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

def load(func):
    @wraps(func)
    def loaded_func(traj,top=None,stride=1,atom_indices=None,*args,**kwargs):
        if isinstance(traj,mdtraj.Trajectory):
            t = traj 
        elif top is None:
            t = mdtraj.load(traj,stride=stride)
        else:
            t = mdtraj.load(traj,top=top,stride=stride)
        return func(t,*args,**kwargs)
    return loaded_func

@load
def compute_dist_from_com(traj,indices = None, com = None):
    """
    Args:
        traj (mdtraj.Trajectory)
        indices (str or list-like): selection to calculate wrt
        com (str or list-like): selection of reference coords. default to protein selection.
    Returns:
        (np.array): displacements
        (np.array): distances
    Note:
        1) [indices], [indices] --> just use mdtraj
        2) selection str, selection str --> just use mdtraj
        3) 

        should select the centers first
        then broadcast and calculate all the distances and displacements in one shot without for loop!
    """
    disps = []
    dists = []
    if indices is None:
        target_selection = traj.topology.select("is_water and element O")
    elif isinstance(com,str):
        target_selection = traj.topology.select(indices)
    else:
        target_selection = indices
        
    if com is None:
        center_selection = traj.topology.select("protein")
    elif isinstance(com,str):
        center_selection = traj.topology.select(com)
    else:
        center_selection = com

    # calculate com
    if isinstance(center_selection,np.ndarray) and len(center_selection.shape) == 3:
        coms = center_selection
    else:
        coms = mdtraj.compute_center_of_mass(traj.atom_slice(center_selection))

    # calculate distance from com
    disp = traj.atom_slice(target_selection).xyz - coms[:,None,:]
    #disp,dist = compute_min_img(disp, traj.unitcell_lengths[:,None,:])
    disp -= np.rint( disp/traj.unitcell_lengths[:,None,:] )
    dist = np.sqrt( (disp**2.0).sum(-1) )

    return disp,dist
    
def histogram(data,geom="sph",*args,**kwargs):
    """Utility for working with histograms
    Args:
        data (array-like)
        geom (str): sph, cyl, or None for now
    Returns:
        array-like: bin mids
        array-like: bin values, with spherical weights if specified
        array-like: bin edges
    """
    h,b = np.histogram(data,*args,**kwargs)
    
    b_center = 0.5*(b[:-1]+b[1:])
    if isinstance(geom,str):
        if geom.lower() == "sph":
            shell_volumes = 4/3*np.pi*(b[1:]**3.-b[:-1]**3.)
        elif geom.lower() == "cyl":
            shell_volumes = np.pi*(b[1:]**2.-b[:-1]**2.)
        else:
            raise ValueError("Unrecognized geometry argument: {}".format(geom)) 
        h = h/shell_volumes
    
    return b_center, h, b 


def compute_min_img(disps,box):
    """
    Note:
        Todo: handle broadcasting... but would need to know providence and shape of data...
              maybe demand that user needs to make sure they're broadcasted properly?
              for distances, also need to know what dimension to sum over. currently assume is last index.
    """
    disps -= np.rint(disps/box)
    dists = np.sqrt( (disps**2.0).sum(-1) )
    return disps, dists 

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
        simplified (bool): whether to use simplified DSSP codes
    Returns:
        array: dssp codes, by residue
        array: dssp counts for each code type, per frame
        array: statistics of dssp counts
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



@load
def hydration(t,hyd_cut=0.31):
    """Calculate hydration for chains in a system
    Args:
        t (mdtraj Trajectory): with decorator, can also take trajfile, topfile, stride, atom_indices
        hyd_cut (float): 1st hydration cutoff
    Returns:
        array: Rg/frame, [ statistics X nframes X nchains]
        array: statistics of each chain (mean, std, err)
        array: statistics of whole system (mean, std, err)
    Note: 
        criteria: water oxygen within 0.31nm of backbone (Yingling 2014)
    """
    #t_chains = t.atom_slice( t.topology.select("protein") )
    hyd_cut = 0.31 #nm

    backbone = t.topology.select("backbone")
    water_o = t.topology.select("is_water and element O")
    #pairs = list(itertools.product(backbone,water_o))
    pairs = t.topology.select_pairs(backbone,water_o)

    dists = mdtraj.compute_distances(t,pairs)
    dists_hyd = (dists<=0.31)
    dists_hyd_count = dists_hyd.sum(1)

    # calculate statistics
    hyd_stats = stats.stats(dists_hyd_count)
    np.savetxt("hydration_count",dists_hyd_count)
    s = stats.format_stats(*hyd_stats,varname="hydration")
    np.savetxt("hydration_stats",np.array(hyd_stats),header="#mean\tstd\terr\tt0\tg\tNeff")     

    # diagnostic, g(r)
    r,gr = mdtraj.compute_rdf(t,pairs,r_range=[0.,2.0])
    rgr  = np.vstack([r,gr])
    np.savetxt("hydration_gr.txt",rgr.T,header="#r\tg(r)")
    

    return dists_hyd_count,hyd_stats,r,gr


@load
def water_from_protein_com(t,SIunits=False,*args,**kwargs):
    """Calculate distribution of water from given c.o.m.
    Args:
        t (mdtraj.Trajectory)
        SIunits (bool): whether to use g/cm^3 SI units instead of #/nm^3
        *args (misc): arguments for histogram
    Note:
        if using mdtraj, returning density in units of #/nm^3
    """
    backbone = t.topology.select("backbone")
    water_o = t.topology.select("is_water and element O")     
    disps,dists = compute_dist_from_com(t,indices = water_o,com = backbone)

    r,gr,edges = histogram(dists,geom="sph",*args,**kwargs)
    gr /= t.n_frames
    if SIunits:
        gr *= 18.01528/6.02214076/100
        header = "#r(nm)\tgr(g/cm^3)"
    else:
        header = "#r(nm)\tgr(#/nm^3)"

    print(dists)    

    data = np.vstack([r,gr]).T
    np.savetxt("gr_water_from_protein_com.txt",data,header=header)

    return r,gr



