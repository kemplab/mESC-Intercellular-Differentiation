import sys
import random as r
from Simulation import *
from simulationObjects import *
import networkx as nx
import numpy as np
import platform
from IntC_FP import IntC_Diff
from CountDirectory import newDirect


def main(x,path,inhib_value,inhib_type,inhib_time,inhib_time1):
    #Start time, end time, time step
    ts = 0
    te = 72
    dt = 1

    #Now make a new simulation
    
    sim = Simulation(x,
                     path,
                     ts,
                     te,
                     dt,
                     inhib_type,
                     inhib_time=inhib_time,
                     inhib_time1=inhib_time1)

    

    

    #Parameters for intercellular diffusion
    #Index 0 represents a pluripotent cell
    #Index 1 represents a differentiated cell
    num_cell_types = 2
  
    # Probability of a channel being open in each cell types
    open_prob=np.zeros([num_cell_types])

    #Number of active channels at interface
    N=1000. ## arbitrary value 
    #Intrinsic permeability of channel to diffusing substance
    Y=5e-9  ## arbitrary value
    Pmax=Y*N
    ## Each cell type has a set of production and degradation parameters
    ## First 3 represent production following an inhibition curve (1-Hill) to represent a decrease in
    ## production when the cell is already at a high level 
    # a,the max_prod or plateau of curve, index 0
    # k, coefficient representing when 50% of max_prod occurs, index 1
    # n, the steepness of the curve (represents cooperation in hill function), index 2
    ## Degradation rate is Index 3, with the cell degrading this percent of its current concentration each time step 
    celltype_param=np.zeros([num_cell_types,4])
    
    #Pluripotent phenotype parameters
    celltype_param[0,0]=4.89e-6
    celltype_param[0,1]=0.03
    celltype_param[0,2]=1.
    celltype_param[0,3]=1.3e-7
     
    open_prob[0]=0.45
    ##Differentiated phenotype parameters
    celltype_param[1,0]=6.0e-6
    celltype_param[1,1]=0.04
    celltype_param[1,2]=1.
    celltype_param[1,3]=1.3e-7
    
    open_prob[1]=0.85
    
    

    fastest_division_time=18.0
    
    D_camp=444 ##um^2/sec diffusion coefficient for cyclic AMP Dworkin & Keller 1977 J. Bio Chem
    
    net=nx.read_gpickle(os.getcwd()+"\\Colony.gpickle")

    num=len(net.nodes())

    maxCells= int(2*num*((te)/fastest_division_time))
    
    division_times=[19,51]
    #Grid dimension (NxNxN), (I'm pretty sure I stopped using this altogether...)
    grid_dim=130.0
    #Grid spacing in um
    ddim=2.0
    # This isn't used
    sync_time=[500,500]
    ## Initializing the intercellular diffusion module
    XI =IntC_Diff("cAMP",maxCells,D_camp,grid_dim,ddim,Pmax,open_prob,
                  celltype_param,division_times,sync_time)
    sim.add_interC(XI)
    nodes=net.nodes()
    dist=range(num)
    r.shuffle(dist)

    
    for i in range(0,num):
        node=nodes[dist[i]]
        ID=i
        point=node.location
        radius=6.    
        div_set=r.random()*19.

          
        sim_obj=StemCell(point,radius,ID,"U",inhib_value=inhib_value,division_set=div_set)
        
        sim.add(sim_obj)
        sim.inc_current_ID()

        try:
            sim.collide()
        except:
            sim.collide_lowDens()

        
    
    sim.run()
  #set the seperator
    if(platform.system() == "Windows"):
        #windows
        sep = "\\"
    else:
        #linux/unix
        sep = "/"   


