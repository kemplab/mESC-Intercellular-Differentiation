import numpy as np
from extractMetrics import *

import os, sys
import networkx as nx




prepath="D:\\Final FP Sims\\Base Model (exp)\\"
prepath=os.getcwd()
ts=30
te=31
dt=1
#Set time range of simulations 
qnum=range(ts,te,dt)
#Set range of simulation folders 
for fold_num in range(1,2,1):
    path=prepath+repr(fold_num)+".0\\"
    path=prepath+"\\"
    da=[]
    for q in qnum:
        temp=path+"network"+repr(q)+".0.gpickle"
        network = nx.read_gpickle(temp)
        dat_anl=extractMetrics(network)
        print dat_anl
        da.append(dat_anl)
        da_temp=np.asarray(da)
        np.save(path+"sim_pca.npy",da)
        print "Network ", q        
    da=np.asarray(da)
##    np.save(path+"sim_pca.npy",da)






