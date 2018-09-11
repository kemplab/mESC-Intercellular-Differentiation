import numpy as np
import networkx as nx
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt 
from simulationMath1 import Distance
from sklearn.preprocessing import scale
from sklearn import preprocessing 
import random as rand
import os
import math


#Set path to folder containing simulations:



np.set_printoptions(suppress=True)

X1=np.load(os.getcwd()+'\\TrainSet.npy')
scaler = preprocessing.StandardScaler().fit(X1)
p_inds=[0,120,240,360,480,600,720,840,960]
IDs=[]
pat_names=['Undifferentiated','Random','Snaked','Globular','Inv Snake','Inside-Out','Outside-In','Differentiated']
col1=['c','m','g','b','r',[0.7,0.7,0.7],'y','k']
plot_labs=[]

for i in range(len(pat_names)):
    for j in range(p_inds[i],p_inds[i+1]):
        IDs.append(col1[i])
        plot_labs.append(pat_names[i])
        




X1=scaler.transform(X1)
pca=PCA(n_components=3)
pca.fit(X1)

X_r=pca.transform(X1)


##Uncomment below for running PCA on simulation files
"""
sim_path="C:\\..."
sims=[]
#The range of simulation folders you wish to plot
for i in range(1,10):
    temp=np.load(sim_path+repr(i)+".0\\sim_pca.npy")
    temp=scaler.transform(temp)
    sims.append(temp)

sims=np.asarray(sims)
S_r=np.zeros([len(sims),len(sims[0]),3])
for i in range(len(sims)):
    temp=pca.transform(sims[i,:])
    S_r[i,:,:]=temp

S_avg=np.mean(S_r,axis=0)
"""

fig=plt.figure(figsize=(21,10),facecolor='w')
##Plot PC1 vs PC2

plt.subplot(1,2,1)
ax=plt.gca()
#Training Set
for i in range(len(X_r)):
    if i in p_inds:
        ax.scatter(X_r[i, 0], X_r[i, 1],c=IDs[i],s=50.0,alpha=.3,label=plot_labs[i])
    else:
        ax.scatter(X_r[i, 0], X_r[i, 1],c=IDs[i],s=50.0,alpha=.3)
plt.legend()
"""
#Plot simulations converted to LVS
ax.plot(S_avg[:,0],S_avg[:,1],c='k',linewidth=3.)
"""

##Plot PC1 vs PC3
plt.subplot(1,2,2)
ax=plt.gca()
#Training Set
for i in range(len(X_r)):
    ax.scatter(X_r[i, 0], X_r[i, 2],c=IDs[i],s=50.0,alpha=.3)
"""
#Plot simulations converted to LVS
ax.plot(S_avg[:,0],S_avg[:,1],c='k',linewidth=3.)
"""
plt.show()




