import numpy as np
import random as r
from math import atan2,cos,sin
import os, sys
import networkx as nx

import matplotlib.pyplot as plt
from Simulation import *
from simulationObjects import *
import simulationMath1 as simMath
from simulationMath1 import Distance as dist
from region import region
import matplotlib.pyplot as plt


"""
Contains the functions for analyzing a .gpickle and calculating
the 7 metrics used for PCA 


"""
def GetCenter(agents):
    #gets the center of a list of agents in a network
    
    cent = [0,0]
    for i in range(0, len(agents)):
        try:
            cent = simMath.AddVec(cent, agents[i].pos)
        except:
            cent = simMath.AddVec(cent, agents[i].location)
    #Normalize
    cent = simMath.ScaleVec(cent, 1.0/len(agents))
    #now return
    return cent
def getAverageRadialDistance(agents, center):
    
    cent = [0,0,0]
    radii = []
    for i in range(0, len(agents)):
        radius = simMath.SubtractVec(agents[i].pos, center)
        #add tot he list
        radii.append(simMath.Mag(radius))
    #now return the average and the stdev
    return np.average(radii), np.std(radii)    

def GetRadiusandCenter(agents):
    #first get the center
    center = GetCenter(agents)
    #gets the maximum radial size of the network
    
    radius = -1
    for i in range(0, len(agents)):
        r = simMath.SubtractVec(center, agents[i].pos)
        radius = max(simMath.Mag(r), radius)
    #now return
    return radius, center
def closest_nbs(network):
    nodes=network.nodes()
##    nbs = np.zeros((len(nodes),len(nodes)))
    nbs=np.zeros((len(nodes),2))
    
    for i in range(len(nodes)):
        nbs[i,0]=i
        for j in range(i+1,len(nodes)):
            
            obj1=nodes[i]
            obj2=nodes[j]
            dist_vec = simMath.SubtractVec(obj2.pos, obj1.pos)
            #get the magnitude
            dist = simMath.Mag(dist_vec)
            if dist<35:
                nbs[i,1]=nbs[i,1]+1
                nbs[j,1]=nbs[j,1]+1
    for i in range(len(nodes)):
        node=nodes[i]
        node.numNbs=nbs[i,1]
    nodes2=sorted(nodes,key=lambda nodes:nodes.numNbs)
    p=int(0.35*len(nodes2))
    outside_nodes=nodes2[0:p]
    
    inside_nodes=nodes2[len(nodes2)-p:len(nodes2)]
   
    return outside_nodes, inside_nodes


def extractMetrics(network):
    limit=0.6
    mu_lim=0.05
    diffCells=0
    nodes=network.nodes()
    subNet=nx.Graph()
    node_ref=dict()
    nP=len(nodes)
    center=GetCenter(nodes)
    
    
    data=np.zeros([nP,2])
    for i in range(nP):
        nb=nodes[i]
        if nb.state=="D":
            diffCells=diffCells+1.
        ind=int(nb.ID)
        node_ref[ind]=i
        reg=region(i)
        subNet.add_node(reg)

    regs=subNet.nodes()
    for i in range(nP):
        node=nodes[i]
        reg=regs[i]
        reg.savePos(node.location[0],node.location[1])
        if node.state=="U":
            val=0
        else:
            val=1
        nbs=network.neighbors(node)
        for j in range(len(nbs)):
            nb=nbs[j]
            
            ind=int(nb.ID)
            
            subNet.add_edge(reg,regs[node_ref[ind]])
            if nb.state=="D":
                val=val+1.
        num=len(nbs)+1
  
        reg.changeValue(val/num,num)
        data[i,0]=num
        data[i,1]=val


    cell_bodies=nx.connected_components(subNet)
    if len(cell_bodies)>1:
        for i in range(1,len(cell_bodies)):
            clust=cell_bodies[i]
            for j in range(len(clust)):
                nb=clust[j]
                ind=int(nb.ID)
                subNet.remove_node(regs[node_ref[ind]])
            

    mu=np.average(data[:,1]/data[:,0],weights=data[:,0])
    

    edges=subNet.edges()
    high_nodes=[]
    low_nodes=[]
    for node in subNet.nodes():
        if node.diffPer>limit or (mu>mu_lim and node.diffPer>mu): 
            high_nodes.append(node)
            node.change_state(1)
        else:
            low_nodes.append(node)
    
    

    D=subNet.copy()
    UD=subNet.copy()
    nodes=subNet.nodes()

    d_edges=D.edges()
    ud_edges=UD.edges()
    for i in range(len(d_edges)):
        edge=d_edges[i]
        edge2=ud_edges[i]
        obj1=edge[0]
        obj2=edge[1]
        obj3=edge2[0]
        obj4=edge2[1]
        if (obj1.state==0 or obj2.state==0):
            D.remove_edge(obj1,obj2)

        if (obj3.state==1 or obj4.state==1):
            UD.remove_edge(obj3,obj4)

    outside_nodes, inside_nodes=closest_nbs(subNet)
    
    inout_ratios=np.zeros([2,1])
    for i in range(len(inside_nodes)):
        node=inside_nodes[i]
        if node.diffPer>limit or (mu>mu_lim and node.diffPer>mu):
            
            inout_ratios[0]=inout_ratios[0]+1.
    for i in range(len(outside_nodes)):
        node=outside_nodes[i]
        if node.diffPer>limit or (mu>mu_lim and node.diffPer>mu): 
            inout_ratios[1]=inout_ratios[1]+1.

    dat_anal=np.zeros(7)

    dat_anal[0]=inout_ratios[0]/len(inside_nodes)
    dat_anal[1]=inout_ratios[1]/len(outside_nodes)
    dat_anal[6]=dat_anal[1]/(0.75+dat_anal[0])

    p= nx.connected_components(D)
    pud=nx.connected_components(UD)
    rad_avg = []
    for j in range(len(pud)):
        if pud[j][0].state==0:
            r, com = GetRadiusandCenter(pud[j])
            #laos get the number of nodes
            avg, std = getAverageRadialDistance(pud[j], com)
            rad_avg.append(avg)


    ud_crd_avg = np.average(rad_avg)
##################
    
    if(np.isnan(ud_crd_avg)):
        ud_crd_avg = 0
    rs = []
    num_nodes = []

    rej_clust=0
    rej_clustUD=0

    clusters=[]
    clustersUD=[]
    for i in range(len(p)):
        group=p[i]
        if len(group)>3:
            clusters.append(group)
        elif group[0].state==1:
            rej_clust=rej_clust+1
            
    for i in range(len(pud)):
        group=pud[i]
        if len(group)>3:
            clustersUD.append(group)
        elif group[0].state==0:
            rej_clustUD=rej_clustUD+1

    
    edges=D.edges()
    edgesUD=UD.edges()
    num_nodes=float(len(subNet.nodes()))
    
    for i in range(len(clustersUD)):
    
        
        clust_UD=nx.Graph()
        clust_UD.add_nodes_from(clustersUD[i])
        group=clustersUD[i]
        
        for j in range(len(group)):
            for k in range(j,(len(group))):
                cell1=group[j]
                cell2=group[k]
                cell_dist=dist(cell1.pos, cell2.pos)
                dat_anal[5]=max(dat_anal[5],cell_dist)
        for j in range(len(edgesUD)):
            edge=edgesUD[j]
            obj1=edge[0]
            obj2=edge[1]

            if obj1 in clustersUD[i] and obj2 in clustersUD[i]:
                clust_UD.add_edge(obj1,obj2)
        x=nx.average_shortest_path_length(clust_UD)
        dat_anal[2]=max(x,dat_anal[2])

    
    dat_anal[2]=dat_anal[2]/nx.average_shortest_path_length(subNet)


    dat_anal[3]=float(rej_clust/(num_nodes))

    dat_anal[4]=float(ud_crd_avg/num_nodes)
    dat_anal[5]=float(dat_anal[5]/num_nodes)
    return dat_anal

    






