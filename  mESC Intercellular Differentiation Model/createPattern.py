from Simulation import *
import sys
import random as r
from simulationObjects import *
import simulationMath1 as simMath
import networkx as nx
import platform
import os
from Analysis import GetRadiusandCenter as grc
from Analysis import ZeroPad
from simulationMath1 import Distance,RandomPointOnSphere
import numpy as np
from math import atan2
    
##path=os.getcwd()
##struct_path=path+"\\network.gpickle"
##path="C:\\Users\\cglen3\\Desktop\\Patterns\\"
path="C:\\Users\\cglen3\\Desktop\\data_gpickles\\"
spath="C:\\Users\\cglen3\\Desktop\\New Pattern Networks\\Asym Patterns 4\\"
net=nx.read_gpickle(path+"TS2.gpickle")
##net=nx.read_gpickle(os.getcwd()+"\\new.gpickle")
##net=nx.read_gpickle("C:\\Users\\cglen3\\Desktop\\Pattern Networks\\undiff01.gpickle")
rad,cent=grc(net)
nodes=net.nodes()
n=len(nodes)
pn=39
ts=76
te=150
while ts<=te:



    ##Random
    thresh=0.2+0.25*r.random()
    sname="rand"
    for i in range(n):
        agent=nodes[i]
        agent.state="U"
        prob=r.random()
        if prob <thresh:
            agent.state="D"

    ##Full Diff
##    sname="diff"
##    var=1.-0.08*r.random()
####    var=0.92
##    for i in range(n):
##        agent=nodes[i]
##        prob=r.random()
##        if prob<var:
##            agent.state="D"

            
    ##No Diff
##    sname="undiff"
##    var=0.02*r.random()
####    var=0.01
##    for i in range(n):
##        agent=nodes[i]
##        agent.state="U"
##        prob=r.random()
##        if prob <var:
##            agent.state="D"


    

##    #globular

##    sname="glob"
##    numGlob=12
##    globWidth=6.5*1.75
##    var=6.5*1.
##    globs=np.zeros([3,numGlob])
##    for i in range(numGlob):
##        cell_ind=int(r.random()*n)
##        agent=nodes[cell_ind]
##        globs[:,i]=agent.location
##    for i in range(n):
##        agent=nodes[i]
##        agent.state="U"
##        point=agent.location
##        for j in range(numGlob):
##            d=Distance(point,globs[:,j])       
##            if d<(globWidth+r.random()*var):
##                agent.state="D"
####    
    def rnd_snakePoint(p1,p2,frac):
        dx=p2[0]-p1[0]
        dy=p2[1]-p1[1]
        angl=atan2(dy,dx)
        #1/2 or less for frac usually 
        theta=angl-(math.pi*frac)+2.*(math.pi*frac)*r.random()
        x = math.cos(theta)
        y = math.sin(theta)
        return np.array((x,y,0))
    

########    #snaked
##    sname="snake"
##    numGlob=30
##    globWidth=6.*2.
##    point=(rad)*RandomPointOnSphere(2)+cent
##    point1=(rad)*RandomPointOnSphere(2)+cent
##    globs=np.zeros([3,numGlob])
##    globs2=np.zeros([3,numGlob])
##    point2=(2.*globWidth)*rnd_snakePoint(point,cent,0.3)+point
##    point3=(2.*globWidth)*rnd_snakePoint(point1,cent,0.3)+point1
##    globs[:,0]=point
##    globs[:,1]=point2
##    globs2[:,0]=point1
##    globs2[:,1]=point3
##    for i in range(2,numGlob):
##        
##        point=(2*globWidth)*rnd_snakePoint(globs[:,0],globs[:,1],r.random()*0.5)
##        globs[:,i]=point+globs[:,i-1]
##        point2=(2*globWidth)*rnd_snakePoint(globs2[:,0],globs2[:,1],r.random()*0.5)
##        globs2[:,i]=point2+globs2[:,i-1]
##    for i in range(n):
##        agent=nodes[i]
##        agent.state="U"
##        prob=r.random()
########        if prob<0.05:
########            agent.state="D"
##        point=agent.location
##        for j in range(numGlob):
##            d=Distance(point,globs[:,j])
##            d1=Distance(point,globs2[:,j])
##            if d<(globWidth+10*r.random()) or d1<(globWidth+10*r.random()):
##                agent.state="D"

    
##    #inverse snaked
##    sname="inv_snake"
##    numGlob=30
##    globWidth=6.*2.
##    point=(rad)*RandomPointOnSphere(2)+cent
##    point1=(rad)*RandomPointOnSphere(2)+cent
##    globs=np.zeros([3,numGlob])
##    globs2=np.zeros([3,numGlob])
##    point2=(2.*globWidth)*rnd_snakePoint(point,cent,0.3)+point
##    point3=(2.*globWidth)*rnd_snakePoint(point1,cent,0.3)+point1
##    globs[:,0]=point
##    globs[:,1]=point2
##    globs2[:,0]=point1
##    globs2[:,1]=point3
##    for i in range(2,numGlob):
##        
##        point=(2*globWidth)*rnd_snakePoint(globs[:,0],globs[:,1],r.random()*0.5)
##        globs[:,i]=point+globs[:,i-1]
##        point2=(2*globWidth)*rnd_snakePoint(globs2[:,0],globs2[:,1],r.random()*0.5)
##        globs2[:,i]=point2+globs2[:,i-1]
##    for i in range(n):
##        agent=nodes[i]
##        agent.state="D"
##        prob=r.random()
####        if prob<0.05:
####            agent.state="D"
##        point=agent.location
##        for j in range(numGlob):
##            d=Distance(point,globs[:,j])
##            d1=Distance(point,globs2[:,j])
##            if d<(globWidth+10*r.random()) or d1<(globWidth+10*r.random()):
##                agent.state="U"
    def closest_nbs(network,nb_dist,perc):
        nodes=network.nodes()
    ##    nbs = np.zeros((len(nodes),len(nodes)))
        nbs=np.zeros((len(nodes),2))
        
        for i in range(len(nodes)):
            nbs[i,0]=i
            for j in range(i+1,len(nodes)):
                
                obj1=nodes[i]
                obj2=nodes[j]
                dist_vec = simMath.SubtractVec(obj2.location, obj1.location)
                #get the magnitude
                dist = simMath.Mag(dist_vec)
                if dist<nb_dist:
                    nbs[i,1]=nbs[i,1]+1
                    nbs[j,1]=nbs[j,1]+1
        for i in range(len(nodes)):
            node=nodes[i]
            node.adh=nbs[i,1]
        nodes2=sorted(nodes,key=lambda nodes:nodes.adh)
        p=int(perc*len(nodes2))
        outside_nodes=nodes2[0:p]
        ## may need to be nodes2 - 1
        inside_nodes=nodes2[len(nodes2)-p:len(nodes2)]
            
        return outside_nodes, inside_nodes

    ##outsidein
##    sname="outin"
##    perc=0.2+0.45*r.random()
##    outside_nodes,inside_nodes=closest_nbs(net,80,perc) 
##    for i in range(n):
##        agent=nodes[i]
##        if agent in outside_nodes:
##            agent.state="D"
##        else:
##            agent.state="U"

    
##    ## inside out
##    sname="inout"
##    perc=0.4+0.45*r.random()
##    print perc
##    outside_nodes,inside_nodes=closest_nbs(net,80,perc) 
##    for i in range(n):
##        agent=nodes[i]
##        if agent in outside_nodes:
##            agent.state="U"
##        else:
##            agent.state="D"
 

    s=ZeroPad(50,ts)          
    name=spath+sname+s+".gpickle"
    print name
    nx.write_gpickle(net,name)
    ts=ts+1
