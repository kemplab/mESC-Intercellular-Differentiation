import math as math
import random as rand
import numpy as np
from scipy import ndimage

def RandomPointOnSphere(dim):
    """ Computes a random point on a sphere
        Returns - a point on a unit sphere [x,y,z] at the origin
    """
    if dim==3:
        u = rand.random()
        theta = rand.random()*2*math.pi
        x = math.sqrt(1-(u*u))*math.cos(theta)
        y = math.sqrt(1 - (u*u))*math.sin(theta)
        z = u
        return np.array((x,y,z))
    else:
        theta = rand.random()*2*math.pi
        x = math.cos(theta)
        y = math.sin(theta)
        
        return np.array((x,y,0))


def AddVec(v1, v2):
    """ Adds two vectors that are in the form [x,y,z]
        Returns - a new vector [x,y,z] as a numpy array
    """
    n=len(v1);
    temp=np.array(v1)
    for i in range(0,n):
        temp[i]+=float(v2[i])
    return temp

def SubtractVec(v1, v2):
    """ Subtracts vector [x,y,z] v2 from vector v1
        Returns - a new vector [x,y,z] as a numpy array
    """
    n=len(v1);
    temp=np.array(v1)
    for i in range(0,n):
        temp[i]-=float(v2[i])
    return temp
    

def ScaleVec(v1, s):
    """ Scales a vector f*[x,y,z] = [fx, fy, fz]
        Returns - a new scaled vector [x,y,z] as a numpy array
    """
    n=len(v1);
    temp=np.array(v1)
    for i in range(0,n):
        temp[i]=temp[i]*s
    return temp

def Mag(v1):
    """ Computes the magnitude of a vector
        Returns - a float representing the vector magnitude
    """
    n=len(v1);
    temp=0.
    for i in range(0,n):
        temp+=(v1[i]*v1[i])
    return math.sqrt(temp)

def NormVec(v1):
    """ Computes a normalized version of the vector v1
        Returns - a normalizerd vector [x,y,z] as a numpy array
    """
    
    mag = Mag(v1)
    temp=np.array(v1)
    if(mag == 0):
        return temp*0
    return temp/mag

def Distance(p1, p2):
    """ Computes the distance between two points [x,y,z]
        Returns - a float representing the distance
    """
    n=len(p1);
    
    temp=0.
    for i in range(0,n):
        temp+=(p1[i]-p2[i])*(p1[i]-p2[i])
    temp=math.sqrt(temp)
    return temp
    


