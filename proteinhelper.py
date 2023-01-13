from vectors import *
import matplotlib.pyplot as plt
from kearsley import Kearsley

k = Kearsley()

class Point:
    def __init__(self,coord,*args): #coord is a 1D array of 3 numbers
        self.coord = Vec(coord)
        if len(args) != 0:
            self.label = args[0]
    def translate(self,vector): #vector by which to translate (1D array)
        self.coord.add(vector)
    def rotate(self,axis,angle): #axis is a vector around which to rotate angle is in radians + when anticlockwise looking down in the direction opp to the axis
        axis_cap = axis.cap()
        A_vec = Vec(self.coord.vec - (self.coord.dot(axis_cap))*axis_cap.vec)
        A1_vec = Vec((A_vec.mod()*mt.cos(angle))*A_vec.cap().vec + (A_vec.mod()*mt.sin(angle))*(A_vec.cross(axis_cap).vec))        
        self.coord = Vec(A1_vec.vec + self.coord.dot(axis_cap)*axis_cap.vec)

class Protein:
    def __init__(self,points): #points is a list of Point objects
        self.points = points
    
    def translate(self,vector):
        for p in self.points:
            p.translate(vector)
    
    def rotate(self,axis,angle):
        for p in self.points:
            p.rotate(axis,angle)
        
    def display(self,col):
        fig = plt.figure()
        ax = plt.axes(projection ='3d')
        x,y,z = [], [], []
        for p in self.points:
            x.append(p.coord.vec[0])
            y.append(p.coord.vec[1])
            z.append(p.coord.vec[2])
        xs = np.array(x)
        ys = np.array(y)
        zs = np.array(z)
        ax.scatter(x,y,z,c=col)
        plt.show()

class Pairings:
    def __init__(self,A,B): #accepts structures A and B
        self.A = A
        self.B = B