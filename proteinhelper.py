from vectors import *
import matplotlib.pyplot as plt
from kearsley import Kearsley
from scipy.spatial.transform import Rotation as R

k = Kearsley()

# class Point:
#     def __init__(self,coord,*args): #coord is a 1D array of 3 numbers
#         self.coord = Vec(coord)
#         if len(args) != 0:
#             self.label = args[0]
#     def translate(self,vector): #vector by which to translate (1D array)
#         self.coord.add(vector)
#     def rotate(self,axis,angle): #axis is a vector around which to rotate angle is in radians + when anticlockwise looking down in the direction opp to the axis
#         axis_cap = axis.cap()
#         A_vec = Vec(self.coord.vec - (self.coord.dot(axis_cap))*axis_cap.vec)
#         A1_vec = Vec((A_vec.mod()*mt.cos(angle))*A_vec.cap().vec + (A_vec.mod()*mt.sin(angle))*(A_vec.cross(axis_cap).vec))        
#         self.coord = Vec(A1_vec.vec + self.coord.dot(axis_cap)*axis_cap.vec)

class Structure:
    def __init__(self,points): #points is a nx3 array representing pointsin the structure
        self.points = points
    
    def translate(self,vector): #vector is a 1x3 array representing translation vector
        n = self.points.shape[0]
        L = [vector]*n
        self.points += np.array(L)

    
    def rotate(self,rotation_vec): #rotation_vec is a vector of axis of ratation whith magnitude = angle or rotation in radians
        r = R.from_rotvec(rotation_vec)
        self.points = r.apply(self.points) 

    def display(self,col):
        fig = plt.figure()
        ax = plt.axes(projection ='3d')
        x = self.points[:,0]
        y = self.points[:,1]
        z = self.points[:,2]
        ax.scatter(x,y,z,c=col)
        plt.show()


class Pairing:
    def __init__(self,permutation,A,B): #permutation is a list of size N containing integers 1 to N, such that L[i] = j means point i of A is paired to point j of B 
        self.perm = permutation
        self.A = A
        self.B = B
    
    def kearsley(self):
        V = []
        for i in range(self.B.points.shape[0]):
            V.append(self.B.points[self.perm[i]])
        V_arr = np.array(V)
        U = self.A.points
        k.fit(U,V_arr)
        V_trans = k.transform(V_arr)
        self.B.points = V_trans

    def display(self,col1,col2):
        fig = plt.figure()
        ax = plt.axes(projection ='3d')
        x = self.A.points[:,0]
        y = self.A.points[:,1]
        z = self.A.points[:,2]
        ax.scatter(x,y,z,c=col1)
        x = self.B.points[:,0]
        y = self.B.points[:,1]
        z = self.B.points[:,2]
        ax.scatter(x,y,z,c=col2)
        plt.show()



    

    
    


