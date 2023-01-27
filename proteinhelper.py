from vectors import *
import matplotlib.pyplot as plt
from kearsley import Kearsley
from scipy.spatial.transform import Rotation as R
from scipy.optimize import linear_sum_assignment

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

def dist(a,b):
    S = 0
    for i in range(len(a)):
        S += (a[i]-b[i])**2
    
    return mt.sqrt(S)


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
        rmsd = k.fit(U,V_arr)
        V_trans = k.transform(V_arr)
        self.B.points = V_trans
        return rmsd

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

    def hungarian(self):
        cost = np.zeros((len(self.A.points),len(self.B.points)))
        i = -1
        for a in self.A.points:
            i += 1
            j = -1
            for b in self.B.points:
                j += 1
                cost[i,j] = dist(a,b)
        
        row_ind,col_ind = linear_sum_assignment(cost)
        self.perm = col_ind
        return col_ind

                


    

    
    


