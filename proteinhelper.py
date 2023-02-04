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

def dist(a: np.ndarray,b: np.ndarray):
    S = 0
    for i in range(len(a)):
        S += (a[i]-b[i])**2
    
    return mt.sqrt(S)

def norm(v: np.ndarray):
    return dist(np.array([0,0,0]),v)

class Structure:
    def __init__(self,*args,**kargs): #points is a nx3 array representing points in the structure
        if kargs == {}:
            kargs['random'] = False
        if kargs['random']:
            L = []
            vec = []
            for i in range(kargs['size']):
                L.append(np.random.randint(-5,5,size=3))
                V3 = [np.random.randint(-5,5,size=3),np.random.randint(-5,5,size=3),np.random.randint(-5,5,size=3)]
                vec.append([V3[0]/norm(V3[0]),V3[1]/norm(V3[1]),V3[2]/norm(V3[2])])
            points = np.array(L)
            vectors = np.array(vec)
        else:
            points = args[0]
            if len(args) >= 2:
                vectors = args[1]
            else:
                vectors = np.array([[[0,0,0],[0,0,0],[0,0,0]]])
        self.points = points
        self.vectors = vectors #vectors is a nx3x3 3D array where vectors[i][j] is the jth vector of the ith point [[[1,1,3],[1,4,2],[2,3,2]],]
    
    def translate(self,vector: np.ndarray): #vector is a 1x3 array representing translation vector
        n = self.points.shape[0]
        L = [vector]*n
        self.points += np.array(L)

    
    def rotate(self,rotation_object,**kargs): #rotation_vec is a vector of axis of ratation whith magnitude = angle or rotation in radians
        r = rotation_object
        if kargs == {}:
            kargs['only_vectors'] = False
        if kargs["only_vectors"]:
            n = self.vectors.shape[0]
            for i in range(n):
                self.vectors[i] = r.apply(self.vectors[i])
        else:
            self.points = r.apply(self.points)
            n = self.vectors.shape[0]
            for i in range(n):
                self.vectors[i] = r.apply(self.vectors[i])

    def display(self,col):
        fig = plt.figure()
        ax = plt.axes(projection ='3d')
        x = self.points[:,0]
        y = self.points[:,1]
        z = self.points[:,2]
        ax.scatter(x,y,z,c=col)
        plt.show()


class Pairing:
    def __init__(self,permutation,A,B): #permutation is a list of size N containing integers 0 to N-1, such that L[i] = j means point i of A is paired to point j of B 
        self.perm = permutation
        self.A = A
        self.B = B
    
    def rmsd(self):
        S = 0
        n = self.A.points.shape[0]
        for i in range(n):
            S += (self.A.points[i] - self.B.points[self.perm[i]])**2
        S = sum(S)
        return mt.sqrt(S/n)
        
    def kearsley(self):
        V = []
        for i in range(self.B.points.shape[0]):
            V.append(self.B.points[self.perm[i]])
        V_arr = np.array(V)
        U = self.A.points
        rmsd = k.fit(U,V_arr)
        V_trans = k.transform(V_arr)
        self.B.points = V_trans
        self.B.rotate(k.rot,only_vectors = True)
        return rmsd, k.trans, k.rot

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

                


def make_neibourhood(a: int,A: Structure,m:int):
    D = {}
    for i in range(A.points.shape[0]):
        D[i] = dist(A.points[a],A.points[i])
    m_closest = sorted(D.items(),key= lambda x:x[1])[1:m+1]
    norbors = [x[0] for x in m_closest]
    neibours = [a]
    neibours.extend(norbors)
    
    return Structure(A.points[neibours],A.vectors[neibours]), neibours





    
    


