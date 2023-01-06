import numpy as np
import math as mt

class Vec:
    def __init__(self,L) -> None:
        self.vec = np.array(L)

    def add(self,V):
        self.vec = self.vec + V.vec

    def dot(self,V): #dot product
        return sum(self.vec*V.vec)

    def mod(self): #returns magnitude of vector
        return mt.sqrt(self.dot(self))
    
    def cap(self):
        return Vec(self.vec/self.mod())
        
    def cross(self,V): #cross product
        [a1,a2,a3] = self.vec
        [b1,b2,b3] = V.vec
        C = Vec([(a2*b3)-(b2*a3),(b1*a3)-(a1*b3),(a1*b2)-(b1*a2)])
        return C

i = Vec([1,0,0])
j = Vec([0,1,0])
k = Vec([0,0,1])