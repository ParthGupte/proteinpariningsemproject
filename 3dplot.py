import numpy as np
import matplotlib.pyplot as plt


fig = plt.figure()
ax = plt.axes(projection ='3d')
x = np.array([0,1,2,3,4,5])
y = np.array([0,1,2,3,4,5])
z = np.array([0,1,2,3,4,5])

ax.scatter(x,y,z,'red')
plt.show()