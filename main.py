from mpl_toolkits import mplot3d
#%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt

# 3d plotting on pyplot: https://jakevdp.github.io/PythonDataScienceHandbook/04.12-three-dimensional-plotting.html

fig = plt.figure()
ax = plt.axes(projection='3d')

"""
# Data for a three-dimensional line
zline = np.linspace(0, 15, 1000)
xline = np.sin(zline)
yline = np.cos(zline)
ax.plot3D(xline, yline, zline, 'gray')

# Data for three-dimensional scattered points
zdata = 15 * np.random.random(100)
xdata = np.sin(zdata) + 0.1 * np.random.randn(100)
ydata = np.cos(zdata) + 0.1 * np.random.randn(100)
ax.scatter3D(xdata, ydata, zdata, c=zdata, cmap='Greens');

fig.savefig("model1.jpg")
"""
# Define wireframe boundary
boundary = np.array([[0,1,1,1,1,0,0,0,0],
                     [0,0,0,1,1,1,1,0,0],
                     [0,0,1,1,0,0,1,1,0]])
xbound = boundary[0]
ybound = boundary[1]
zbound = boundary[2]
# Plot wireframe
ax.plot3D(xbound, ybound, zbound, 'gray')
print("Plotted: wire frame")

# Expand boundary

# Number of points along each length of boundary
n_side = 10

boundary0 = np.array([[0],[0],[0]])
for i in range(len(boundary[0])-1):
  #print(boundary[:,i][np.newaxis,:].T)
  #print([np.linspace(boundary[0,i], boundary[0,i+1], n_side)])
  #print([np.linspace(boundary[r:r+1,i], boundary[r:r+1,i+1], n_side) for r in range(3)])
  linB = np.concatenate([np.linspace(boundary[r:r+1,i], boundary[r:r+1,i+1], n_side) for r in range(3)],axis=1).T
  
  #print(linB,end = "\n\n")
  
  boundary0 = np.concatenate((boundary0, linB), axis = 1)
  
  #boundary0 = np.concatenate((boundary0, boundary[:,i][np.newaxis,:].T), axis = 1)
# boundary[:][i]

# print(boundary0)
ax.scatter3D(boundary0[0], boundary0[1], boundary0[2], c=[0 if i==0 else 1 for i in range(len(boundary0[0]))], cmap='Greens');
print("Finished establishing soap bounds")

# Define initial soap bubble
soap0 = np.array([[0],[0],[0]])
for i in range(1*n_side,4*n_side):
  linB = np.concatenate([np.linspace(boundary0[r:r+1,i+1], boundary0[r:r+1,9*n_side - i], n_side) for r in range(3)],axis=1).T
  
  #print(linB,end = "\n\n")
  
  soap0 = np.concatenate((soap0, linB), axis = 1)

ax.scatter3D(soap0[0], soap0[1], soap0[2], c=[0 if i==0 else 1 for i in range(len(soap0[0]))], cmap='Reds');
print("Finished generating initial soap bubble")

fig.savefig("model2.jpg")