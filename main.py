from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

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
                     [0,0,1,1,0,0,1,1,0]]).T
xbound = boundary.T[0]
ybound = boundary.T[1]
zbound = boundary.T[2]

"""
Print coordinates along the boundary:
for ptx, pty, ptz in boundary.T:
  print(f"({ptx},{pty},{ptz})")
"""

# Plot wireframe
ax.plot3D(xbound, ybound, zbound, 'gray')
print("Plotted: wire frame")

# Number of points along each length of boundary
n_side = 10

# Expand boundary
boundary0 = np.array([[0,0,0]])

for i in range(len(boundary)-1):
  linB = np.concatenate([np.linspace(boundary[i,r:r+1], boundary[i+1,r:r+1], n_side) for r in range(3)], axis=1)
  boundary0 = np.append(boundary0, linB, axis=0)
#print(boundary0)

# Plot soap bounds
ax.scatter3D(boundary0.T[0], boundary0.T[1], boundary0.T[2], c=[0 if i==0 else 1 for i in range(len(boundary0))], cmap='Greens');
print("Finished establishing soap bounds")

######### CAN IT BE 3 INSTEAD????

# Define initial soap bubble
soap0 = np.array([[0,0,0]])
for i in range(1*n_side,4*n_side):
  linB = np.concatenate([np.linspace(boundary0[i+1,r:r+1], boundary0[9*n_side - i, r:r+1], n_side) for r in range(3)], axis=1)
  soap0 = np.append(soap0, linB, axis=0) #, axis=1)

# Plot initial soap bubble
# ax.scatter3D(soap0.T[0], soap0.T[1], soap0.T[2], c=[0 if i==0 else 1 for i in range(len(soap0))], cmap='Reds');
print("Finished generating initial soap bubble")

# Number of neighbors that influence a given point's behavior
attentionN = 8
# Minimum distance beyond which to stop paying attention to neighbors
dMin = 0.05
# Multiplier to bias towards points along wire
mult = 40
# Step-length
step = 0.0001
# Number of iterations to perform
REPS = 1
# Display it as a surface or collection of points?
surf = False#True

def find_nearest(bubblepts0):
  dirs = []
  for pt in bubblepts0:
    #print(pt)
    # Find each point's closest neighbors:
    # For n smallest items in a numpy array: https://stackoverflow.com/questions/33623184/fastest-method-of-getting-k-smallest-numbers-in-unsorted-list-of-size-n-in-pytho
    #if pt in boundary0.tolist():
    if any(np.equal(boundary0,pt).all(1)):
      dirs.append([0,0,0])
      continue

    dists = [sum([(p1[i] - pt[i])**2 if (p1[i] - pt[i])**2 > dMin else 100 for i in range(3)]) for p1 in bubblepts0]
    #print(f"dists: {dists}")
    neighborhood = np.partition(dists, attentionN)[:attentionN]
    #print(neighborhood)
    vec_sum = [0,0,0]
    neighbors = []
    for k in range(len(bubblepts0)):
      if dists[k] in neighborhood:
        neighbor = bubblepts0[k]
        neighbors.append(neighbor)
        mag = sqrt(sum([(neighbor[c]-pt[c])**2 for c in range(3)]))
        #### sqrt dists[k] #####################
        
        for coord in range(3):
          if neighbor in boundary0:
            vec_sum[coord] += mult * (neighbor[coord] - pt[coord]) / mag * step
          else:
            vec_sum[coord] += (neighbor[coord] - pt[coord]) / mag * step
    
    dirs.append(vec_sum)
  return(dirs)

def bubble_evolve(bubblepts0, reps):
  soap = bubblepts0
  for rep in range(reps):
    soap = np.add(soap, find_nearest(soap0))
    print(f"Iteration: {rep + 1}")
  return soap

film = bubble_evolve(soap0, REPS)
#print(find_nearest(soap0))

if surf:
  ax.plot_trisurf(film.T[0], film.T[1], film.T[2], cmap='viridis', edgecolor='none');
else:
  ax.scatter3D(film.T[0], film.T[1], film.T[2], c=[0 if i==0 else 1 for i in range(len(film))], cmap='Reds');
# Save resulting figure as an image
fig.savefig("model3_0.jpg")