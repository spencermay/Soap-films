from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
import timeit
#from scipy.spatial import Delaunay, ConvexHull
#Delaunay no longer needed?
#import matplotlib.tri as mtri
#import open3d as o3d
#import pygalmesh
#import pygmsh


# 3d plotting on pyplot: https://jakevdp.github.io/PythonDataScienceHandbook/04.12-three-dimensional-plotting.html

#Deprecated version:
#fig = plt.figure()
#ax = fig.gca(projection='3d')

fig = plt.figure()
ax = plt.axes(projection='3d')

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


#_, soap_idx = np.unique(soap0, return_index=True)
#soap0 = soap0[np.sort(soap_idx)]

#soap0 = np.unique(soap0, axis=0)



#####################################################




# Plot initial soap bubble
# ax.scatter3D(soap0.T[0], soap0.T[1], soap0.T[2], c=[0 if i==0 else 1 for i in range(len(soap0))], cmap='Reds');
print("Finished generating initial soap bubble")

# Number of neighbors that influence a given point's behavior
attentionN = 8
# Minimum distance beyond which to stop paying attention to neighbors
dMin = 0.05
# Multiplier to bias towards points along wire
mult = 2
# Step-length
step = 0.01
# Number of iterations to perform
REPS = 4
# Display it as a surface or collection of points?
surf = True
# Number of points remains constant?
constN = True#False
# Used refining technique?
sRefined = False#True

def find_nearest(bubblepts0):
  dirs = []
  for pt in bubblepts0:
    #print(pt)
    # Find each point's closest neighbors:
    # For n smallest items in a numpy array: https://stackoverflow.com/questions/33623184/fastest-method-of-getting-k-smallest-numbers-in-unsorted-list-of-size-n-in-pytho
    #if pt in boundary0.tolist():
    if any(np.equal(boundary0,pt).all(1)) and constN:
      dirs.append([0,0,0])
      continue

    dists = [sum([(p1[i] - pt[i])**2 if (p1[i] - pt[i])**2 > dMin else 100 for i in range(3)]) for p1 in bubblepts0]
    #print(f"dists: {dists}")
    neighborhood = np.partition(dists, attentionN)[:attentionN]
    #print(neighborhood)
    vec_sum = [0,0,0]
    # neighbors = []
    for k in range(len(bubblepts0)):
      if dists[k] in neighborhood:
        neighbor = bubblepts0[k]
        # neighbors.append(neighbor)
        mag = sqrt(dists[k]) # sqrt(sum([(neighbor[c]-pt[c])**2 for c in range(3)]))
        
        for coord in range(3):
          if any(np.equal(boundary0,neighbor).all(1)):
            vec_sum[coord] += mult * (neighbor[coord] - pt[coord]) / mag * step
          else:
            vec_sum[coord] += (neighbor[coord] - pt[coord]) / mag * step
    
    dirs.append(vec_sum)

  if constN:
    return(dirs)
  else:
    return(dirs + len(boundary0)*[[0,0,0]])


def find_tris(bubblepts):
  return np.array([[[n,n+1,n+11], [n,n+10,n+11]] for n in range(len(bubblepts)-11) if n%10 != 0]).reshape(-1,3)



def area(tri):
  return(0.5 * np.linalg.norm(np.cross(tri[0]-tri[1],tri[2]-tri[1])))

def refine(bubblepts, tris): 
  # Find the triangulation of the surface -- use find_tris()

  
  """
  # Article on surface reconstruction with python: https://towardsdatascience.com/5-step-guide-to-generate-3d-meshes-from-point-clouds-with-python-36bad397d8ba
  print("Beginning")
  pcd = o3d.geometry.PointCloud()
  pcd.points = o3d.utility.Vector3dVector(bubblepts)
  pcd.normals = o3d.geometry.estimate_normals(pcd,search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=0.1, max_nn=30))

  # Calculate mesh
  poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd, depth=8, width=0, scale=1.1, linear_fit=False)[0]

  # Crop mesh to size of original point cloud
  bbox = pcd.get_axis_aligned_bounding_box()
  p_mesh_crop = poisson_mesh.crop(bbox)
  print("p_mesh:", p_mesh_crop)
  """
  #tris = Delaunay(bubblepts)
  ##### tris = ConvexHull(bubblepts)
  #tris = mtri.Triangulation(bubblepts)
  # Create an empty list in which to store the adjustments that need to be made to the mesh points
  adjustments = [[] for j in range(len(bubblepts))]
  # For each triangle, calculate the adjustments to its vertices based solely on that triangle
  
  """
  ax.plot_trisurf(bubblepts[:,0], bubblepts[:,1], bubblepts[:,2])#, tris.simplices)
  fig.savefig("model3_2.jpg")

  """
  for tri in tris:
    #tris.triangles:#tris.simplices[:10]:
    # Convert the indices stored in tris.simplices to actual points
    triangle = [bubblepts[x] for x in tri]
    # Calculate areas using cross product
    # Triangle needs to be a list of numpy arrays
    triArea = area(triangle)
    
    # Adjust in direction of that area for each vertex (sum two vectors and multiply by area swept out between them)
    
    # For each of the indices stored in tri (for each vertex)
    for i,p in enumerate(tri):
      # If the point is on the boundary and it's set to not increase the number of points:
      print(triangle[i])
      if any(np.equal(boundary0,triangle[i]).all(1)) and constN:
        adjustments[p] = np.array([0], ndmin=1)
        #print(p)
        continue

      #print("T", np.array(triArea * sum( \
      #  [triangle[(g + i + 1) % 3] for g in range(2)]), ndmin=1))
      # Add to the list of adjustments
      adjustments[p].append(np.array(triArea * sum( \
        [triangle[(g + i + 1) % 3] for g in range(2)]), ndmin=1))
  #"#""
  print(adjustments[0:11])
  #for k in adjustments:
  #  print(k)
  #print(np.array([len(k) for k in adjustments]))
  #print([np.ndarray.size(n) for n in adjustments])
  #"#""
  
  #print([np.ndarray.size(n) for n in adjustments])
  #print(len([np.sum(n) for n in adjustments]))
  #print(([n for n in range(len(adjustments)) if len(adjustments[n])==0]))
  #print(len([(step * (np.sum(n)) / len(n)) for n in adjustments]))

  #print(np.array([(step * (np.sum(n)) / len(n)) for n in adjustments]))
  print(np.array([len(n) for n in adjustments]))
  return(np.array([(step * (np.sum(n)) / len(n)) for n in adjustments]))
  """"""
  #return([step * (np.sum(n)) / np.ndarray.size(n) for n in adjustments])
  #plt.triplot(bubblepts[:,0], bubblepts[:,1], bubblepts[:,2], tri.simplices)

def bubble_evolve(bubblepts0, reps):
  soap = bubblepts0
  
  if sRefined == True:
    tris = [[n,n+1,n+10] for n in range(len(bubblepts0))]
    plt.triplot(bubblepts[:,0], bubblepts[:,1], bubblepts[:,2], tris)

  
  for rep in range(reps):
    #print(np.array(refine(soap))[:40])

    
    if sRefined == True:
      

      
      if constN:
        soap = np.add(soap, refine(soap))
      else:
        soap = np.add(np.concatenate([soap, boundary0], axis=0), refine(soap))

    else:
      if constN:
        soap = np.add(soap, find_nearest(soap))
      else:
        soap = np.add(np.concatenate([soap, boundary0], axis=0), find_nearest(soap))
    print(f"Iteration: {rep + 1}")
  return soap



#"""
############### film = bubble_evolve(soap0, REPS)
#print(find_nearest(soap0))

"""
if surf:
  ax.plot_trisurf(film.T[0], film.T[1], film.T[2], cmap='viridis', edgecolor='none');
else:
  ax.scatter3D(film.T[0], film.T[1], film.T[2], c=[0 if i==0 else 1 for i in range(len(film))], cmap='Reds');
# Save resulting figure as an image
fig.savefig("model3_1.jpg")
#"#""
"""





def find_nearest0():
  return(find_nearest(soap0))

def find_nearest_time():
  SETUP_CODE = '''
from __main__ import find_nearest0
import numpy as np'''

  TEST_CODE = '''
find_nearest0()'''

  print (timeit.timeit(setup = SETUP_CODE,
                     stmt = TEST_CODE,
                     number = 4))
#find_nearest_time()


def bubble_evolve0():
  return(bubble_evolve(soap0, REPS))
  
def bubble_evolve_time():
  SETUP_CODE = '''
from __main__ import bubble_evolve0
import numpy as np'''

  TEST_CODE = '''
bubble_evolve0()'''

  print (timeit.timeit(setup = SETUP_CODE,
                     stmt = TEST_CODE,
                     number = 1))
# bubble_evolve_time()

# bubble_evolve with 4 iterations takes 34.5 sec, find_nearest * 4 takes 20 sec

from mpl_toolkits.mplot3d import Axes3D

bubblepts=soap0

#tris = np.array([[[n,n+1,n+10], [n,n+9,n+10]] for n in range(len(bubblepts)-10) if n%10 == 9][:40]).reshape(-1,3)

#print(tris)
ax.plot_trisurf(bubblepts[:,0], bubblepts[:,1], bubblepts[:,2], triangles = tris)
fig.savefig("begin0.jpg")

#Commented out np.unique line