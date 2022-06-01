from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from math import sqrt, ceil, floor
import timeit
# from scipy.spatial import Delaunay, ConvexHull
# Delaunay no longer needed?
# import matplotlib.tri as mtri
# import open3d as o3d
# import pygalmesh
# import pygmsh


# 3d plotting on pyplot: https://jakevdp.github.io/PythonDataScienceHandbook/04.12-three-dimensional-plotting.html

# Deprecated version:
# fig = plt.figure()
# ax = fig.gca(projection='3d')


# Number of points along each length of boundary
n_side = 10
# Number of neighbors that influence a given point's behavior
attentionN = 8
# Minimum distance beyond which to stop paying attention to neighbors
dMin = 0.1
# Multiplier to bias towards points along wire
mult = 1#2
# Step-length. # 0.01 #5 for n=5; 70 for n=10
step = 0.005  # 0.015
# Number of iterations to perform # 4
REPS = 15
# Display it as a surface or collection of points?
surf = True
# Number of points remains constant?
constN = True  # False
# Used refining technique?
sRefined = False  # True
# Border precision -- number of decimal places to check to figure out if a point is on the border:
bPrecision = 10
# Should surface area be tracked?
surfArea = True
# List of steps:
stepsList = [["find_nearest"]] + [["refine", 0.2]] + [["add", 30]] #+ [["refine", 0.3]]  # 5 * [["refine",1.7]] + [["add",3]] #+ 2 * [["refine",1.3]]
# Weightings on final refining steps
Z = [0.000000001, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1] + 12 * [0.15]
# Power to which to raise the norm of the vector in the refine function
refpow = 1.0 # 1.6

# fig = plt.figure()
# ax = plt.axes(projection='3d')
# Calculate the number of plots to make:
for rows in range(floor(sqrt(REPS)),1,-1):
  if REPS % rows == 0:
    break
fig = plt.figure()
axs = [fig.add_subplot(rows, REPS//rows, rep + 1, projection='3d') for rep in range(REPS)]

# Define wireframe boundary
boundary = np.array([[0, 1, 1, 1, 1, 0, 0, 0, 0],
                     [0, 0, 0, 1, 1, 1, 1, 0, 0],
                     [0, 0, 1, 1, 0, 0, 1, 1, 0]]).T
xbound = boundary.T[0]
ybound = boundary.T[1]
zbound = boundary.T[2]

"""
Print coordinates along the boundary:
for ptx, pty, ptz in boundary.T:
  print(f"({ptx},{pty},{ptz})")
"""

# Plot wireframe
for ax in axs:
  ax.plot3D(xbound, ybound, zbound, 'gray')
print("Plotted: wire frame")


# Expand boundary
boundary0 = np.array([[0,0,0]])

for i in range(len(boundary)-1):
  linB = np.concatenate([np.linspace(boundary[i,r:r+1], boundary[i+1,r:r+1], n_side) for r in range(3)], axis=1)
  boundary0 = np.append(boundary0, linB, axis=0)
#print(boundary0)

# Plot soap bounds
for ax in axs:
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



def find_nearest(bubblepts0):
  # Create an empty list that will be filled with the directions in which to adjust each point
  dirs = []
  # Go through each point one at a time
  for pt in bubblepts0:
    # Find that point's closest neighbors:
    # For n smallest items in a numpy array: https://stackoverflow.com/questions/33623184/fastest-method-of-getting-k-smallest-numbers-in-unsorted-list-of-size-n-in-pytho
    # Check if the point, "pt", is on the boundary
    #if any(np.equal(boundary0,pt).all(1)) and constN:
    if any(np.equal(np.around(boundary0,bPrecision),                np.around(pt,bPrecision)).all(1)):
      # If the point is on the boundary it shouldn't be adjusted
      # (the adjustment for that point is zero)
      dirs.append([0,0,0])
      # Skip the rest of the adjustment calculation
      continue

    # Find the squares of the distances from "pt" to the other points
    dists = [sum([(p1[i] - pt[i])**2 if (p1[i] - pt[i])**2 > dMin else 100 for i in range(3)]) for p1 in bubblepts0]
    # Find the N nearest neighbors to point "pt"
    neighborhood = np.partition(dists, attentionN)[:attentionN]
    # Start the vector sum off at zero:
    vec_sum = [0,0,0]
    # For each of the other points:
    for k in range(len(bubblepts0)):
      # Check if the other point is in the point neighborhood of "pt":
      if dists[k] in neighborhood:
        # If it is in the neighborhood, use it in the vector sum
        neighbor = bubblepts0[k]

        # Calculate magnitude
        mag = sqrt(dists[k])

        # Iterate through the coordinates (x, y, and z) one by one and find the coordinates of the vector sum
        for coord in range(3):
          if any(np.equal(boundary0,neighbor).all(1)):
            vec_sum[coord] += mult * (neighbor[coord] - pt[coord]) / mag * step
          else:
            vec_sum[coord] += (neighbor[coord] - pt[coord]) / mag * step

    dirs.append(vec_sum)

  # If the program is set to keep the same number of points, return the adjustments
  if constN:
    return(dirs)
  # Otherwise, return the adjustments, with the appropriate number of zeros added to the end (representing the zero adjustments made to the new points)
  else:
    return(dirs + len(boundary0)*[[0,0,0]])

# Find the triangulation for the initial set of points
def find_tris(bubblepts, tris0=[]):
  if not tris0:
    return np.array([[[n, n + 1, n + n_side + 1], [n, n + n_side, n + n_side + 1]] for n in range(len(bubblepts) - n_side - 1) if n % n_side != 0]).reshape(-1,3)
  else:
    pass


def area(tri):
  return(0.5 * np.linalg.norm(np.cross(tri[0]-tri[1],tri[2]-tri[1])))

def refine(bubblepts, tris, stp = step):
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
  if surfArea:
    surfA = 0
  for tri in tris: #[0:1]
    #tris.triangles:#tris.simplices[:10]:
    # Convert the indices stored in tris.simplices to actual points
    triangle = [bubblepts[x] for x in tri]
    ###print("tri",tri)
    ###print("triangle",triangle)
    # Calculate areas using cross product
    # Triangle needs to be a list of numpy arrays
    triArea = area(triangle)
    if surfArea:
      surfA += triArea
    ###print("triArea", triArea)

    # Adjust in direction of that area for each vertex (sum two vectors and multiply by area swept out between them)

    # For each of the indices stored in tri (for each vertex)
    for i,p in enumerate(tri):
      # If the point is on the boundary:
      #no longer needed for refine() --> "and it's set to not increase the number of points"
      #print(triangle[i])
      #if any(np.equal(boundary0,triangle[i]).all(1)):
      if any(np.equal(np.around(boundary0,bPrecision),                np.around(triangle[i],bPrecision)).all(1)):
        #adjustments[p] = np.array([0], ndmin=1)
        adjustments[p] = [[0,0,0]]
        ###print("p",p)
        continue

      #print("T", np.array(triArea * sum( \
      #  [triangle[(g + i + 1) % 3] for g in range(2)]), ndmin=1))
      # Add to the list of adjustments
      #adjustments[p].append(np.array(triArea * sum( \
      #  [triangle[(g + i + 1) % 3] for g in range(2)]), ndmin=1))
      #print(triArea * sum( \
      #  [triangle[(g + i + 1) % 3] for g in range(2)]))

      # Calculate vector sum for adjustment:

      #adjustments[p].append(list(triArea * (sum( \
      #  [triangle[(i + g + 1) % 3] for g in range(2)]) - triangle[i])))
      adjustments[p].append(list(triArea * (triangle[(i+1)%3] + triangle[(i+2)%3] - 2*triangle[i])))
      ###print(adjustments[p])

  ###print("adjustments", adjustments[0:15])
  #if not adjustments[0]:
  #  adjustments[0]=[[0,0,0]]

  # If for some reason one of the points doesn't appear in the mesh, don't adjust it
  adjustments = [[[0,0,0]] if y==[] else y for y in adjustments]
  #"#""

  #for k in adjustments:
  # print(k)
  # print(np.array([len(k) for k in adjustments]))
  # print([np.ndarray.size(n) for n in adjustments])
  #"#""

  # print([np.ndarray.size(n) for n in adjustments])
  # print(len([np.sum(n) for n in adjustments]))
  # print(([n for n in range(len(adjustments)) if len(adjustments[n])==0]))
  # print(len([(step * (np.sum(n)) / len(n)) for n in adjustments]))

  # print(np.array([(step * (np.sum(n)) / len(n)) for n in adjustments]))
  # print(np.array([len(n) for n in adjustments]))
  # print([(step * (np.sum(n)) / len(n)) for n in adjustments])
  # print("adj2", [[j[k] for k in range(3)] for n in adjustments for j in n][:10])
  # print("adj3", [[sum([j[k] for j in n]) for k in range(3)] for n in adjustments][:10])
  # print([step * ([sum(n[:,k]) for k in range(3)]) / len(n) for n in adjustments])
  adj = [[sum([j[k] for j in n]) for k in range(3)] for n in adjustments]
  norms = [np.linalg.norm(u) for u in adj]
  # print(norms)
  ### print("adj final", np.array([[1 * step * sum([j[k] for j in n]) / len(n) for k in range(3)] for n in [[[0,0,0]] if y==[] else y for y in adjustments]][:15]))
  if surfArea:
    return([np.array([[1 * stp * sum([j[k] for j in n]) / len(n) / [o ** refpow if o else 1 for o in norms][N] for k in range(3)] for N,n in enumerate(adjustments)]),surfA])
    return([np.array([[1 * stp * sum([j[k] for j in n]) / len(n) / [o if o else 1 for o in norms][N] for k in range(3)] for N,n in enumerate(adjustments)]),surfA])

  return(np.array([[1 * stp * sum([j[k] for j in n]) / len(n) / [o if o else 1 for o in norms][N] for k in range(3)] for N,n in enumerate(adjustments)]))

  return(np.array([[1 * stp * sum([j[k] for j in n]) / len(n) for k in range(3)] for n in adjustments]))

  return(np.array([[1 * stp * sum([j[k] for j in n]) / len(n) / max(norms) for k in range(3)] for n in adjustments]))

  # return(np.array([(step * (np.sum(n)) / len(n)) for n in adjustments]))
  """"""

  # return([step * (np.sum(n)) / np.ndarray.size(n) for n in adjustments])
  # plt.triplot(bubblepts[:,0], bubblepts[:,1], bubblepts[:,2], tri.simplices)


###########################################################
def smooth(bubblepts, tris):
  # Convert the list of indices to a list of points (point-cloud)
  triangles = [bubblepts[tri] for tri in tris]
  for i, tripts in enumerate(triangles):
    pass



# Number of points to add:
attentionN1 = 20

def addpts(bubblepts, tris, N_to_add = 20):
  # New triangles in triangulation that need to be added
  new_tris = []
  # Old triangles in triangulation that need to be removed
  old_tri_indices = []
  # List of new points
  new_pts = []
  # Number of points
  npts = len(bubblepts)-1

  # Convert the indices stored in the triangulation to actual points
  triangles = [bubblepts[tri] for tri in tris]
  # Calculate the area of each triangle
  areas = [area(triangle) for triangle in triangles]
  # print("areas", areas[:10])
  # Find the N_to_add largest triangles to split up
  largest = max(areas)-np.partition([max(areas)-j for j in areas], N_to_add)[:N_to_add]
  # print("largest", largest)

  # Iterate through each triangle, finding the largest ones
  for i in range(len(tris)):
    # If the triangle is one of the largest ones:
    if areas[i] in largest:
      # Add the triangle to the list of triangles to be removed
      old_tri_indices.append(i)

      # Add a point at the center of the triangle to the list of points
      new_pts.append(np.ndarray.tolist(np.sum(triangles[i], axis=0)/3))
      # print("newpts",new_pts)
      # print("triangle", triangles[i])

      # Now there is one more point
      npts += 1
      # print(f"new_pts addition: {new_pts[-1]}")

      # Add triangles to the triangulation
      for pti in range(3):
        # Add a triangle excluding the pti'th point of the triangle, replacing this point with the new point.
        # print("Tris stuff", tris[i][:pti], tris[i][pti+1:], [npts])
        #new_tris.append(tris[i][:pti]+tris[i][pti+1:]+[npts])
        new_tris.append(np.array(np.ndarray.tolist(tris[i][:pti])+np.ndarray.tolist(tris[i][pti+1:])+[npts]))
        #new_tris.append(np.ndarray.tolist(tris[i][:pti])+np.ndarray.tolist(tris[i][pti+1:])+[npts])
  # print(f"bubblepts: {bubblepts[:10]}")
  # print(f"new_pts: {new_pts[:10]}")
  # print("tris removed",[tris[k] for k in range(len(bubblepts)) if k not in old_tri_indices][:10])
  # print("bubblepts out", np.array(np.ndarray.tolist(bubblepts) + new_pts)[:10])
  return(np.array(np.ndarray.tolist(bubblepts) + new_pts), [tris[k] for k in range(len(tris)) if k not in old_tri_indices]+new_tris)
  return(np.array(np.ndarray.tolist(bubblepts) + new_pts), [np.ndarray.tolist(tris[k]) for k in range(len(bubblepts)) if k not in old_tri_indices]+new_tris)



def bubble_evolve(bubblepts0, reps, steps):
  soap = bubblepts0

  # Triangulate the surface if using the mesh strategy:
  if True:#sRefined == True:
    tris = find_tris(bubblepts0)

  # Repeat the specified number of times:
  for rep in range(reps):
    axs[rep].plot_trisurf(soap.T[0], soap.T[1], soap.T[2], triangles=tris, cmap='viridis', edgecolor='none')

    # If using meshing strategy
    if sRefined == True:
      # Add the adjustments to the soap point cloud
      soap = np.add(soap, refine(soap, tris))

    # If using the boring strategy
    else:
      # If the number of points remains constant
      if constN:
        # Complete all of the steps
        for stepn in steps:
          # If the step uses "find_nearest"
          if stepn[0] == "find_nearest":
            # Add the adjustments to the soap point cloud
            soap = np.add(soap, find_nearest(soap))

          # Refine
          if stepn[0] == "refine":
            radj=refine(soap, tris, stp=stepn[1] * step)
            soap = np.add(soap, radj[0])
            print(radj[1])

          # Add points to mesh
          if stepn[0] == "add":
            # If a number of points to add is specified:
            if len(stepn) == 2:
              mesh = addpts(soap, tris, N_to_add=stepn[1])
            # Otherwise, add the default number of points
            else:
              mesh = addpts(soap, tris)
            soap = mesh[0]
            tris = mesh[1]
            # print("\n\nsoap",soap[-10:])

        #soap = np.add(soap, refine(soap, tris, stp=50*step))
      # If more points are added each time
      else:
        # Add the adjustments to the soap point cloud, and append the boundary points once more
        soap = np.add(np.concatenate([soap, boundary0], axis=0), find_nearest(soap))
    # Print number of iterations completed to keep track of progress
    print(f"Iteration: {rep + 1}")

  # Return the final result
  return [soap,tris]


vers = "3_3"
#"""

film, tris = bubble_evolve(soap0, REPS, stepsList)
for i in range(len(Z)):
  radj=refine(film, tris, stp = Z[i] * step)
  film = np.add(film, radj[0])
  print(radj[1])
# print(find_nearest(soap0))

# print("triArea", area([np.array([0,0,0]), np.array([1,0,0]), np.array([1,1,1])]), 1/sqrt(2))
#"#""
# print("film", film[0][-30:])
fig2 = plt.figure()
ax2 = plt.axes(projection='3d')
if surf:
  ax2.plot_trisurf(film.T[0], film.T[1], film.T[2], triangles=tris, cmap='viridis', edgecolor='none');
  #ax.plot_trisurf(film.T[0], film.T[1], film.T[2], triangles=find_tris(soap0), cmap='viridis', edgecolor='none');
else:
  ax2.scatter3D(film.T[0], film.T[1], film.T[2], c=[0 if i==0 else 1 for i in range(len(film))], cmap='Reds');
# Save resulting figure as an image
fig.savefig(f"model{vers}_end.jpg")

#"#""
#"#""
# print(refine([np.array([0.5, 0.5, 0.5]), np.array([1.5, 0.5, 0.5]), np.array([1.5, 1.5, 1.5])], [[0, 1, 2]]))


with open(f"out{vers}.txt", "w") as f:
  f.truncate(0)
  for pt in film:
    f.write(f"{pt[0]},{pt[1]},{pt[2]}\n")





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
# find_nearest_time()


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

# from mpl_toolkits.mplot3d import Axes3D

bubblepts=soap0

# tris = np.array([[[n,n+1,n+10], [n,n+9,n+10]] for n in range(len(bubblepts)-10) if n%10 == 9][:40]).reshape(-1,3)


# print(tris)
# ax.plot_trisurf(bubblepts[:,0], bubblepts[:,1], bubblepts[:,2], triangles = find_tris(bubblepts))
# fig.savefig("begin0.jpg")



# Commented out np.unique line
