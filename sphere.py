'''
Modified from:
 sphere.cpp
     Class that implements direction quantization as described in
     Jeltsch, Dalitz, Pohle-Froehlich: "Hough Parameter Space
     Regularisation for Line Detection in 3D." VISAPP, pp. 345-352, 2016 

 Author:  Manuel Jeltsch, Tilman Schramke, Christoph Dalitz
 Date:    2017-03-16
 License: see LICENSE-BSD2
'''

### Optimization notes
# Make saved files for tesselation / subdividing (i.e. save numerical values rather than regenerating them) - will test runtime
# Remove unecessary triangles early on... - i.e. if all three points will be cut eventually be makeUnique(), they might as well be cut early
# Removing triangles in makeUnique() seems COMPLETELY unecessary - this would only be needed if you wanted to implement something like the line directly above, but even still this is not quite what you'd want
# Cont. - this is not quite right because you want to remove triangles with all vertices to be cut, but this will still be only something like a factor of 2x faster

import numpy as np
import pickle as pkl
import os

# also optimizable

def sphereVertices(subDivisions=4,returnVertsAndTriangles=False,saveToFile=True,saveToFile_name='sphere.pkl'):
  if saveToFile and os.path.exists(saveToFile_name):
    with open(saveToFile_name,'rb') as f:
      sphere_dict = pkl.load(f)
    vertices,triangles = sphere_dict['vertices'][-1], sphere_dict['triangles'][-1]
    granularity = len(sphere_dict['vertices']) - 1
    if granularity >= subDivisions:
      if returnVertsAndTriangles:
        return sphere_dict['vertices'][subDivisions],sphere_dict['triangles'][subDivisions]
      else:
        return sphere_dict['uniqueVertices'][subDivisions]
  else:
    tau = 1.61803399; # golden_ratio
    norm = (1 + tau * tau)**.5
    v = 1 / norm
    tau = tau / norm
    vertices = [[-v,tau,0],[v,tau,0],[0,v,-tau],[0,v,tau],[-tau,0,-v],[tau,0,-v],[-tau,0,v],[tau,0,v],[0,-v,-tau],[0,-v,tau],[-v,-tau,0],[v,-tau,0]]
    triangles = [[0,1,2],[0,1,3],[0,2,4],[0,4,6],[0,3,6],[1,2,5],[1,3,7],[1,5,7],[2,4,8],[2,5,8],[3,6,9],[3,7,9],[4,8,10],[8,10,11],[5,8,11],[5,7,11],[7,9,11],[9,10,11],[6,9,10],[4,6,10]]
    if saveToFile:
      uniqueVertices = [vert for vert in vertices if not (vert[2] < 0 or (vert[2] == 0 and vert[0] < 0) or (vert[2] == 0 == vert[0] and vert[1] == -1))]
      sphere_dict = {'uniqueVertices':[uniqueVertices],'vertices':[vertices],'triangles':[triangles]}
    granularity = 0
  while granularity < subDivisions:
    granularity += 1
    triangles2 = []
    for ai,bi,ci in triangles:
      # this section can likely be optimized
      a=vertices[ai]
      b=vertices[bi]
      c=vertices[ci]
      d = [a[j]+b[j]for j in range(3)]
      e = [b[j]+c[j]for j in range(3)]
      f = [c[j]+a[j]for j in range(3)]
      d = [d[j]/np.linalg.norm(d)for j in range(3)]
      e = [e[j]/np.linalg.norm(e)for j in range(3)]
      f = [f[j]/np.linalg.norm(f)for j in range(3)]
      if d in vertices:
        di = vertices.index(d)
      else:
        di = len(vertices)
        vertices += [d]
      if e in vertices:
        ei = vertices.index(e)
      else:
        ei = len(vertices)
        vertices += [e]
      if f in vertices:
        fi = vertices.index(f)
      else:
        fi = len(vertices)
        vertices += [f]
      triangles2 += [[ai,di,fi],[di,bi,ei],[fi,ei,ci],[fi,di,ei]]
    triangles = triangles2[:]
    if saveToFile:
      sphere_dict['vertices'] += [vertices]
      sphere_dict['triangles'] += [triangles]
      sphere_dict['uniqueVertices'] += [[vert for vert in vertices if not (vert[2] < 0 or (vert[2] == 0 and vert[0] < 0) or (vert[2] == 0 == vert[0] and vert[1] == -1))]]
  if saveToFile:
    with open(saveToFile_name,'wb') as f:
      pkl.dump(sphere_dict,f)
    if returnVertsAndTriangles:
      return sphere_dict['vertices'][subDivisions],sphere_dict['triangles'][subDivisions]
    else:
      return sphere_dict['uniqueVertices'][subDivisions]
  else:
    if returnVertsAndTriangles:
      return vertices,triangles
    else:
      return [vert for vert in vertices if not (vert[2] < 0 or (vert[2] == 0 and vert[0] < 0) or (vert[2] == 0 == vert[0] and vert[1] == -1))]

if __name__ == "__main__":
  # demonstrates sphere.py when ran
  # to be added - runLine arguments to allow for customization of granularity
  # to be added - runLine arguments to allow for customization of out file name

  import matplotlib.pyplot as plt
  fig = plt.figure()
  ax = fig.add_subplot(projection='3d')
  vertices = sphereVertices(6)
  print(len(vertices))
  # vertices, triangles = sphereVertices(7,True)
  # vertices = np.array(vertices)
  # # print(len(triangles))
  # print(len(vertices))
  # # ax.
  # # for i in range(3): 
  # #   ax=fig.add_subplot(1,3,i+1,projection='3d')
  # # ax.scatter(vertices[:,0])
  # # for a,b,c in vertices:
  # #   ax.scatter(a,b,c,alpha=0.1)
  # for triangle in triangles:
  #   verts = np.array(vertices)[triangle]
  #   # print(verts)
  #   ax.scatter(verts[:,0],verts[:,1],verts[:,2],c='r',s=1)
  #   ax.plot(verts[:,0],verts[:,1],verts[:,2],c='r')
  # # ax.scatter(vertices[:,0],vertices[:,1],vertices[:,2],alpha=1)
  # plt.savefig('test.png')
  # # for i in [0,2,4]: