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

def sphereVertices(subDivisions=4,returnTriangles=False):
  tau = 1.61803399; # golden_ratio
  norm = (1 + tau * tau)**.5
  v = 1 / norm
  tau = tau / norm
  vertices = [[-v,tau,0],[v,tau,0],[0,v,-tau],[0,v,tau],[-tau,0,-v],[tau,0,-v],[-tau,0,v],[tau,0,v],[0,-v,-tau],[0,-v,tau],[-v,-tau,0],[v,-tau,0]]
  triangles = triangles2 = [[0,1,2],[0,1,3],[0,2,4],[0,4,6],[0,3,6],[1,2,5],[1,3,7],[1,5,7],[2,4,8],[2,5,8],[3,6,9],[3,7,9],[4,8,10],[8,10,11],[5,8,11],[5,7,11],[7,9,11],[9,10,11],[6,9,10],[4,6,10]]
  for i in range(subDivisions):
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
  if returnTriangles:
    return vertices,triangles2
  for vert in vertices[:]:
    if vert[2] < 0:
      vertices.remove(vert)
    elif vert[2] == 0:
      if vert[0] < 0:
        vertices.remove(vert)
      elif vert[0] == 0 and vert[1] == -1:
        vertices.remove(vert)
  return vertices

if __name__ == "__main__":

  # demonstrates sphere.py when ran
  # to be added - runLine arguments to allow for customization of granularity
  # to be added - runLine arguments to allow for customization of out file name

  import matplotlib.pyplot as plt
  fig = plt.figure()
  ax = fig.add_subplot(projection='3d')
  vertices = sphereVertices(4,False)
  -print(len(vertices))
  vertices, triangles = sphereVertices(subDivisions=2,returnTriangles=True)
  vertices = np.array(vertices)
  # print(len(triangles))
  # print(vertices)
  # ax.
  # for i in range(3): 
  #   ax=fig.add_subplot(1,3,i+1,projection='3d')
  # ax.scatter(vertices[:,0])
  # for a,b,c in vertices:
  #   ax.scatter(a,b,c,alpha=0.1)
  for triangle in triangles:
    verts = np.array(vertices)[triangle]
    # print(verts)
    ax.plot3D(verts[:,0],verts[:,1],verts[:,2],c='r')
  # ax.scatter(vertices[:,0],vertices[:,1],vertices[:,2],alpha=1)
  plt.savefig('test.png')
  # for i in [0,2,4]:

