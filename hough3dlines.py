'''
MODIFIED FROM:
 hough3dlines.cpp
     Main program implementing the iterative Hough transform

 Author:  Tilman Schramke, Christoph Dalitz
 Date:    2018-03-24
 License: see License-BSD2
'''

# import pointcloud
# import hough
import math
import matplotlib.pyplot as plt
import numpy as np
# from hough import Hough
# from pointcloud import PointCloud
from scipy.linalg import eigh
from sphere import sphereVertices

# NOTE: It may be better to fix opt_dx based on a singular application of the TPC - right now it is based on the extremes of the pointCloud but it may be better served to 
# NOTE: It is more robust to allow the caller of hough3d to plot the points / grouped points / lines manually based on the return rather than having a robust set of plotting parameters


def hough3d(pointCloud=[],opt_dx=0,sphereGranularity=4,opt_nlines=0,opt_minvotes=10,opt_verbose=0,infile_name='',retLineParams=False,retNumPoints=True,retIndices=True,retPoints=False,saveTestVisFig=False,testVisFigFile_name = 'testVis.png'):
  
  # Plausibility checks
  if pointCloud and infile_name:print('Manually points and a file name cannot both be passed to hough3d!');return
  if opt_dx < 0:print('opt_dx < 0');return
  if opt_nlines < 0:print('opt_nlines < 0');return
  if opt_minvotes < 0:print('opt_minvotes < 0');return
  if 0 < opt_minvotes < 2:print('opt_minvotes < 2 (must be >= 2)');return
  if infile_name:
    try:
      pointCloud = np.loadtxt(infile_name,delimiter=',')
      if opt_verbose >= 1:
        print('Loaded %i points from file'%len(pointCloud))
    except:
      print('Cannot open infile %s'%infile_name)
      return
  if len(pointCloud) < 2:print('point cloud has less than two points');return
  
  X = pointCloud
  pointCloud = list(pointCloud)

  B = sphereVertices(sphereGranularity)
  num_b = len(B)
  if opt_verbose >= 2:print('Successfully loaded %i vertices'%num_b)

  line_params = []
  indices = []
  points =[]

  def orthogonal_LSQ(pc):
    if len(pc):a = np.mean(pc,0)
    else:a = np.zeros(3)
    centered = np.matrix(pc - a)
    scatter = centered.getH() * centered
    eigenVal, eigenVec = eigh(scatter)
    b = eigenVec[:,2]
    rc = eigenVal[2]
    return rc,a,b

  maxBoundPoint = np.max(X,0) # max corner # maxPshifted
  minBoundPoint = np.min(X,0) # min corner # minPshifted

  shift = (maxBoundPoint + minBoundPoint)/2

  X -= shift

  maxBoundPointshifted = np.max(X,0) # max corner # maxPshifted
  minBoundPointshifted = np.min(X,0) # min corner # minPshifted

  extentBoundBox = np.linalg.norm(maxBoundPointshifted-minBoundPointshifted) # d

  if opt_dx == 0:
    dx = extentBoundBox / 64
  elif opt_dx > extentBoundBox:
    -print('dx too large')
  else:
    dx = opt_dx

  indexShift_n = math.ceil((extentBoundBox/dx-1)/2) # N = 1 + 2n, n >= (d/dx-1)/2 [from (2n+1)dx >= d] # Side length of lattice in number of cells # Half of extent[rounded up to odd multiple of dx] over dx
  latticeSideLengthAnchorPoints = 1+2*indexShift_n # Number
  anchorPointShift = indexShift_n * dx # End to origin (Anchor Points)
  latticeBoundingBoxSideLength = (latticeSideLengthAnchorPoints+1) * dx # End to end length of x' space (= previous + dx)

  def roundToNearest(num):
    return math.floor(num+0.5) if num>0 else math.ceil(num-0.5)
  def xPrime_yPrime_from_index(xPrime_ind,yPrime_ind):
    return dx*(xPrime_ind - indexShift_n), dx*(yPrime_ind - indexShift_n)
  def index_from_xPrime_yPrime_and_round(xPrime,yPrime):
    return roundToNearest(xPrime/dx) + indexShift_n,roundToNearest(yPrime/dx) + indexShift_n
  def pointsCloseToLine(X,a,b,dx):
    b = np.array(b)
    return X[np.where(np.linalg.norm(X - (a + np.array([i*b for i in np.dot(X-a,b)])),axis=1) <= dx)]

  A = np.zeros((len(B),latticeSideLengthAnchorPoints,latticeSideLengthAnchorPoints)) # Accumulator Array A, |B| X |X'| X |Y'| 

  for x in X:
    for b_ind,b in enumerate(B):
      beta = 1/(1+b[2])
      xP = (1-b[0]**2*beta)*x[0] - b[0]*b[1]*beta*x[1] - b[0]*x[2]
      yP = -b[0]*b[1]*x[0] + (1-b[1]**2*beta)*x[1] -b[1]*x[2]
      xPrime_ind, yPrime_ind = index_from_xPrime_yPrime_and_round(xP,yP)
      A[b_ind,xPrime_ind,yPrime_ind] += 1

  nlines = 0
  Y = []

  if saveTestVisFig:
    fig = plt.figure(figsize=(15,15))
    ax=plt.subplot(projection='3d')
    # p11 = ax.scatter(X[:,0]+shift[0],X[:,1]+shift[1],X[:,2]+shift[2],c=X[:,2]+shift[2],cmap='gnuplot')
    # ax.plot([],[],[])

  while len(X) > 1 and (opt_nlines == 0 or opt_nlines > nlines):
    
    for x in Y:
      for b_ind,b in enumerate(B):
        beta = 1/(1+b[2])
        xP = (1-b[0]**2*beta)*x[0] - b[0]*b[1]*beta*x[1] - b[0]*x[2]
        yP = -b[0]*b[1]*x[0] + (1-b[1]**2*beta)*x[1] -b[1]*x[2]
        xPrime_ind, yPrime_ind = index_from_xPrime_yPrime_and_round(xP,yP)
        A[b_ind,xPrime_ind,yPrime_ind] -= 1

    nvotes = np.max(A)
    b_ind, xPrime_ind, yPrime_ind = np.unravel_index(A.argmax(),A.shape)
    xP, yP = xPrime_yPrime_from_index(xPrime_ind,yPrime_ind)
    b = B[b_ind]
    a = xP*np.array([1-b[0]**2/(1+b[2]),-b[0]*b[1]/(1+b[2]),-b[0]]) + yP*np.array([-b[0]*b[1]/(1+b[2]),1-b[1]**2/(1+b[2]),-b[1]])
    if opt_verbose >= 2: # print voting information
      p = a + shift
      print("info: highest number of Hough votes is %i for the following line:\info a=(%f,%f,%f), b=(%f,%f,%f)"%(nvotes,p[0],p[1],p[2],b[0],b[1],b[2]))

    Y = pointsCloseToLine(X,a,b,dx)
    rc, a, b = orthogonal_LSQ(Y)
    if rc == 0:print('rc == 0');break

    Y = pointsCloseToLine(X,a,b,dx)
    nvotes = len(Y)
    if nvotes < opt_minvotes:print('nvotes < opt_minvotes',nvotes,opt_minvotes);break

    rc,a,b = orthogonal_LSQ(Y)

    for i in Y:X=X[np.where(1-np.all(X==i,1))]
    nlines += 1

    a += shift
    if retLineParams:
      if retNumPoints:line_params += [[len(Y),a,b]]
      else:line_params += [[a,b]]
    if retPoints:
      points += [Y]
    if retIndices:
      indices += [[pointCloud.index(i)for i in Y]]
    if opt_verbose >= 1:
      print("npoints=%i, a=(%f,%f,%f), b=(%f,%f,%f)"%(len(Y),a[0],a[1],a[2],b[0],b[1],b[2]))#,file=outfile)
    if saveTestVisFig:
      ax.autoscale(True)
      T0 = -extentBoundBox*5
      T1 = extentBoundBox*5
      if b[0]:
        t0,t1 = sorted([(minBoundPoint[0]-a[0])/b[0],(maxBoundPoint[0]-a[0])/b[0]])
        T0,T1 = max([T0,t0]),min([T1,t1])
      if b[1]:
        t0,t1 = sorted([(minBoundPoint[1]-a[1])/b[1],(maxBoundPoint[1]-a[1])/b[1]])
        T0,T1 = max([T0,t0]),min([T1,t1])
      if b[2]:
        t0,t1 = sorted([(minBoundPoint[2]-a[2])/b[2],(maxBoundPoint[2]-a[2])/b[2]])
        T0,T1 = max([T0,t0]),min([T1,t1])
      p1 = np.array([a+b*T0,a+b*T1])
      ax.plot(p1[:,0],p1[:,1],p1[:,2])
      p1 = Y + shift
      ax.scatter(p1[:,0],p1[:,1],p1[:,2])
  X += shift
  if saveTestVisFig:
    # print(len(X))
    # print(X[:2])
    # ax.scatter(X[:,0],X[:,1],X[:,2],c='black')
    # cax = fig.add_axes([0.9, 0.5, 0.05, 0.2])
    # p = fig.colorbar(p11,cax=cax)
    # p.update_ticks()
    ax.set_xlim(minBoundPoint[0],maxBoundPoint[0])
    ax.set_ylim(minBoundPoint[1],maxBoundPoint[1])
    ax.set_zlim(minBoundPoint[2],maxBoundPoint[2])
    ax.set_proj_type('ortho')
    plt.savefig(testVisFigFile_name)

if __name__ == '__main__':
  hough3d(infile_name='testdata.dat',retLineParams=True,retNumPoints=True,retIndices=False,saveTestVisFig=True,testVisFigFile_name='testVis6_0.png',sphereGranularity=6,opt_verbose=0)