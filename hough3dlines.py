'''
MODIFIED FROM:
 hough3dlines.cpp
     Main program implementing the iterative Hough transform

 Author:  Tilman Schramke, Christoph Dalitz
 Date:    2018-03-24
 License: see License-BSD2
'''

import math

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eigh

from sphere import sphere_vertices

# NOTE: It may be better to fix opt_dx based on a singular application of the TPC - right now it is based on the extremes of the pointCloud but it may be better served to 
# NOTE: It is more robust to allow the caller of hough3d to plot the points / grouped points / lines manually based on the return rather than having a robust set of plotting parameters

def hough3d(point_cloud=[],opt_dx=0,sphere_granularity=4,opt_nlines=0,opt_minvotes=10,opt_verbose=0,infile_name='',return_num_points=True,return_line_params=False,return_indices=True,return_points=False,save_test_vis_fig=False,test_vis_fig_file_name = 'testVis.png'):
  """Main Run Call"""
  # Plausibility checks
  point_cloud = np.array(point_cloud,np.float64)
  if len(point_cloud) and infile_name:print('Manual points and a file name cannot both be passed to hough3d!');return{}
  if opt_dx < 0:print('opt_dx < 0');return{}
  if opt_nlines < 0:print('opt_nlines < 0');return{}
  if opt_minvotes < 0:print('opt_minvotes < 0');return{}
  if 0 < opt_minvotes < 2:print('opt_minvotes < 2 (must be >= 2)');return{}
  if infile_name:
    try:
      point_cloud = np.loadtxt(infile_name,delimiter=',')
      if opt_verbose >= 1:print(f'Loaded {len(point_cloud)} points from file')
    except:print(f'Cannot open infile {infile_name}');return{}
    if len(point_cloud) < 2:print('point cloud has less than two points');return{}
  
  X = point_cloud # Initiatize action PointCloud X
  # point_cloud = list(point_cloud) # Will later allow for pointCloud.index() to be used easily

  B = np.array(sphere_vertices(sphere_granularity)) # Load direction vectors using sphere.py
  num_b = len(B)
  if opt_verbose >= 2:
    print(f'Successfully loaded {num_b} vertices')

  num_points = []
  line_params = []
  indices = []
  points =[]

  def orthogonal_lsq(pc):
    if len(pc): # do not take mean of empty array
      a = np.mean(pc,0)
    else:a = np.zeros(3)
    centered = np.matrix(pc - a)
    scatter = centered.getH() * centered
    eigen_val, eigen_vec = eigh(scatter)
    b = eigen_vec[:,2]
    rc = eigen_val[2]
    return rc,a,b

  max_bound_point = np.max(X,0) # max corner # maxPshifted
  min_bound_point = np.min(X,0) # min corner # minPshifted

  shift = (max_bound_point + min_bound_point)/2

  X -= shift

  max_bound_pointshifted = np.max(X,0) # max corner # maxPshifted
  min_bound_pointshifted = np.min(X,0) # min corner # minPshifted

  extent_bound_box = np.linalg.norm(max_bound_pointshifted-min_bound_pointshifted) # d

  if opt_dx == 0:
    dx = extent_bound_box / 64
  elif opt_dx > extent_bound_box:
    -print('dx too large')
  else:
    dx = opt_dx

  index_shift_n = math.ceil((extent_bound_box/dx-1)/2) # N = 1 + 2n, n >= (d/dx-1)/2 [from (2n+1)dx >= d] # Side length of lattice in number of cells # Half of extent[rounded up to odd multiple of dx] over dx
  lattice_side_length_anchor_points = 1+2*index_shift_n # Number
  anchor_point_shift = index_shift_n * dx # End to origin (Anchor Points)
  lattice_bounding_box_side_length = (lattice_side_length_anchor_points+1) * dx # End to end length of x' space (= previous + dx)

  def x_p_y_p_from_index(x_p_ind,y_p_ind):
    return dx*(x_p_ind - index_shift_n), dx*(y_p_ind - index_shift_n)
  def points_close_to_line(X,a,b,dx):
    return X[np.where(np.linalg.norm(X - (a + np.array([i*b for i in np.dot(X-a,b)])),axis=1) <= dx)]
  def point_vote(X,add=True):
    c = [-1,1][add]
    for x in X:
      for b_ind,b in enumerate(B):
        beta = 1/(1+b[2])
        x_p = (1-b[0]**2*beta)*x[0] - b[0]*b[1]*beta*x[1] - b[0]*x[2]
        y_p = -b[0]*b[1]*x[0] + (1-b[1]**2*beta)*x[1] -b[1]*x[2]
        x_p_ind, y_p_ind = round(x_p/dx)+index_shift_n, round(y_p/dx)+index_shift_n
        A[b_ind,x_p_ind,y_p_ind] += c

  A = np.zeros((len(B),lattice_side_length_anchor_points,lattice_side_length_anchor_points)) # Accumulator Array A, |B| X |X'| X |Y'| 

  point_vote(X,True)

  nlines = 0
  Y = []

  if save_test_vis_fig:
    fig = plt.figure(figsize=(15,15))
    ax=plt.subplot(projection='3d')
    # p11 = ax.scatter(X[:,0]+shift[0],X[:,1]+shift[1],X[:,2]+shift[2],c=X[:,2]+shift[2],cmap='gnuplot')
    # ax.plot([],[],[])

  while len(X) > 1 and (opt_nlines == 0 or opt_nlines > nlines):
    
    point_vote(Y,False) # Subtract Y votes from A

    nvotes = np.max(A) # Max vote in A
    b_ind, x_p_ind, y_p_ind = np.unravel_index(A.argmax(),A.shape) # Index array of max vote in A
    x_p, y_p = x_p_y_p_from_index(x_p_ind,y_p_ind) # Gets discretized (x_j',y_k') from indices
    b = B[b_ind] # Direction Vector
    a = x_p*np.array([1-b[0]**2/(1+b[2]),-b[0]*b[1]/(1+b[2]),-b[0]]) + y_p*np.array([-b[0]*b[1]/(1+b[2]),1-b[1]**2/(1+b[2]),-b[1]]) # Anchor point
    if opt_verbose >= 2: # print voting information
      p = a + shift
      print(f"info: highest number of Hough votes is {nvotes} for the following line:\info a=({p[0]},{p[1]},{p[2]}), b=({b[0]},{b[1]},{b[2]})")

    Y = points_close_to_line(X,a,b,dx) # Select points close to Heuristic line 1
    rc, a, b = orthogonal_lsq(Y) # Perform OLS fit on selected points
    if rc == 0: # ?
      if opt_verbose >= 1:
        print('rc == 0')
      break

    Y = points_close_to_line(X,a,b,dx) # Select points close to Heuristic line 2
    nvotes = len(Y)
    if nvotes < opt_minvotes:print('nvotes < opt_minvotes',nvotes,opt_minvotes);break # End loop if nvotes is less than opt_minvotes

    rc,a,b = orthogonal_lsq(Y) # Perform OLS fit on selected points to best describe the line

    for i in Y:X=X[np.where(1-np.all(X==i,1))] # Remove points in Y from X
    nlines += 1

    a += shift # shift anchor point before storing line information
    if return_num_points:
      num_points += [len(Y)]
    if return_line_params:
      line_params += [[a,b]]
    if return_points:
      points += [Y+shift]
    if return_indices:
      indices += [[np.where(np.all(point_cloud==i,1))[0][0]for i in Y]]
    if opt_verbose >= 1:
      print(f"npoints={len(Y)}, a=({a[0]},{a[1]},{a[2]}), b=({b[0]},{b[1]},{b[2]})")
    if save_test_vis_fig:
      ax.autoscale(True)
      T0 = -extent_bound_box*5
      T1 = extent_bound_box*5
      if b[0]:
        t0,t1 = sorted([(min_bound_point[0]-a[0])/b[0],(max_bound_point[0]-a[0])/b[0]])
        T0,T1 = max([T0,t0]),min([T1,t1])
      if b[1]:
        t0,t1 = sorted([(min_bound_point[1]-a[1])/b[1],(max_bound_point[1]-a[1])/b[1]])
        T0,T1 = max([T0,t0]),min([T1,t1])
      if b[2]:
        t0,t1 = sorted([(min_bound_point[2]-a[2])/b[2],(max_bound_point[2]-a[2])/b[2]])
        T0,T1 = max([T0,t0]),min([T1,t1])
      p1 = np.array([a+b*T0,a+b*T1])
      ax.plot(p1[:,0],p1[:,1],p1[:,2])
      p1 = Y + shift
      ax.scatter(p1[:,0],p1[:,1],p1[:,2])
  X += shift
  output = {}
  if return_num_points:
    output['num_points'] = num_points
  if return_line_params:
    output['line_params'] = line_params
  if return_points:
    output['points'] = points
  if return_indices:
    output['indices'] = indices
  if save_test_vis_fig:
    ax.scatter(X[:,0],X[:,1],X[:,2],c='black')
    # cax = fig.add_axes([0.9, 0.5, 0.05, 0.2])
    # p = fig.colorbar(p11,cax=cax)
    # p.update_ticks()
    ax.set_xlim(min_bound_point[0],max_bound_point[0])
    ax.set_ylim(min_bound_point[1],max_bound_point[1])
    ax.set_zlim(min_bound_point[2],max_bound_point[2])
    ax.set_proj_type('ortho')
    plt.savefig(test_vis_fig_file_name)
  return output

if __name__ == '__main__':
  output = hough3d(infile_name='test.dat',return_num_points=True,return_line_params=True,return_points=True,return_indices=True,save_test_vis_fig=True,test_vis_fig_file_name='testVis5_0.png',sphere_granularity=4,opt_verbose=0)
  if len(output.keys())==0:print('output is empty');exit()
  for i in output.keys():
    print(i)
    [print(j)for j in output[i]]
    print()