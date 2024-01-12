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
import numpy as np
from hough import Hough
from pointcloud import PointCloud
from scipy.linalg import eigh
# orthogonal least squares fit?
# rc = largest eienvalue...
# c = np.matrix()
# c.getH()
# a = np.array([[0,1,2],[3,4,5],[6,7,8]],np.float64)
# print(a.mean(0))
# a -= np.mean(a,0)
# print(a)
def orthogonal_LSQ(pc):
  rc= 0.0
  # a = pc.meanValue() # i.e. weighted mean center point
  a = pc.meanValue()
  # a = np.random.randint(0,100,(6,3),np.float64)
  # points = np.zeros((len(pc.points),3),np.float64)
  points = np.zeros((len(pc.points),3),np.float64)
  for i,j in enumerate(pc.points):
  # for i,j in enumerate(pc.points):
    points[i] = j
  centered = np.matrix(points - a)
  # centered2 = np.matrix(points - points.mean(0))
  # print()
  # print(centered)
  # print(centered2)
  scatter = centered.getH() * centered
  # print()
  # print(scatter)
  eigenVal, eigenVec = eigh(scatter)
  # print()
  # print(eigenVal)
  # print()
  # print(eigenVec)
  b = eigenVec[:,2]
  rc = eigenVal[2]
  return rc,a,b
# pc = np.random.randint(-99,100,(6,3))#,np.float64)
# # print(pc)
# rc,a,b = orthogonal_LSQ(pc)
# print()
# print(rc)
# print(a)
# print(b)
opt_dx = 0
opt_nlines = 0
opt_minvotes = 10
opt_verbose = 0
infile_name = 'test.dat'
outfile_name = ""
granularity = 4
num_directions = [12,21,81,321,1281,5121,20481]
infile = ""
outfile = ""
# bounding box of point cloud
minP,maxP,minPshifted,maxPshifted = np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0])
d = 0 # diagonal length of point cloud
# plausibility checks
if infile_name=="":
  -print('no infile name')
if opt_dx < 0:
  -print('dx < 0')
if opt_nlines < 0:
  -print('nlines < 0')
if opt_minvotes < 0:
  -print('minvotes < 0')
if opt_minvotes < 2:
  opt_minvotes = 2
if outfile_name:
  try:
    outfile_name = open(outfile_name,'w')
  except:
    -print('cannot open outfile %s'%outfile_name)
X = PointCloud()
X.readFromFile(infile_name)
# try:X.readFromFile(infile_name)
# except:-print('cannot open infile %s'%infile_name)
if len(X.points) < 2:
  -print('point cloud has less than two points')
minP, maxP = X.getMinMax3D()
d = np.linalg.norm(maxP-minP)
if d == 0:
  -print('all points in point cloud identical')
X.shiftToOrigin()
# minPshifted, maxPshifted = X. 
minPshifted, maxPshifted = X.getMinMax3D() # this seems over complicated?
if opt_dx == 0:
  opt_dx = d / 64.0
elif opt_dx >= d:
  -print('dx too large')
num_x = math.floor(d / opt_dx + 0.5)
num_cells = num_x * num_x * num_directions[granularity]
'''
#   if (opt_verbose) {
#     printf("info: x'y' value range is %f in %.0f steps of width dx=%f\n",
#            d, num_x, opt_dx);
#     printf("info: Hough space has %.0f cells taking %.2f MB memory space\n",
#            num_cells, num_cells * sizeof(unsigned int) / 1000000.0);
#   }
# #ifdef WEBDEMO
#   if (num_cells > 1E8) {
#     fprintf(stderr, "Error: program was compiled in WEBDEMO mode, "
#             "which does not permit %.0f cells in Hough space\n", num_cells);
#     return 2;
#   }
# #endif
'''
hough = Hough(minPshifted, maxPshifted, opt_dx, granularity)
hough.add(X)

Y= PointCloud()
rc = 0
nvotes = 0
nlines = 0
while len(X.points) > 1 and (opt_nlines == 0 or opt_nlines > opt_nlines):
  a = np.zeros(3)
  b = np.zeros(3)
  hough.subtract(Y)
  nvotes,a,b = hough.getLine()
  if opt_verbose > 1:
    p = a + X.shift
    print("info: highest number of Hough votes is %i for the following line:\info a=(%f,%f,%f), b=(%f,%f,%f)"%(nvotes,p[0],p[1],p[2],b[0],b[1],b[2]))
  
  Y.points = X.pointsCloseToLine(a,b,opt_dx)

  rc,a,b = orthogonal_LSQ(Y) # optimizes a,b... somehow?
  if rc == 0:break
  
  Y.points = X.pointsCloseToLine(a,b,opt_dx)
  nvotes = len(Y.points)
  if nvotes < opt_minvotes:break

  rc,a,b = orthogonal_LSQ(Y)
  if rc == 0:break

  a += X.shift

  nlines += 1
  print("npoints=%i, a=(%f,%f,%f), b=(%f,%f,%f)"%(len(Y.points),a[0],a[1],a[2],b[0],b[1],b[2]))#,file=outfile)

  X.removePoints(Y)
  # if (opt_outformat == format_normal) {
  #   fprintf(outfile, "npoints=%lu, a=(%f,%f,%f), b=(%f,%f,%f)\n",
  #           Y.points.size(), a.x, a.y, a.z, b.x, b.y, b.z);
  # }
  # else if (opt_outformat == format_gnuplot) {
  #   fputs(", \\\n    ", outfile);
  #   fprintf(outfile, "%f + u * %f, %f + u * %f, %f + u * %f "
  #           "with lines notitle lc rgb 'black'",
  #           a.x, b.x, a.y, b.y, a.z, b.z);
  # }
  # else {
  #   fprintf(outfile, "%f %f %f %f %f %f %lu\n",
  #           a.x, a.y, a.z, b.x, b.y, b.z, Y.points.size());
  # }