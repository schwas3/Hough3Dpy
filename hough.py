'''
MODIFIED FROM hough.cpp / hough.h
    Author:  Tilman Schramke, Manuel Jeltsch, Christoph Dalitz
    Date:    2017-03-16
    License: see LICENSE-BSD2
'''

# def test(a):
#     a = 5
# a = 1
# print(a)
# test(a)
# print(a)

import math
import numpy as np
# import sphereVertices from sphere
from sphere import sphereVertices

def roundToNearest(num):
    return math.floor(num+.5) if num>0 else math.ceil(num-.5)

class Hough:
    def __init__(self, minP=np.array([0,0,0]), maxP=np.array([0,0,0]), var_dx=0, sphereGranularity=0):
        self.minP = minP
        self.maxP = maxP
        self.var_dx = var_dx
        self.vertices = sphereVertices(sphereGranularity)
        # self.sphere = sphere.Sphere(sphereGranularity)
        # self.sphere.fromIcosahedron(sphereGranularity)
        
        self.num_b = len(self.vertices)

        self.max_x = max(np.linalg.norm(maxP), np.linalg.norm(minP))
        range_x = 2 * self.max_x
        dx = var_dx
        if dx==0:
            dx = range_x / 64.0
        self.dx = dx
        self.num_x = roundToNearest(range_x / dx)

        self.VotingSpace = np.zeros(self.num_x*self.num_x*self.num_b)

        self.a=np.array([0,0,0])
        self.b=np.array([0,0,0])

        # ~Hough: delete self.sphere ## Necessary? I don't get why this is done...

    def add(self,pc=np.zeros((3,3))):
        for i in pc.points:
            self.pointVote(i,True)
    def subtract(self,pc=np.zeros((3,3))):
        for i in pc.points:
            self.pointVote(i,False)
    def pointVote(self, point = np.array([0,0,0]), add = True):
        for j,b in enumerate(self.vertices):
            # b = self.vertices[j]
            beta = 1 / (1 + b[2])
            x_new = (((1 - (beta * b[0]**2))) * point[0]) - ((beta * (b[0] * b[1])) * point[1]) - (b[0] * point[2])
            y_new = ((-beta * (b[0] * b[1])) * point[0]) + ((1 - (beta * (b[1] * b[1]))) * point[1]) - (b[1] * point[2])
            x_i = roundToNearest((x_new + self.max_x) / self.dx)
            y_i = roundToNearest((y_new + self.max_x) / self.dx)
            index = (x_i * self.num_x * self.num_b) + (y_i * self.num_b) + j
            if index < len(self.VotingSpace):
                if add:
                    self.VotingSpace[index] += 1
                else:
                    self.VotingSpace[index] -= 1
    def getLine(self):
        votes = 0
        index = 0
        for i,j in enumerate(self.VotingSpace):
            if j > votes:
                votes = j
                index = i
        x = int(index / (self.num_x * self.num_b))
        index -= (int) (x * self.num_x * self.num_b)
        x = x * self.dx - self.max_x
        y = int(index / self.num_b)
        index -= (int)(y * self.num_b)
        y = y * self.dx - self.max_x
        b = self.vertices[index]
        a = np.array([
            x * (1 - ((b[0] * b[0]) / (1 + b[2]))) - y * ((b[0] * b[1]) / (1 + b[2])),
            x * (-((b[0] * b[1]) / (1 + b[2]))) + y * (1 - ((b[1] * b[1]) / (1 + b[2]))),
            -x * b[0] - y * b[1]
        ]
        )
        return votes,a,b

{# static double roundToNearest(double num) {
#   return (num > 0.0) ? floor(num + 0.5) : ceil(num - 0.5);
# }

# Hough::Hough(const Vector3d& minP, const Vector3d& maxP, double var_dx,
#              unsigned int sphereGranularity) {

#   // compute directional vectors
#   sphere = new Sphere();
#   sphere->fromIcosahedron(sphereGranularity);
#   num_b = sphere->vertices.size();

#   // compute x'y' discretization
#   max_x = std::max(maxP.norm(), minP.norm());
#   double range_x = 2 * max_x;
#   dx = var_dx;
#   if (dx == 0.0) {
#     dx = range_x / 64.0;
#   }
#   num_x = roundToNearest(range_x / dx);

#   // allocate voting space
#   VotingSpace.resize(num_x * num_x * num_b);
# }

# Hough::~Hough() {
#   delete sphere;
# }
# // add all points from point cloud to voting space
# void Hough::add(const PointCloud &pc) {
#   for (std::vector<Vector3d>::const_iterator it = pc.points.begin();
#        it != pc.points.end(); it++) {
#     pointVote((*it), true);
#   }
# }

# // subtract all points from point cloud to voting space
# void Hough::subtract(const PointCloud &pc) {
#   for (std::vector<Vector3d>::const_iterator it = pc.points.begin();
#        it != pc.points.end(); it++) {
#     pointVote((*it), false);
#   }
# }

# // add or subtract (add==false) one point from voting space
# // (corresponds to inner loop of Algorithm 2 in IPOL paper)
# void Hough::pointVote(const Vector3d& point, bool add){

#   // loop over directions B
#   for(size_t j = 0; j < sphere->vertices.size(); j++) {

#     Vector3d b = sphere->vertices[j];
#     double beta = 1 / (1 + b.z);	// denominator in Eq. (2)

#     // compute x' according to left hand side of Eq. (2)
#     double x_new = ((1 - (beta * (b.x * b.x))) * point.x)
#       - ((beta * (b.x * b.y)) * point.y)
#       - (b.x * point.z);

#     // compute y' according to right hand side Eq. (2)
#     double y_new = ((-beta * (b.x * b.y)) * point.x)
#       + ((1 - (beta * (b.y * b.y))) * point.y)
#       - (b.y * point.z);

#     size_t x_i = roundToNearest((x_new + max_x) / dx);
#     size_t y_i = roundToNearest((y_new + max_x) / dx);

#     // compute one-dimensional index from three indices
# 	// x_i * #planes * #direction_Vec + y_i * #direction_Vec + #loop
#     size_t index = (x_i * num_x * num_b) + (y_i * num_b) + j;

#     if (index < VotingSpace.size()) {
#       if(add){
#         VotingSpace[index]++;
#       } else {
#         VotingSpace[index]--;
#       }
#     }
#   }
# }
}

# // returns the line with most votes (rc = number of votes)
# unsigned int Hough::getLine(Vector3d* a, Vector3d* b){
#   unsigned int votes = 0;
#   unsigned int index = 0;

#   for(unsigned int i = 0; i < VotingSpace.size(); i++){
#     if (VotingSpace[i] > votes) {
#       votes = VotingSpace[i];
#       index = i;
#     }
#   }

#   // retrieve x' coordinate from VotingSpace[num_x * num_x * num_b]
#   double x = (int) (index / (num_x * num_b));
#   index -= (int) x * num_x * num_b;
#   x = x * dx - max_x;

#   // retrieve y' coordinate from VotingSpace[num_x * num_x * num_b]
#   double y = (int) index / num_b;
#   index -= (int) y * num_b;
#   y = y * dx - max_x;

#   // retrieve directional vector
#   *b = sphere->vertices[index];

#   // compute anchor point according to Eq. (3)
#   a->x = x * (1 - ((b->x * b->x) / (1 + b->z)))
#     - y * ((b->x * b->y) / (1 + b->z));
#   a->y = x * (-((b->x * b->y) / (1 + b->z)))
#      + y * (1 - ((b->y * b->y) / (1 + b->z)));
#   a->z = - x * b->x - y * b->y;

#   return votes;
# }
