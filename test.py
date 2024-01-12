import numpy as np
a = np.loadtxt('testdata.dat',delimiter=',')
# a = np.fromfile('./testdata.dat')
# a = np.transpose([a[::3],a[1::3],a[2::3]])
print(a)
for i in a[-5:]:a = a[np.where(1-np.all(a==i,1))]
print(a)
# print(a[np.where(np.all(a==a[-5:]))])
print(a.shape)

# a = np.array([[0,1,2],[3,4,5],[6,7,8],[5,5,5]])
# b = [[0,1,2]]
# b=b[0]
# # print(np.where(a==np.array([1,2,1])))
# print(a.min(0))