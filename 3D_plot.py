from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pandas

points = pandas.read_csv('3dplot.csv')
print(points['x'].values)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

y = points['z'].values
x = points['y'].values
z = points['x'].values

#  x   y   z  
# 1.1,1.2,1.3
# 2.1,2.2,2.3
# 3.1,3.2,3.3
# 4.1,4.2,4.3
ax.set_xlabel('dp/dp_b ')
ax.set_ylabel('Q/Q_b')
ax.set_zlabel('Reward')
ax.scatter(x, y, z, c='b', marker='o')
plt.show()
plt.savefig('foo.png')
