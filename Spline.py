import numpy as np
from scipy import interpolate
from matplotlib import pyplot as plt

x = np.array([1, 2, 3, 4, 5])
y = np.array([1, 2, 3, 2, 1.5])

# getting starting points from x , y 
x = np.r_[x, x[0]]
y = np.r_[y, y[0]]

# fit splines to x=f(u) and y=g(u) 
# s=0 for fit to point
tck, u = interpolate.splprep([x, y], s=0, per=True)

# How many points path? 100
xi, yi = interpolate.splev(np.linspace(0, 1, 100), tck)

# plot the result
fig, ax = plt.subplots(1, 1)
ax.plot(x, y, 'or')
ax.plot(xi, yi, '-b')
plt.show()
plt.savefig('goo.png')
