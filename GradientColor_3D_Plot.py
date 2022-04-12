import matplotlib
import pandas

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D

points = pandas.read_csv('3dplot.csv')

print(points['x'].values)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

y = points['z'].values
x = points['y'].values
z = points['x'].values

cs = range(1,499)


def scatter3d(x,y,z, cs, colorsMap='jet'):
    cm = plt.get_cmap(colorsMap)
    cNorm = matplotlib.colors.Normalize(vmin=min(cs), vmax=max(cs))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(x, y, z, c=scalarMap.to_rgba(cs))
    scalarMap.set_array(cs)
    fig.colorbar(scalarMap)
    plt.show()
    plt.savefig('goo.png')

if __name__=="__main__":
    scatter3d(x,y,z,cs)
