
from shapes_utils import *
from meshes_utils import *


def generateim():
    n_pts          = 5
    n_sampling_pts = 10
    radius         = 0.5*np.ones([n_pts])
    edgy           = 0.5*np.ones([n_pts])
    plot_pts       = True
    filename       = 'shape_{}.csv'.format(i)
    cylinder       = False

    # Generate and mesh shape
    shape = Shape(filename,None,n_pts,n_sampling_pts,radius,edgy)
    #shape = Shape()
    shape.read_csv(filename)
    shape.generate()
    #shape.generate(cylinder=cylinder)
    
    shape.generate_image(plot_pts=plot_pts)
    shape.write_csv()

for i in range(1,186):
    if __name__ == "__main__":
        generateim()
