from dolfin import * 
from dolfin_utils import meshconvert

# Convert to XML
meshconvert.convert2xml(shape_0, shape_0, iformat=mesh)
# Order mesh
os.system("dolfin-order %s" % shape_0)
