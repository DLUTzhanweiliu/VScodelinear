import numpy as np

p1 = np.array([65, 150, 78.3])
p2 = np.array([70, 150, 84.0])
p3 = np.array([65, 300, 113.1])

# These two vectors are in the plane
v1 = p3 - p1
v2 = p2 - p1

# the cross product is a vector normal to the plane
cp = np.cross(v1, v2)
a, b, c = cp

# This evaluates h * x3 + w * y3 + g * z3 which equals d
d = np.dot(cp, p3)
e = c
a, b, c, d = -a/e, -b/e , c/e, d/e

print('The equation is g = [{0} {1}][h, w]^T + {2}'.format(a, b, d))