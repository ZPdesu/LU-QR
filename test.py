from numpy import *
from scipy.linalg import *

U = zeros((4, 4), float)
print type(U)

a = array([[0.6,0.64,-0.48],[-20.0, 27.0,11.0],[-20.0, 27.0,11.0]])
b = array([[-20.0, 27.0],[-20.0, 27.0],[-20.0, 27.0]])
print float(0.6) * float(-20) + float(0.64)*float(27) + float(-0.48) *float(11)
c = around(a,4)

c[0,0] = 5
print a
print c