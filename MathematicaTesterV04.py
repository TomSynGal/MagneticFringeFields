from MathematicaFunctionsV04 import u, v, zeta, hj, G, f, g, bj, cj, bu, bv, bx, by
import numpy as np


# Updates from previous talk on Thursday

# Now making code in a newer/neater file with the functions in a seperate folder.

# Function G now has an if statement to allow values that are invalid to output
# the number 0.5. This completes the graphs created.

# Changed how the bj values are called and created an initial bj which can be
# acted upon by a single function bj to determine the element bjInitiate[[0]].
# All changes to bj for each multipole type are now done in function and are all
# automated correctly including the correct display of complex numbers.

uTest = u(1,1)
# Complex 128, size 1, complex number.
print('The value of u test, ', uTest)
# This has passed the test.


vTest = v(1,1)
# Complex 128, size 1, complex number.
print('The value of v test, ', vTest)
# This has passed the test.


zetaTest = zeta(1)
# Float 64, size 1.
print('The value of zeta test, ', zetaTest)
# This has passed the test.


GTest = G(0,1)
# MPF object, size 1, complex number but.. single float 0.5 when first parameter xi == 0.
print('The value of G test, ', GTest)
# This has passed the test.


fTest = f(7,3)
# MPF object, size 1, complex number but.. single float 0 when first parameter eta == 0.
print('The value of f test, ', fTest)
# This has passed the test.


gTest = g(7,3)
# MPF object, size 1, complex number but.. single float 0 when first parameter eta == 0.
print('The value of g test, ', gTest)
# This has passed the test.


bjArray = np.array([0, 0.95, 0.90, 0.85], dtype = np.complex_)
# Float 64 Numpy array. Specified np.complex to show the values as complex numbers.


bjTest = bj(bjArray, 1)
# Complex 128, size (2,).
print('The value of bj test, ', bjTest)
# This has passed the test.


hjTest = hj(1, 1, bjArray, 1, 1)
# Complex 128, size 1.
print('The value of hj test, ', hjTest)
# This has passed the test.


cjTest = cj(bjArray, 1)
# Complex 128, size (2,).
print('The value of cj test, ', cjTest)
# This has passed the test.


buTest = bu(1, 1, 1, bjArray, 1)
# Python MPC object, size 1.
print('The value of bu test, ', buTest)
# This has passed the test.


bvTest = bv(1, 1, 1, bjArray, 1)
# Python MPC object, size 1.
print('The value of bv test, ', bvTest)
# This has passed the test.


bxTest = bx(1, 1, 1, bjArray, 1)
# Python MPC object, size 1.
print('The value of bx test, ', bxTest)
# This has passed the test.


byTest = by(1, 1, 1, bjArray, 1)
# Python MPC object, size 1.
print('The value of by test, ', byTest)
# This has passed the test.