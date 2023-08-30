import optionpricing
import numpy as np
from timeit import timeit


bm1 = optionpricing.BinomialMesh1D(48, 0.0)
bm1.set_initial_condition()
print(bm1.calculate_parallel())


# bm2 = optionpricing.BinomialMesh1D(24, 0.0)
# bm2.set_initial_condition()
# print(bm2.calculate_serial())
