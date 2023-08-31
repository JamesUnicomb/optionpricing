import optionpricing
import numpy as np
from timeit import timeit


bm1 = optionpricing.BinomialMesh1D(48, 0.0)
bm1.set_initial_condition()
bm1.calculate_parallel()
