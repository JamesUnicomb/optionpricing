import optionpricing
import numpy as np
from timeit import timeit



bm2 = optionpricing.BinomialMesh1D(48, 0.0)
bm2.set_initial_condition()
bm2.calculate_serial()
