import numpy as np

import optionpricing

N = 48
bm1 = optionpricing.BinomialMesh(N, 0.0)
bm1.set_initial_condition()
print(bm1.calculate_parallel())
bm1.print()

bm2 = optionpricing.BinomialMesh(N, 0.0)
bm2.set_initial_condition()
print(bm2.calculate_serial())
bm2.print()
