import numpy as np

import optionpricing

N = 9
bm1 = optionpricing.BinomialMesh(N, 0.0)
bm1.set_initial_condition()
bm1.calculate_parallel()
bm1.print()

bm2 = optionpricing.BinomialMesh(N, 0.0)
bm2.set_initial_condition()
bm2.calculate_serial()
bm2.print()