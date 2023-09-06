import optionpricing
import numpy as np
from timeit import timeit




for i in range(1000):

    s = int(np.power(2, np.random.randint(8, 14)))
    k = 2*np.random.randint(2, 8)
    
    bm1 = optionpricing.BinomialMesh1D(s, 0.0)
    bm1.set_initial_condition()

    bm2 = optionpricing.BinomialMesh1D(s, 0.0)
    bm2.set_initial_condition()

    v1 = bm1.calculate_parallel()
    v2 = bm2.calculate_serial()

    print(s,k,np.round(v1, 3), np.round(v2, 3))

    if not np.isclose(v1,v2, atol=1e-3):
        raise Exception
   