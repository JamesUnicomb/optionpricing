import optionpricing
import numpy as np
from timeit import timeit


def test_mesh_calculate():
    bmp = optionpricing.BinomialMesh(1000, 0.0)
    bmp.set_initial_condition()

    ptime = timeit(bmp.calculate_parallel, number=1000) / 1000
    print("parallel: ", ptime)

    bms = optionpricing.BinomialMesh(1000, 0.0)
    bms.set_initial_condition()
    stime = timeit(bms.calculate_serial, number=1000) / 1000
    print("serial: ", stime)
    print("speedup: ", stime / ptime, "x")

    bms = optionpricing.BinomialMesh(1200, 0.0)
    bms.set_initial_condition()
    bms.calculate_serial()

    bmp = optionpricing.BinomialMesh(1200, 0.0)
    bmp.set_initial_condition()
    bmp.calculate_parallel()

    for i in range(10):
        for j in range(i + 1):
            assert np.isclose(bmp.get(i, j), bms.get(i, j))

    # bm1 = optionpricing.BinomialMesh(12, 0.0)
    # bm1.set_initial_condition()
    # bm1.calculate_parallel()
    # bm1.print()

    # bm2 = optionpricing.BinomialMesh(12, 0.0)
    # bm2.set_initial_condition()
    # bm2.calculate_serial()
    # bm2.print()
