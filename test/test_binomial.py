import optionpricing
import numpy as np
from timeit import timeit


def test_print():
    bm1 = optionpricing.BinomialMesh(12, 0.0)
    bm1.set_initial_condition()
    bm1.calculate_parallel()
    out1 = bm1.print()

    bm2 = optionpricing.BinomialMesh(12, 0.0)
    bm2.set_initial_condition()
    bm2.calculate_serial()
    out2 = bm2.print()

    assert out1 == out2


def test_mesh_calculate():
    N = 1500

    bmp = optionpricing.BinomialMesh(N, 0.0)
    bmp.set_initial_condition()

    ptime = timeit(bmp.calculate_parallel, number=1000) / 1000
    print("parallel: ", ptime)

    bms = optionpricing.BinomialMesh(N, 0.0)
    bms.set_initial_condition()
    stime = timeit(bms.calculate_serial, number=1000) / 1000
    print("serial: ", stime)
    print("speedup: ", stime / ptime, "x")

    bms = optionpricing.BinomialMesh(N, 0.0)
    bms.set_initial_condition()
    bms.calculate_serial()

    bmp = optionpricing.BinomialMesh(N, 0.0)
    bmp.set_initial_condition()
    bmp.calculate_parallel()

    for i in range(N):
        assert np.isclose(bmp.get(i, 0), bms.get(i, 0))
