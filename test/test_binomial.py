import optionpricing
import numpy as np
from timeit import timeit


def test_print():
    # bm1 = optionpricing.BinomialMesh1D(24, 0.0)
    # bm1.set_initial_condition()
    # # print(bm1.calculate_parallel())
    # bm1.calculate_parallel()

    bm2 = optionpricing.BinomialMesh1D(24, 0.0)
    bm2.set_initial_condition()
    # print(bm2.calculate_serial())
    bm2.calculate_serial()

    # bm2 = optionpricing.BinomialMesh1D(12, 0.0)
    # bm2.set_initial_condition()
    # bm2.calculate_serial()
    # out2 = bm2.print()

    # assert out1 == out2


# def test_print():
#     bm1 = optionpricing.BinomialMesh(12, 0.0)
#     bm1.set_initial_condition()
#     bm1.calculate_parallel()
#     out1 = bm1.print()

#     bm2 = optionpricing.BinomialMesh(12, 0.0)
#     bm2.set_initial_condition()
#     bm2.calculate_serial()
#     out2 = bm2.print()

#     assert out1 == out2


def test_mesh_calculate():
    N = 600
    num = 1

    # bmp = optionpricing.BinomialMesh(N, 0.0)
    # bmp.set_initial_condition()
    # ptime = timeit(bmp.calculate_parallel, number=num) / num
    # print("parallel: ", ptime)

    # bms = optionpricing.BinomialMesh(N, 0.0)
    # bms.set_initial_condition()
    # stime = timeit(bms.calculate_serial, number=num) / num
    # print("serial: ", stime)

    bmp1d = optionpricing.BinomialMesh1D(N, 0.0)
    bmp1d.set_initial_condition()
    v1 = bmp1d.calculate_parallel()
    ptime = timeit(bmp1d.calculate_parallel, number=num) / num
    print("parallel: ", ptime)

    bms1d = optionpricing.BinomialMesh1D(N, 0.0)
    bms1d.set_initial_condition()
    v2 = bms1d.calculate_serial()

    stime = timeit(bms1d.calculate_serial, number=num) / num
    print("serial: ", stime)

    # print("speedup: ", stime / ptime, "x")

    print(v1, v2)

