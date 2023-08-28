import optionpricing

from timeit import timeit


def test_timeit():
    bp = optionpricing.BinomialPricing(10000000)

    print("parallel: ", timeit(bp.sum_parallel, number=1000) / 1000)
    print("serial: ", timeit(bp.sum_serial, number=1000) / 1000)


def test_mesh_calculate():
    bm = optionpricing.BinomialMesh(1000, 0.0)

    bm.set_initial_condition()

    print("parallel: ", timeit(bm.calculate_parallel1, number=1000) / 1000)
    print("parallel: ", timeit(bm.calculate_parallel2, number=1000) / 1000)
    print("serial: ", timeit(bm.calculate_serial, number=1000) / 1000)
