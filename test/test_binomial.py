import optionpricing

from timeit import timeit


# def test_timeit():
#     bp = optionpricing.BinomialPricing(10000000)

#     print("parallel: ", timeit(bp.sum_parallel, number=1000) / 1000)
#     print("serial: ", timeit(bp.sum_serial, number=1000) / 1000)


def test_mesh_calculate():
    bm = optionpricing.BinomialMesh(1000, 0.0)
    bm.set_initial_condition()

    # print("parallel: ", timeit(bm.calculate_parallel, number=1000) / 1000)
    # print("serial: ", timeit(bm.calculate_serial, number=1000) / 1000)

    n = 10

    bm = optionpricing.BinomialMesh(n, 0.0)
    bm.set_initial_condition()
    print(bm.calculate_serial())

    for i in range(n):
        print(i, 0, bm.get(i, 0))

    bm.print()

    # bm1.print()

    bm = optionpricing.BinomialMesh(10, 0.0)
    bm.set_initial_condition()
    print(bm.calculate_serial())
    bm.print()
