import optionpricing

from timeit import timeit

def test_timeit():
    bp = optionpricing.BinomialPricing(10000000)

    print("serial: ", timeit(bp.sum_serial, number=1000)/1000)
    print("parallel: ", timeit(bp.sum_parallel, number=1000)/1000)

def test_mesh_init():
    bm = optionpricing.BinomialMesh(10, 0.0)

    print(bm.get(1,1))

    bm.print()