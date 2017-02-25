from pathos.multiprocessing import Pool


def trier(a):
    return 1,'2',[(5,6),(7,8)]

p=Pool(4)

print zip(*p.map(trier, range(4)))