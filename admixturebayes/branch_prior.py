from scipy.stats import gamma, expon
from math import log, exp
from numpy import roots, iscomplex, real
from random import random

class double_expon(object):

    @staticmethod
    def dispersion_to_scale(dispersion):
        r=roots([-1,1,3,-2-1.0/dispersion,1])
        for r_val in r:
            if r_val>0 and r_val<1:
                assert not iscomplex(r_val)
                return real(r_val)


    @staticmethod
    def logpdf(x, scale):
        expon.pdf(x=x,)

    @staticmethod
    def rvs(scale):
        if random()<1.0-scale:
            vscale=(1.0-scale**2)/(1.0-scale)
            return expon.rvs()


if __name__=='__main__':
    print double_expon.dispersion_to_lambda(0.04)