from scipy.stats import expon, geom
from Rtree_operations import get_all_branch_lengths, get_all_admixture_proportions

def prior(tree, p=0.5):
    admixtures=get_all_admixture_proportions(tree)
    if not all(prop>=0 and prop <=1 for prop in admixtures):
        return -float('inf')
    branches=get_all_branch_lengths(tree)
    logsum=sum(map(expon.logpdf, branches))
    logsum+= geom.logpmf(len(admixtures), 1.0-p)
    return logsum
    