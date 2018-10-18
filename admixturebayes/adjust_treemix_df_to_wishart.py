from scipy.stats import wishart, norm
from Rtree_to_covariance_matrix import make_covariance

def adjust_treemix_df(wishart_df, starting_tree):
    cov=make_covariance(starting_tree)
    lmax=wishart.logpdf(wishart, scale=wishart/wishart_df, df=wishart_df)
    