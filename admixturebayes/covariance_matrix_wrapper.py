from numpy import zeros, diag

import _covariance_matrix

class Covariance_Matrix2():
    
    def __init__(self, nodes_to_index):
        self.ni=nodes_to_index
        self.covmat=zeros((len(nodes_to_index), len(nodes_to_index)))
    
    def get_indices(self, nodes):
        return [self.ni[n] for n in nodes]
    
    def update(self, branch_length, population):
        self.covmat=_covariance_matrix.update_matrix(self.covmat, branch_length, population.weights, self.get_indices(population.members))
        
    def get_matrix(self):
    	return self.covmat+self.covmat.T-diag(self.covmat.diagonal())