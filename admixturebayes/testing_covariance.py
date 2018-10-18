import numpy as np

class Covariance_Matrix():
    
    def __init__(self, nodes_to_index):
        self.nodes_to_index=nodes_to_index
        self.covmat=np.zeros((len(nodes_to_index), len(nodes_to_index)))
    
    def update(self, branch_length, members, weights):
        indices=[self.nodes_to_index[n] for n in members]
        self.covmat[np.ix_(indices,indices)]+=branch_length*np.outer(weights, weights)


from random import random
def call_update_realistically():
    weights=[random()]
    members=["s0"]
    covmat=Covariance_Matrix({str("s"+str(i)):i for i in range(40)})
    for i in range(1,40):
        covmat.update(random(), members, weights)
        new_member="s"+str(i)
        new_weight=random()
        covmat.update(random(), [new_member],[new_weight])
        members.append(new_member)
        weights.append(new_weight)
        
import cProfile

cProfile.run("[call_update_realistically() for _ in range(20)]")
    




