#this file contains functions that work on trees, mostly to extract information.


def get_number_of_admixes(tree):
    tot_length=sum(map(len, tree))
    no_admixs=(tot_length-len(tree)*3)/4
    return no_admixs

def get_all_branch_lengths(tree):
    branch_lengths={}
    for branch in tree:
        times=branch[2::2]
        branch_lengths[branch[0]+branch[1]]=sum(times)
    return branch_lengths
    