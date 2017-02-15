#this file contains functions that work on trees, mostly to extract information.


def get_number_of_admixes(tree):
    tot_length=sum(map(len, tree))
    no_admixs=(tot_length-len(tree)*3)/4
    return no_admixs

def get_number_of_admixes_on_branch(branch):
    return (len(branch)-3)/2

def get_all_branch_lengths(tree):
    branch_lengths={}
    for branch in tree:
        times=branch[2::2]
        branch_lengths[branch[0]+branch[1]]=sum(times)
    return branch_lengths

def make_flat_list_no_admix(size=2):
    step_size=1.0/size
    res=[["s1s2","s1",step_size], ["s1s2","s2",step_size]]
    rooter="s1s2"
    for i in range(3,size+1):
        old_rooter=rooter
        si="s"+str(i)
        rooter+=si
        res.append([rooter,si, step_size*(i-1)])
        res.append([rooter,old_rooter,step_size])
    res[-2][0]="r"
    res[-1][0]="r"
    return res

def get_number_of_leaves(tree):
    pass