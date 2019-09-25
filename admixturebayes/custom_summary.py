from Rtree_operations import node_is_admixture, get_parents, get_sibling_on_parent_side



def first_admixture_proportion_custom_summary(Rtree, **kwargs):
    '''
    This is an example of a custom function. It returns the admixture proportion of the immediate ancestors of population s1 and s3.
    The function is being read by downstream_analysis_parser because it ends with '_custom_summary'.
    :param Rtree: A tree
    :param kwargs: all the parameters passed to all summary functions that are not used by this function
    :return: a float. the admixture proportion.
    '''
    #print 'entering first admixture proportion'
    for key, node in Rtree.items():
        if node_is_admixture(node): #checking for the pattern s1/ nX \ a1 \s2 /a1 / nY \ s3
            p1,p2=get_parents(node)
            p1_other_children=get_sibling_on_parent_side(Rtree, p1, key)
            p2_other_children=get_sibling_on_parent_side(Rtree, p2, key)
            if 'MA1' in p1_other_children and 'LBK' in p2_other_children:
                return node[2]
            if 'LBK' in p1_other_children and 'MA1' in p2_other_children:
                return 1.0-node[2]
    else:
        return '.'






if __name__ == '__main__':

    ##The tree is written in the format
    ##  {node_name: [ parent_node_1_name,
                    # parent_node_2_name,
                    # admixture_proportion_between_parent_1 and parent_2,
                    # branch_length_of_branch_to_parent_1,
                    # branch_length_of_branch_to_parent_2,
                    # child_node_1_name,
                    # child_node_2_name]
    #If a node does not have 2 parents or 2 children the corresponding entries will be None.
    Rtree={'s1':['n1',None, None, 1, None, None, None],
            's2':['a1', None, None,0.5, None, None,None],
            's3':['n2', None, None, 3, None, None, None],
            's4':['n3', None, None, 2, None, None, None],
            'n1':['n4', None, None, 4,None,'s1','a2'],
            'a1':['a2', 'n2',0.3, 0.2,0.3,'s2',None],
            'n2':['n4', None,None, 2.5, None, 's3','a1'],
            'n3':['r', None, None, 1.8, None, 's4','a2'],
            'a2':['n3', 'n1',0.5,0.2,0.18,'a1',None],
            'n4':['r', None, None, 0.05, None, 'n2','n1']}
    add=0
    inputs={'Rtree':Rtree, 'add':add}
    print first_admixture_proportion_custom_summary(**inputs)















##Wrap functions
def wrap_function(func, func_name):
    def class_initialization_wrap():
        def new_function(**kwargs):
            return_val= func(**kwargs)
            return {func_name:return_val},False
        return new_function
    return class_initialization_wrap


##wraps all functions ending in _custom_summary to return to downstream analysis parser for making custom functions.
def all_custom_summaries():
    commands = {}
    for name, value in globals().items():
        if callable(value) and name.endswith("_custom_summary"):
            command_name = name.rsplit("_custom_summary", 1)[0]
            commands[command_name] = wrap_function(value, command_name)
    return commands


