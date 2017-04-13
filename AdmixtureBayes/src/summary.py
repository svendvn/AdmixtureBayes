import matplotlib.pyplot as plt
from Rtree_operations import get_number_of_admixes, get_all_branch_lengths
from tree_statistics import unique_identifier
from data_wrangling_functions import values_to_numbers, count_strings, count_strings2, thin_out_nans
from numpy import isfinite, array, cumsum, linspace


class Summary(object):
       
    def __init__(self, name, pandable=True, output='double'):
        self.name=name
        self.pandable=pandable
        self.output=output
    
    def __call__(self, **kwargs):
        pass
    
    def pretty_print(self, output):
        return str(output)
    
    def summary_of_phylogeny(self, tree):
        return None
    
    def make_trajectory(self, index, x, **kwargs):
        if isinstance(x[0], float):
            x, index=thin_out_nans(x,index)
            plt.plot(index,x,**kwargs)
        else:
            numbers=values_to_numbers(x)
            plt.plot(index, numbers, **kwargs)
        plt.title(self.name)
    
    

            
    def make_histogram(self, x, a=None, **kwargs):
        if isinstance(x[0], float):
            x,_=thin_out_nans(x,x)
            _,bins,_= plt.hist(x, fc=(1, 0, 0, 0.5), normed=True, **kwargs)
            if a is not None:
                a,_=thin_out_nans(a,a)
                if 'bins' in kwargs:
                    plt.hist(a,fc=(0, 1, 0, 0.5), normed=True,  **kwargs)
                else:
                    plt.hist(a,fc=(0, 1, 0, 0.5), normed=True, bins=bins, **kwargs)
        else:
            if a is None:
                labels, counts1 = count_strings(x)
                counts2=None
            else:
                labels, counts1, counts2 = count_strings2(x,a)
            print 'labels', labels
            plt.bar(range(len(labels)), counts1, width=1.0, alpha=0.5, color='r', label='MCMC')
            if counts2 is not None:
                plt.bar(range(len(labels)), counts2, width=1.0, alpha=0.5, color='g', label='MCMC')
        plt.title(self.name)
    
class s_no_admixes(Summary):
    
    def __init__(self):
        super(s_no_admixes,self).__init__('no_admixes', output='integer')

    def __call__(self, **kwargs):
        old_tree=kwargs['tree']
        return get_number_of_admixes(old_tree)
    
    def summary_of_phylogeny(self, tree):
        return get_number_of_admixes(tree)
    
    def make_histogram(self, x,a):
        maxval=max(x)
        if a is not None:
            maxval=max(maxval, max(a))
        return Summary.make_histogram(self, x, a, bins=range(maxval+1))
    

class s_total_branch_length(Summary):

    def __init__(self):
        super(s_total_branch_length,self).__init__('total_branch_length')

    def __call__(self, **kwargs):
        tree=kwargs['tree']
        return sum(get_all_branch_lengths(tree))
    
    def summary_of_phylogeny(self, tree):
        return sum(get_all_branch_lengths(tree))
    
class s_average_branch_length(Summary):
    
    def __init__(self):
        super(s_average_branch_length,self).__init__('average_branch_length')

    def __call__(self, **kwargs):
        tree=kwargs['tree']
        all_branch=get_all_branch_lengths(tree)
        return sum(all_branch)/len(all_branch)
    
    def summary_of_phylogeny(self, tree):
        all_branch=get_all_branch_lengths(tree)
        return sum(all_branch)/len(all_branch)
    
class s_variable(Summary):
    
    def __init__(self, variable, pandable=True, output='double'):
        super(s_variable, self).__init__(variable, pandable, output)

    def __call__(self, **kwargs):
        if self.name not in kwargs:
            return None
        return kwargs[self.name]
    
    
class s_tree_identifier(Summary):
    
    def __init__(self):
        super(s_tree_identifier,self).__init__('tree_identifier', output='string')
        
    def __call__(self, **kwargs):
        old_tree=kwargs['tree']
        return unique_identifier(old_tree)
    
    def summary_of_phylogeny(self, tree):
        return unique_identifier(tree)
    
    def make_histogram(self, x, a=None, **kwargs):
        if a is None:
            labels, counts1 = count_strings(x)
            counts2=None
        else:
            labels, counts1, counts2 = count_strings2(x,a)
            labels, counts1, counts2 = zip(*[(x,y,z) for (x,y,z) in sorted(zip(labels,counts1,counts2), key=lambda pair: pair[1], reverse=True)])
            print 'labels', labels
            xaxis=range(0,20,2)+range(20,30,1)+list(linspace(30,40, len(labels)-20))
            yaxis=counts1
            zaxis=[2]*10+[1]*10+[10.0/(len(labels)-20)]*(len(labels)-20)
            print len(xaxis), xaxis
            print len(yaxis), yaxis
            print len(zaxis), zaxis
            print labels[:20]
            plt.bar(xaxis, counts1, width=zaxis, alpha=0.5, color='r', label='MCMC')
            if counts2 is not None:
                plt.bar(xaxis, counts2, width=zaxis, alpha=0.5, color='g', label='MCMC')
        plt.title(self.name)
    

class s_tree_identifier_new_tree(s_tree_identifier):
    
    def __init__(self):
        super(s_tree_identifier,self).__init__('tree_identifier_new_tree', output='string')
        
    def __call__(self, **kwargs):
        tree=kwargs['proposed_tree']
        return unique_identifier(tree)
        
        
 