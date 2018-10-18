from subprocess import call
import os


class stop_criteria(object):
    
    TOPOLOGICAL_MIN=200
    CONTINUOUS_MIN=200
    BURN_IN=0.5
    
    def __init__(self, 
                        frequency=20000, 
                        summaries=['no_admixes','average_branch_length','add','descendant_sets'], 
                        outfile='tmp_stop_criteria.txt', 
                        topological=False, 
                        verbose_level='normal'):
        self.counter=0
        self.frequency=frequency
        self.summaries=summaries
        self.non_topological_summaries=len(summaries)
        self.outfile=outfile
        self.verbose_level=verbose_level
        if topological:
            self.summaries.extend(['Zero_Ntree','random_Ntree','mode_Ntree'])
        
        
    def __call__(self, cum_iterations, filename):
        if self.check_yet(cum_iterations):
            return self.stop_yet(filename)
        return False
        
    def check_yet(self, cum_iterations):
        if int(cum_iterations/self.frequency)>self.counter:
            self.counter+=1
            return True
        return False
    
    def stop_yet(self, filename):
        dir = os.path.dirname(__file__)
        command=['Rscript',os.path.join(dir, 'ESS.R'), filename, str(stop_criteria.BURN_IN), self.outfile,  str(self.verbose_level)]+self.summaries
        if self.verbose_level=='normal':
            print command
        call(command)
        return self.check_outfile()
    
    def check_outfile(self):
        with open(self.outfile, 'r') as f:
            f.readline()
            for n,lin in enumerate(f.readlines()):
                ess=float(lin.split()[-1])
                name=lin.split()[0]
                if self.verbose_level=='normal':
                    print n, name, ess
                if n<self.non_topological_summaries:
                    if ess<stop_criteria.CONTINUOUS_MIN:
                        return False
                else:
                    if ess<stop_criteria.TOPOLOGICAL_MIN:
                        return False
        return True
    
if __name__=='__main__':
    sc=stop_criteria(topological=False)
    sc2=stop_criteria(summaries=['posterior','descendant_sets'])
    print 'sc(10000)', sc2(100000, filename='result_mc3.csv')
    #print 'sc2(10000)', sc2(10000, filename='result_mc3_localcopy.csv')
    #print 'sc(100000)', sc(100000, filename='result_mc3_localcopy.csv')
    #print 'sc2(100000)', sc2(100000, filename='result_mc3_localcopy.csv')
