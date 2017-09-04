from subprocess import call
import os


class stop_criteria(object):
    
    def __init__(self, frequency=20000, summaries=['no_admixes','average_branch_length','add'], outfile='tmp_stop_criteria.txt'):
        self.counter=0
        self.frequency=frequency
        self.summaries=summaries
        self.outfile=outfile
        
        
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
        command=['Rscript',os.path.join(dir, 'ESS.R'), filename, '0.5', self.outfile]+self.summaries
        print command
        call(command)
        return self.check_outfile()
    
    def check_outfile(self):
        with open(self.outfile, 'r') as f:
            f.readline()
            for lin in f.readlines():
                ess=float(lin.split()[-1])
                name=lin.split()[1]
                print name, ess
                if ess<200:
                    return False
        return True
    
if __name__=='__main__':
    sc=stop_criteria()
    sc2=stop_criteria(summaries=['posterior'])
    print 'sc(10000)', sc(10000, filename='result_mc3_localcopy.csv')
    print 'sc2(10000)', sc2(10000, filename='result_mc3_localcopy.csv')
    print 'sc(100000)', sc(100000, filename='result_mc3_localcopy.csv')
    print 'sc2(100000)', sc2(100000, filename='result_mc3_localcopy.csv')