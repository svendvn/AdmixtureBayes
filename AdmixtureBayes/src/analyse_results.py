import pandas as pd
import seaborn as sns
import os
import matplotlib.pyplot as plt
from subprocess import call

def save_to_csv(list_of_tuples, summaries):
    df=pd.DataFrame.from_items([('iteration',list_of_tuples[0])]+[(summ_object.name,summ_col) for summ_object,summ_col in zip(summaries,list_of_tuples[1:])])
    df['origin']=1
    df['layer']=1
    df.to_csv('results.csv')
    
    
def full_analysis(summaries,
                  csv_file='results.csv', 
                  trajectories=True, 
                  trajectories_for_all_temperatures=True,
                  trajectories_for_all_chains=True,
                  burn_in=0.0,
                  plot_distribution=True,
                  prior_distribution={},
                  plot_prefix='plots'):
    df=pd.DataFrame.from_csv(csv_file)
    main_chain=df.query('layer == 1')
    print main_chain
    if plot_distribution and prior_distribution:
        df_priors=pd.DataFrame.from_dict(prior_distribution)
    for summary in summaries: #the last two are origin and layer which is not what we want to plot of course. 
        print summary.name
        if trajectories:
            if trajectories_for_all_temperatures:
                for layer,data in df.groupby('layer'):
                    filename=os.path.join(plot_prefix, summary.name+'_trail_layer'+str(layer)+'.png')
                    fig= plt.figure()
                    summary.make_trajectory(data[summary.name])
                    fig.savefig(filename, bbox_inches='tight')
            if trajectories_for_all_chains:
                for chain,data in df.groupby('origin'):
                    filename=os.path.join(plot_prefix, summary.name+'_trail_chain'+str(chain)+'.png')
                    fig= plt.figure()
                    summary.make_trajectory(data[summary.name])
                    fig.savefig(filename, bbox_inches='tight')
            if not (trajectories_for_all_chains or trajectories_for_all_temperatures):
                filename=os.path.join(plot_prefix, summary.name+'_trail_layer1.png')
                fig= plt.figure()
                summary.make_trajectory(main_chain[summary.name])
                fig.savefig(filename, bbox_inches='tight')
        if plot_distribution:
            filename=os.path.join(plot_prefix, summary.name+'_distrub.png')
            fig= plt.figure()
            if summary.name in prior_distribution:
                summary.make_histogram(main_chain[summary.name], a=prior_distribution[summary.name])
            else:
                summary.make_histogram(main_chain[summary.name], a=None)
            fig.savefig(filename, bbox_inches='tight')
    #call_notebook()
    
def call_notebook():
    ## DOESNT WORK IN THE MAC 
    dir_path = os.path.dirname(os.path.realpath(__file__))
    cmd=['Rscript', dir_path+os.path.sep+'order_report.R']
    print cmd
    call(cmd)
    
                
                
if __name__=='__main__':

    import summary
    summaries=[summary.s_variable('posterior'), 
               summary.s_variable('mhr'), 
               summary.s_no_admixes(), 
               summary.s_tree_identifier(),
               summary.s_average_branch_length(),
               summary.s_total_branch_length(),
               summary.s_tree_identifier_new_tree()]+[summary.s_variable(s) for s in ['backward_choices',
                                                                                      'backward_density',
                                                                                      'forward_density',
                                                                                      'forward_choices',
                                                                                      'proposal_type',
                                                                                      'prior',
                                                                                      'branch_prior',
                                                                                      'no_admix_prior',
                                                                                      'top_prior']]
    from generate_prior_trees import get_distribution_under_prior
    prior_distribution=get_distribution_under_prior(leaves=4, sim_length=1000, list_of_summaries=[summaries[2]])
    print prior_distribution
    full_analysis(summaries,
                  trajectories_for_all_temperatures=False,
                  trajectories_for_all_chains=False,
                  prior_distribution=prior_distribution)
            
    
    
    
    
    