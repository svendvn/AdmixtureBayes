import pandas as pd
import seaborn as sns
import os

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
    for summary in summaries: #the last two are origin and layer which is not what we want to plot of course. 
        if trajectories:
            if trajectories_for_all_temperatures:
                g=sns.FacetGrid(df, row='layer')
                g=summary.make_trajectory(g)
                filename=os.path.join(plot_prefix, summary.name+'_trail_per_temp.png')
                g.savefig(filename)
            if trajectories_for_all_chains:
                g=sns.FacetGrid(df, row='origin')
                g=summary.make_trajectory(g)
                filename=os.path.join(plot_prefix, summary.name+'_trail_per_chain.png')
                g.savefig(filename)
            if not (trajectories_for_all_chains or trajectories_for_all_temperatures):
                g=sns.FacetGrid(main_chain, row='layer')
                filename=os.path.join(plot_prefix, summary.name+'_trail.png')
                g=summary.make_trajectory(g)
                g.savefig(filename)
        if plot_distribution:
            g=sns.FacetGrid(main_chain, row='layer')
            filename=os.path.join(plot_prefix, summary.name+'_distrub.png')
            summary.make_histogram(g)
            if summary.name in prior_distribution:
                summary.add_histogram()
            g.savefig(filename)
                
                
if __name__=='__main__':
    from summary import *
    summaries=[s_variable('posterior'), s_variable('mhr'), s_no_admixes()]
    full_analysis(summaries,
                  trajectories_for_all_temperatures=False,
                  trajectories_for_all_chains=False)
            
    
    
    
    
    