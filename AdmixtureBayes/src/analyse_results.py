import pandas as pd

def save_to_csv(list_of_tuples, summaries):
    df=pd.DataFrame.from_items(((summ_object.name,summ_col) for summ_object,summ_col in zip(summaries,list_of_tuples)))
    df['origin']=1
    df['layer']=1
    df.to_csv('results.csv')
    
def full_analysis(csv_file='results.csv', 
                  trajectories=True, 
                  trajectories_for_all_temperatures=True,
                  burn_in=0.0,
                  plot_distribution=True,
                  prior_distribution={},
                  plot_prefix='plots/'):
    
    
    
    
    