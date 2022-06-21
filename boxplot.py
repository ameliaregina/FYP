import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns
import sys 

def run():
    df_ori = pd.read_csv(input_dataset, sep='\t')
    df_ori['Origin'] = np.where(df_ori.Sample.str.contains('MB-',case=True), 'METABRIC', 'TCGA')
    # proj_name = '_'.join(proj_list)
    # for project in proj_list:
    if project == 'METABRIC':
        df = df_ori[df_ori['Origin'] == 'METABRIC']
    if project == 'TCGA':
        df = df_ori[df_ori['Origin'] == 'TCGA']

    df2 = pd.read_csv(input_geneinfo, sep='\t')
    merge_df = pd.merge(df,df2)

    # driver_genes = ['PIK3CA', 'TP53', 'GATA3', 'KMT2C', 'MAP3K1', 'CDH1', 'TBX3', 
    #                 'ARID1A', 'CBFB', 'PTEN', 'NCOR1', 'AKT1', 'RUNX1', 'NF1', 'MAP2K4']

    grp1 = merge_df[merge_df['Gene']== wanted_gene]
    grp1.drop_duplicates(subset ='Sample',keep = 'first', inplace = True)
    # values not in grp1
    grp2= pd.merge(df,grp1, indicator=True, how='outer', on='Sample').query('_merge=="left_only"').drop('_merge', axis=1)
    grp2.columns = grp2.columns.str.strip('_x')

    grp1_new = grp1[[GI_score]].assign(Location=1)
    grp2_new = grp2[[GI_score]].assign(Location=2)


    cdf = pd.concat([grp1_new, grp2_new])    
    mdf = pd.melt(cdf, id_vars=['Location'], var_name=['Letter']) 

    fig1, ax = plt.subplots(figsize=(4,6))
    ax = sns.boxplot(x="Location", y="value", data=mdf, showfliers=False, showmeans=True, palette='Set2', linewidth=1)  
    ax.set_ylabel(GI_score)  
    ax.set_xlabel('Breast cancer')
    ax.set_title(wanted_gene)
    ax.set_xticks([0,1])
    ax.set_xticklabels(['MUT','WT'])

    plt.tight_layout()
    
    output_png = 'res2/'+ project + '_' + wanted_gene + '_'+ GI_score +'.png'
    plt.savefig(output_png, format='png', bbox_inches='tight')

if __name__ == '__main__':
    input_dataset = sys.argv[1]
    input_geneinfo = sys.argv[2]
    wanted_gene = sys.argv[3]
    GI_score = sys.argv[4]
    project = sys.argv[5]
    
    run()
