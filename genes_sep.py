import pandas as pd
import numpy as np
import scipy.stats as stats
import statsmodels.stats.multitest as sm_multi
import sys 
import os

def comparison_between_two_groups(grp1, grp2, comparison):
    U_stats, p_value = stats.mannwhitneyu(grp1[comparison].tolist(), grp2[comparison].tolist())
    return p_value

def run():
    df_ori = pd.read_csv(input_dataset, sep='\t')
    df_ori['Origin'] = np.where(df_ori.Sample.str.contains('MB-',case=True), 'METABRIC', 'TCGA')
    proj_name = '_'.join(proj_list)

    for project in proj_list:
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
        # print(grp2)
        # get GI score names
        cols = []
        for column in df:
            cols.append(column)
        GI_list = cols[7:-1]

        pval_list=[]
        means_dict= {'MUT':[],'WT':[]}
        # get pairwise p values
        for i in GI_list:
            p_value = comparison_between_two_groups(grp1, grp2, i)
            curr_mut_mean = grp1[i].mean()
            curr_wt_mean = grp2[i].mean()
            means_dict['MUT'].append(curr_mut_mean)
            means_dict['WT'].append(curr_wt_mean)
            pval_list.append(p_value)

        # correct for false discovery rate
        corrected = sm_multi.multipletests(pval_list, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
        corrected_p_values = corrected[1].tolist()

        # make new df of pvalues
        new_df = {'Gene': wanted_gene, 'Project': project,'Group': "WT_vs_MUT", 'GI_Score': GI_list, 'p-value': pval_list,
                'corrected_p_values': corrected_p_values, 'MUT_means': means_dict['MUT'], 'WT_means': means_dict['WT']}
        pval_df = pd.DataFrame(new_df, columns= ['Gene', 'Project', 'Group', 'GI_Score','p-value','corrected_p_values','MUT_means', 'WT_means'])
        # print(pval_df)

        # add columns
        diff_list = []
        for i, row in pval_df.iterrows():
            pval = row['p-value']
            group1 = row['WT_means']
            group2 = row['MUT_means']

            if (pval <0.05) & (group1 > group2):
                diff_list.append('MUT lower than WT')
            elif (pval <0.05) & (group1 < group2):
                diff_list.append('MUT higher than WT')
            else:
                diff_list.append('N.A.')
        pval_df['Difference'] = diff_list
        pval_df['log10_corrected_pvalues'] = -np.log10(corrected_p_values)
        
        output_pvals = 'temp/'+ proj_name + '_' + wanted_gene +'_pvals.txt'

        # pval_df.to_csv(output_pvals, sep='\t', index=False, index_label=False)
        pval_df.to_csv(output_pvals, mode='a', sep='\t', index=False, index_label=False, header=not os.path.exists(output_pvals))



if __name__ == '__main__':
    input_dataset = sys.argv[1]
    input_geneinfo = sys.argv[2]
    wanted_gene = sys.argv[3]
    proj_list = sys.argv[4:]

    # input_dataset = 'GI_CN.txt'
    # input_geneinfo = 'compiled_somatic_mutation_170522.txt'
    # wanted_gene = 'TP53'
    # proj_list = ['TCGA','METABRIC']

    run()
    