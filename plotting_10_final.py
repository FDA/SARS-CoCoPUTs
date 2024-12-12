import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import warnings
import sys
import numpy as np
import seaborn as sns
import codonbias as cb
from statsmodels.graphics.tsaplots import plot_acf
from statsmodels.tsa.stattools import adfuller
from collections import Counter
import re
from concurrent.futures import ProcessPoolExecutor
import subprocess
import os
import dask.dataframe as dd
from dask.distributed import Client, LocalCluster, wait
from matplotlib.colors import ListedColormap
from dask import config as cfg
import matplotlib.dates as mdates

#this will prevent heartbeat errors from dask
cfg.set({'distributed.scheduler.worker-ttl': None})

warnings.filterwarnings("ignore")

def vait_pie_charts(ss_df, vait_pie_chart_path):
    vait_spike = 'TCTGCTTTACTAATGTCTATGCAGATT'
    spike_vait_regex_string = ''
    for i in vait_spike:
        spike_vait_regex_string += i
        spike_vait_regex_string += '-*'
    spike_vait_regex_string = spike_vait_regex_string[0:-2]
    # wt_genome_sequence = ss_df.loc['referenceB','sequence']
    wt_genome_sequence = list(ss_df[ss_df['terminal_lineage'].str.contains('referenceB')]['sequence'])[0]

    spike_search = re.search(spike_vait_regex_string, wt_genome_sequence)

    ss_df['spike_vait_forgi_struct'] = ss_df['msa_forgi_output'].str[spike_search.start():spike_search.end()].str.replace("-","")
    
    spike_vait_value_counts = dict(Counter(ss_df['spike_vait_forgi_struct']))
    largest_category = max(spike_vait_value_counts, key=spike_vait_value_counts.get)
    
    modified_spike_vait_value_counts = {"Wild Type": spike_vait_value_counts[largest_category]}
    modified_spike_vait_value_counts["Other"] = sum(value for key, value in spike_vait_value_counts.items() if key != largest_category)

    spikevaitdf = pd.DataFrame(list(modified_spike_vait_value_counts.items()), columns=['Sequence', 'Value'])

    # create data
    names = list(spikevaitdf['Sequence'])
    size = list(spikevaitdf['Value'])

    # Create a circle at the center of the plot
    my_circle = plt.Circle( (0,0), 0.7, color='white')

    # Custom wedges
    patches, texts, autotexts = plt.pie(size, labels=names,autopct='%1.1f%%', textprops={'fontsize': 14})
    for autotext in autotexts:
        autotext.set_color('black')
        autotext.set_size(16)
        # Move the text towards the center of the pie
        autotext.set_position((autotext.get_position()[0] * 0.8, autotext.get_position()[1] * 0.8))

    #finish the plot
    p = plt.gcf()
    p.gca().add_artist(my_circle)
    plt.title('Spike VAIT region secondary structure composition\n(position 22731 to 22757)', pad=-10, fontsize = 18)
    plt.tight_layout()
    plt.savefig(vait_pie_chart_path)
    plt.clf()
def find_vaits(sequence):
    query = 's+i{2,}s{2,}h{4,}s{2,}i{2,}s+'
    r = re.compile(query)
    result = [[m.start(),m.end()] for m in r.finditer(sequence)]
    return len(result)
def run_multiprocessing(function, input_list, workers):
    workers = workers
    with ProcessPoolExecutor(workers) as executor:
        results = executor.map(function, input_list)
        return list(results)
def load_tsv(file_path):
    return pd.read_csv(file_path, sep='\t')
def generate_variants_over_time(gisaid_metadata, cleaned_metadata, counts_metadata):
    # clean the big tsv first
    command = f'awk -F\'\t\' \'{{print $6"\\t"$7"\\t"$14}}\' {gisaid_metadata} > {cleaned_metadata}'
    subprocess.run(command, shell = True)
    
    #now load smaller tsv into pandas
    df = pd.read_csv(cleaned_metadata, sep = '\t')
    df.rename({'Collection date':'date', 'Location':'location', 'Pango lineage':'clade'}, axis = 1, inplace = True)
    df['location'] = df['location'].str.split(' / ').str[1]
    df = df[df['date'].astype(str).apply(len) > 4]
    df['date'] = pd.to_datetime(df['date'], errors = 'coerce')
    df['date'] = df['date'].dt.strftime('%Y-%m')
    
    count_df = df.groupby(['date', 'location', 'clade']).size().reset_index(name='sequences')
    count_df = count_df[(count_df['date'] >= '2019-12') & (count_df['date'] < '2024-08')]
    count_df = count_df[count_df['clade'] != 'Unassigned']
    count_df.to_csv(counts_metadata, sep = '\t', index = False)
def deduplicate_matching(group):
    group = group.sort_values("match_type")
    # Set all "GISAID" columns to NaN for rows not the first
    for col in group.columns:
        if "GISAID" in col:
            group.iloc[1:, group.columns.get_loc(col)] = np.nan
    return group

def generate_variants_over_time_combined(combined_metadata, gisaid_accessions,out_counts_metadata_combined,out_counts_metadata_gisaid,out_counts_metadata_ncbi, monthly_plot_outfile):

    with LocalCluster(n_workers = 20, threads_per_worker = 2, memory_limit = '100GB') as cluster:
        with Client(cluster, timeout="30s") as client:  

            df = dd.read_csv(combined_metadata, delimiter='\t', dtype=str)

            ##deduplicate multiple matches between gisaid and ncbi
            category_order = [
                'seq_and_meta_match',
                'seq_match_and_likely_meta_match',
                'seq_match_and_attempted_closest_meta_date_match',
                'seq_mismatch_and_meta_mismatch',
                'seq_mismatch_but_meta_match',
            ]
            df['match_type'] = df['match_type'].astype(pd.CategoricalDtype(categories=category_order, ordered=True))
            df_nulls = df[df['Accession ID_GISAID'].isnull()]
            df = df.dropna(subset=['Accession ID_GISAID'])
            value_counts = df['Accession ID_GISAID'].dropna().value_counts().to_frame().reset_index()
            value_counts.columns = ['Accession ID_GISAID','gisaid_acc_counts']
            
            dask_divisions = df.set_index('Accession ID_GISAID').divisions
            unique_divisions = list(dict.fromkeys(list(dask_divisions)))
            df = df.set_index('Accession ID_GISAID', divisions = unique_divisions)
            value_counts = value_counts.set_index('Accession ID_GISAID', divisions = unique_divisions)
            df = df.merge(value_counts, how = 'left', left_index = True,right_index = True)

            df = df[(df['gisaid_acc_counts'] == 1) | (df['Accession_NCBI'].notnull())]
            df_dups = df[df['gisaid_acc_counts'] > 1]
            df_dups = df_dups.reset_index()
            df = df[df['gisaid_acc_counts'] == 1]
            df_dups = df_dups.groupby('Accession ID_GISAID',group_keys = False).apply(deduplicate_matching)
            
            
            df = df.reset_index()
            
            df = dd.concat([df, df_dups, df_nulls])
            wait(df)
            df.persist()
            wait(df)
            df.to_csv('../other_data/deduplicated_combined_metadata.tsv', sep = '\t', index = False, single_file=True)
            wait(df)

            df['clade'] = df['Virus Pangolin Classification_NCBI'].where(df['Virus Pangolin Classification_NCBI'].notnull(), df['Pango lineage_GISAID'])
            df['Geographic Location_NCBI'] =df['Geographic Location_NCBI'].str.split(':').str[0]
            df['country_GISAID'] = df['Location_GISAID'].str.split(' / ').str[1]
            df['continent_GISAID'] = df['Location_GISAID'].str.split(' / ').str[0]
            df['location'] =  df['Geographic Location_NCBI'].where(df['Geographic Location_NCBI'].notnull(), df['country_GISAID'])
            df['continent'] = df['Geographic Region_NCBI'].where(df['Geographic Region_NCBI'].notnull(), df['continent_GISAID'])
            #again use this for pie chart in excel
            print(df['match_type'].value_counts().compute())
            print(df['Location_GISAID'].str.split(' / ').str[0].dropna().value_counts().compute())
            print(df['Geographic Region_NCBI'].dropna().value_counts().compute())
            print(df[df['match_type'].isin(['seq_and_meta_match','seq_match_and_likely_meta_match'])]['continent'].value_counts().compute())
            print(df['continent'].value_counts().compute())

            #only include samples with complete date data for time course plotting
            df['date'] = df['Isolate Collection date_NCBI'].where(df['Isolate Collection date_NCBI'].notnull(), df['Collection date_GISAID'])
            df = df[df['date'].astype(str).apply(len) > 4]
            df['date'] = dd.to_datetime(df['date'], errors = 'coerce')
            df['date'] = df['date'].dt.strftime('%Y-%m')
            
            #monthly plotting
            monthly_df_gisaid = df[df['Pango lineage_GISAID'].notnull()].groupby('date').size().compute().reset_index(name='GISAID Counts')
            monthly_df_gisaid = monthly_df_gisaid[(monthly_df_gisaid['date'] >= '2019-12') & (monthly_df_gisaid['date'] < '2024-08')]

            monthly_df_ncbi = df[df['Virus Pangolin Classification_NCBI'].notnull()].groupby('date').size().compute().reset_index(name='NCBI Counts')
            monthly_df_ncbi = monthly_df_ncbi[(monthly_df_ncbi['date'] >= '2019-12') & (monthly_df_ncbi['date'] < '2024-08') ]
            
            monthly_df_combined = df[df['match_type'].isin(['seq_and_meta_match','seq_match_and_likely_meta_match'])].groupby('date').size().compute().reset_index(name='Overlap Counts')
            monthly_df_combined = monthly_df_combined[(monthly_df_combined['date'] >= '2019-12') & (monthly_df_combined['date'] < '2024-08')]

            merge_monthly_df = pd.merge(monthly_df_gisaid, monthly_df_ncbi, on= 'date', how = 'outer')
            merge_monthly_df = pd.merge(merge_monthly_df, monthly_df_combined, on = 'date', how = 'outer')
            merge_monthly_df.set_index('date', inplace = True)

            merge_monthly_df.index = pd.to_datetime(merge_monthly_df.index, format='%Y-%m')
            merge_monthly_df = merge_monthly_df.sort_index()
            merge_monthly_df.index = merge_monthly_df.index.strftime('%Y-%m')

            merge_monthly_df.plot(kind = 'bar',figsize = (16, 6), width=0.9)

            plt.title('Monthly Submissions for GISAID, NCBI and the Overlap', fontsize = 16)
            plt.xlabel('Date', fontsize = 14)
            plt.ylabel('Counts', fontsize = 14)
            plt.legend(loc='upper left', fontsize = 16)
            plt.xticks(fontsize = 14)
            plt.tight_layout()
            plt.savefig(monthly_plot_outfile)
            plt.clf()


            df_combined = df[['date','location', 'clade']].compute()
            df_combined = df_combined.groupby(['date', 'location', 'clade']).size().reset_index(name='sequences')
            df_combined = df_combined[(df_combined['date'] >= '2019-12') & (df_combined['date'] < '2024-08')]
            df_combined = df_combined[~df_combined['clade'].isin(['Unassigned','unclassifiable'])]
            df_combined.to_csv(out_counts_metadata_combined, sep = '\t', index = False)

            df_gisaid = df[df['Pango lineage_GISAID'].notnull()][['date','location', 'clade']].compute()
            df_gisaid = df_gisaid.groupby(['date', 'location', 'clade']).size().reset_index(name='sequences')
            df_gisaid = df_gisaid[(df_gisaid['date'] >= '2019-12') & (df_gisaid['date'] < '2024-08')]
            df_gisaid = df_gisaid[~df_gisaid['clade'].isin(['Unassigned','unclassifiable'])]
            df_gisaid.to_csv(out_counts_metadata_gisaid, sep = '\t', index = False)

            df_ncbi = df[df['Virus Pangolin Classification_NCBI'].notnull()][['date','location', 'clade']].compute()
            df_ncbi = df_ncbi.groupby(['date', 'location', 'clade']).size().reset_index(name='sequences')
            df_ncbi = df_ncbi[(df_ncbi['date'] >= '2019-12') & (df_ncbi['date'] < '2024-08')]
            df_ncbi = df_ncbi[~df_ncbi['clade'].isin(['Unassigned','unclassifiable'])]
            df_ncbi.to_csv(out_counts_metadata_ncbi, sep = '\t', index = False)

def calc_type_muts(mutations_string, syn_or_nonsyn):
    mutations_string = mutations_string.strip("[]").replace("'","")
    mutations_list = mutations_string.split(',')
    mut_type_count= [i for i in mutations_list if i.find('_' + syn_or_nonsyn) > -1]
    return len(mut_type_count)
def gather_inputs(result_dir, seq_type):

    codon_data = os.path.join(result_dir, f'codon_usage_{seq_type}_seqs.tsv')    
    dn_data = os.path.join(result_dir, f'dinucleotide_usage_{seq_type}_seqs.tsv')
    jdn_data = os.path.join(result_dir, f'jn_dinucleotide_usage_{seq_type}_seqs.tsv')
    enc_cai_data = os.path.join(result_dir, f'enc_cai_{seq_type}_seqs.tsv')
    mut_data = os.path.join(result_dir, f'mutation_data_{seq_type}_seqs.tsv')
    cds_data = os.path.join(result_dir, f'split_genes_{seq_type}_seqs.tsv')
    ss_data = os.path.join(result_dir, f'all_lineages_added_aligned_{seq_type}_seqs_final_secondary_structure_file.tsv')
    syn_removed_ss_data = os.path.join(result_dir, f'syn_removed_{seq_type}_seqs_final_secondary_structure_file.tsv')
    return codon_data, dn_data, jdn_data, enc_cai_data, mut_data, cds_data, ss_data, syn_removed_ss_data

def gather_outputs(outdir, seq_type):
    #output filenames
    #image filenames
    vait_pie_chart_path = os.path.join(outdir, f'vait_pie_chart_{seq_type}.png')
    vait_plot = os.path.join(outdir, f'vait_plot_{seq_type}.png')
    autocorrelate_out = os.path.join(outdir, f'autocorrelate_{seq_type}.png')
    cai_heatmap_out = os.path.join(outdir, f'CAI_heatmap_{seq_type}.png')
    cai_line_graph_out = os.path.join(outdir, f'cai_plot_{seq_type}.png')
    cpai_heatmap_out = os.path.join(outdir, f'CPAI_heatmap_{seq_type}.png')
    dinucleotide_heatmap_out = os.path.join(outdir, f'dinucleotide_heatmap_{seq_type}.png')
    enc_heatmap_out = os.path.join(outdir, f'ENC_heatmap_{seq_type}.png')
    enc_line_graph_out = os.path.join(outdir, f'enc_plot_{seq_type}.png')
    gc_heatmap_out = os.path.join(outdir, f'gc_abundance_heatmap_{seq_type}.png')
    junc_heatmap_out = os.path.join(outdir, f'junction_dinucleotide_heatmap_{seq_type}.png')
    codon_heatmap_out = os.path.join(outdir, f'codon_abundance_heatmap_{seq_type}.png')
    codon_line_graph_out = os.path.join(outdir, f'codon_line_plot_{seq_type}.png')
    free_energy_plot = os.path.join(outdir, f'free_energies_{seq_type}.png')
    mutations_plot = os.path.join(outdir, f'mutations_{seq_type}.png')
    secondary_structures_plot = os.path.join(outdir, f'structures_{seq_type}.png')
    mutation_counts_file = os.path.join(outdir, f'mutation_counts_{seq_type}.csv')
    structure_table_file = os.path.join(outdir, f'structure_table_{seq_type}.csv')
    Table_S1_codon_usage = os.path.join(outdir, f'Table_S1_codon_usage_{seq_type}.csv')
    Table_S2_dn_usage = os.path.join(outdir, f'Table_S2_dn_usage_{seq_type}.csv')
    Table_S3_jdn_usage = os.path.join(outdir, f'Table_S3_jdn_usage_{seq_type}.csv')
    Table_S4_gc_usage = os.path.join(outdir, f'Table_S4_gc_usage_{seq_type}.csv')
    Table_S5_cai = os.path.join(outdir, f'Table_S5_cai_{seq_type}.csv')
    Table_S6_enc = os.path.join(outdir, f'Table_S6_enc_{seq_type}.csv')
    Table_S7_cpai = os.path.join(outdir, f'Table_S7_cpai_{seq_type}.csv')
    return vait_pie_chart_path,vait_plot,autocorrelate_out,cai_heatmap_out,cai_line_graph_out,cpai_heatmap_out,dinucleotide_heatmap_out,enc_heatmap_out,enc_line_graph_out,gc_heatmap_out,junc_heatmap_out,codon_heatmap_out,codon_line_graph_out,free_energy_plot,mutations_plot,secondary_structures_plot,mutation_counts_file,structure_table_file,Table_S1_codon_usage,Table_S2_dn_usage,Table_S3_jdn_usage,Table_S4_gc_usage,Table_S5_cai,Table_S6_enc,Table_S7_cpai
def create_variant_over_time_df(variant_data):
    variant_df = pd.read_csv(variant_data, sep = '\t')
    variant_df['date'] = pd.to_datetime(variant_df['date'])
    variant_df.sort_values(by='date', inplace = True)
    variant_df['date'] = variant_df['date'].dt.to_period('M')
    variant_df['date'] = variant_df['date'].dt.to_timestamp()
    variant_df = variant_df.groupby(['date', 'clade'])['sequences'].sum().reset_index()

    variant_df = variant_df.sort_values(by=['date', 'sequences'], ascending=[True, False])
    # Resetting the index to make 'date' a column again
    variant_df.reset_index(inplace=True)
    variant_df = variant_df.T
    variant_df.columns = variant_df.loc['date']
    variant_df.drop(index = ['date'], inplace= True)
    variant_df.index.name = None
    variant_df.columns.name = None
    return variant_df
def get_top_variants(group, threshold=0.25):
    total = group['sequences'].sum()
    group = group.sort_values(by='sequences', ascending=False)
    group['cumulative_percent'] = group['sequences'].cumsum() / total
    group['percent'] = group['sequences'] / total *100
    group = group.head(2)
    group['clade_and_percent'] = group['clade'] + ' (' + group['percent'].round(1).astype(str) + ')'
    print(list(group['clade_and_percent']))
    return ' | '.join(list(group['clade_and_percent']))
    

if __name__ == '__main__':
    #this script generates all the figures and supplemental data in our accompanying manuscript. 
    #time series data is based only on the combined dataset. 
    #input filenames
    combined_metadata = '../comparison_data/check_ncbi_in_gisaid_final.tsv'
    # #after deduplication fix
    global_population_filename = '../other_data/global_population_UN.csv'
    gisaid_accessions = '../other_data/accession_ids_11-19-2024.csv'

    #variant proportion files, derived from combined metadata
    counts_metadata_combined = '../other_data/lineage_counts_combined.tsv'
    counts_metadata_gisaid = '../other_data/lineage_counts_gisaid.tsv'
    counts_metadata_ncbi = '../other_data/lineage_counts_ncbi.tsv'

    #outfiles based solely on metadata
    monthly_plot_outfile = '../other_data/GISAID_vs_NCBI_counts_per_month.png'
    country_submission_outfile = '../other_data/country_submission_plot.png'
    variant_monthly_summary = '../other_data/variant_monthly_summary_top_25_percent.csv'


    input_dirs = [
        '../comparison_data/final_files',
        '../gisaidData/final_files',
        '../ncbiData/final_files'
        
        ]
    seq_types = [
        'ancestral',
        'representative',
        'consensus'
    ]

    # #clean metadata
    #somtimes this results in a heartbeat scheduling error from Dask - this can be ignored as the output file
    #is still generated and is complete. You can comment out this line for the subsequent iteration and just use 
    #the output file (counts metadata). 
    #this function prints to screen the Data required for all of our pie charts in the manuscript
    # The actual plotting is done in excel for these
    generate_variants_over_time_combined(combined_metadata,
                                         gisaid_accessions,
                                         counts_metadata_combined,
                                         counts_metadata_gisaid,
                                         counts_metadata_ncbi,
                                         monthly_plot_outfile)

    # # Get Country Distribution
    variant_counts_by_country_df_combined = pd.read_csv(counts_metadata_combined, sep = '\t')
    variant_only_per_month = variant_counts_by_country_df_combined.groupby(['date', 'clade'])['sequences'].sum().reset_index().sort_values(by=['date','sequences'], ascending = False)
    variant_only_per_month.to_csv('../other_data/variant_count_per_month.csv', index = False)
    country_submissions_per_month = variant_counts_by_country_df_combined.groupby(['date', 'location'])['sequences'].sum().reset_index().sort_values(by=['date','sequences'], ascending = [True, False])
    country_submissions_per_month.to_csv('../other_data/country_submissions_per_month.csv', index = False)
    country_totals = variant_counts_by_country_df_combined.groupby(['location'])['sequences'].sum().reset_index().sort_values(by=['sequences'], ascending = [False])
    total_seqs = country_totals['sequences'].sum()
    country_totals['cumulative_percent'] = country_totals['sequences'].cumsum() / total_seqs
    countries_to_achieve_80_percent = pd.concat([country_totals[country_totals['cumulative_percent'] < 0.8], country_totals[country_totals['cumulative_percent'] > 0.8].head(1)])
    country_totals['location'] = country_totals['location'].apply(lambda x: x if x in list(countries_to_achieve_80_percent['location']) else 'All Others')
    country_totals = country_totals.groupby(['location'])['sequences'].sum()
    country_totals = (country_totals / country_totals.sum()) * 100

    #get global population
    global_pop_df = pd.read_csv(global_population_filename)
    global_pop_df['Population_millions'] = global_pop_df['Population_millions'].str.replace(',','').astype(float)
    global_pop_df['Country'] = global_pop_df['Country'].where(global_pop_df['Country'].isin(list(country_totals.index)), 'All Others')
    global_pop_df = global_pop_df.groupby('Country')['Population_millions'].sum()
    global_pop_df = (global_pop_df / global_pop_df.sum()) * 100
    #copy the below into excel for plotting
    print(country_totals)
    print(global_pop_df)

    #plot submissions over time by top countries
    country_submissions_per_month['location'] =  country_submissions_per_month['location'].apply(lambda x: x if x in list(countries_to_achieve_80_percent['location']) else 'All Others')
    country_submissions_per_month = country_submissions_per_month.groupby(['date', 'location'])['sequences'].sum().reset_index().sort_values(by=['date','sequences'], ascending = [True, True])
    country_submissions_per_month = country_submissions_per_month.pivot(index = 'date', columns = 'location', values = 'sequences').fillna(0)
    country_submissions_per_month = country_submissions_per_month.divide(country_submissions_per_month.sum(axis = 1), axis = 0) * 100
    country_submissions_per_month.columns = ['All Others'] + [i for i in country_submissions_per_month.columns if i != 'All Others']
    cmap = ListedColormap(sns.color_palette("bright"))
    default_colors = plt.cm.get_cmap(cmap, len(country_submissions_per_month.columns))
    custom_colors = {'All Others':'black',
                     'USA': 'olive'}
    colors = [custom_colors.get(col, default_colors(i)) for i, col in enumerate(country_submissions_per_month.columns)]
    country_submissions_per_month.plot(kind = 'bar',stacked = True, figsize = (16, 6), color = colors)
    plt.title('Percentage of Sequences Submitted per Country Over Time', fontsize = 16)
    plt.xlabel('Date', fontsize = 14)
    plt.ylabel('% Submissions', fontsize = 14)
    plt.xticks(fontsize = 14)
    plt.legend( bbox_to_anchor=(1.05, 1), loc='upper left', fontsize = 16)
    plt.tight_layout()
    plt.savefig(country_submission_outfile, dpi = 600)
    plt.clf()

    variant_counts_by_country_df_gisaid = pd.read_csv(counts_metadata_gisaid, sep = '\t')
    variant_counts_by_country_df_ncbi = pd.read_csv(counts_metadata_ncbi, sep = '\t')

    variant_counts_by_country_df_combined = variant_counts_by_country_df_combined.groupby(['date', 'clade'])['sequences'].sum().reset_index(name='sequences').sort_values(by=['date','sequences'], ascending = [True, False])
    top_25_percent_variants_combined = variant_counts_by_country_df_combined.groupby('date').apply(get_top_variants)

    variant_counts_by_country_df_gisaid = variant_counts_by_country_df_gisaid.groupby(['date', 'clade'])['sequences'].sum().reset_index(name='sequences').sort_values(by=['date','sequences'], ascending = [True, False])
    top_25_percent_variants_gisaid = variant_counts_by_country_df_gisaid.groupby('date').apply(get_top_variants)

    variant_counts_by_country_df_ncbi = variant_counts_by_country_df_ncbi.groupby(['date', 'clade'])['sequences'].sum().reset_index(name='sequences').sort_values(by=['date','sequences'], ascending = [True, False])
    top_variants_overall = list(variant_counts_by_country_df_combined.groupby('clade')['sequences'].sum().sort_values(ascending=False).head(20).index)
    top_25_percent_variants_ncbi = variant_counts_by_country_df_ncbi.groupby('date').apply(get_top_variants)
    variant_summary_per_month = pd.DataFrame(pd.concat([top_25_percent_variants_gisaid, top_25_percent_variants_ncbi, top_25_percent_variants_combined], axis =1))
    variant_summary_per_month.columns = ['GISAID Counts', 'NCBI Counts','Combined Counts']
    variant_summary_per_month.to_csv('../other_data/variant_monthly_summary.csv')
    #calculate standard deviation for codon usages and secondary structure
    for seq_type in seq_types:
        print(seq_type)
        seq_type_codon_dfs = []
        seq_type_ss_dfs = []
        syn_removed_ss_dfs = []
        for input_dir in input_dirs:
            codon_data, dn_data, jdn_data, enc_cai_data, mut_data, cds_data, ss_data, syn_removed_ss_data = gather_inputs(input_dir, seq_type)
            outdir = os.path.join(os.path.dirname(input_dir), 'plots_and_tables')

            #codon loading
            codon_gc_df = pd.read_csv(codon_data, sep = '\t')
            codon_gc_df = codon_gc_df[codon_gc_df['gene'] == 'whole_cds'].drop(['parent_lineages', 'gene', 'GC1%', 'GC2%', 'GC3%', 'GC%', 'adjusted_cds_length'], axis = 1).set_index('terminal_lineage')
            codon_gc_df = codon_gc_df.subtract(codon_gc_df.loc['referenceB']) 
            seq_type_codon_dfs.append(codon_gc_df)

            #ss loading
            ss_df = pd.read_csv(ss_data, sep = '\t')
            ss_df = ss_df[['terminal_lineage','ensemble_free_energy']].set_index('terminal_lineage')
            ss_df = ss_df.subtract(ss_df.loc['referenceB'])
            seq_type_ss_dfs.append(ss_df)


        combined_codon_df = pd.concat(seq_type_codon_dfs, keys = [os.path.dirname(i) for i in input_dirs])
        stdev_df_codon = combined_codon_df.groupby(level =1).std()
        print(stdev_df_codon.mean().mean())
        combined_ss_df = pd.concat(seq_type_ss_dfs, axis = 1).std(axis = 1).mean()
        print(combined_ss_df)
        combined_syn_removed_ss_df = pd.concat(syn_removed_ss_dfs, axis = 1).std(axis = 1).mean()
        print(combined_syn_removed_ss_df)

    
    
    for input_dir in input_dirs:
        for seq_type in seq_types:
            print(f'---------NOW RUNNING: {input_dir}, {seq_type}----------')
            #inputs
            codon_data, dn_data, jdn_data, enc_cai_data, mut_data, cds_data, ss_data, syn_removed_ss_data = gather_inputs(input_dir, seq_type)
            #outputs
            outdir = os.path.join(os.path.dirname(input_dir), 'plots_and_tables')
            
            vait_pie_chart_path,vait_plot,autocorrelate_out,cai_heatmap_out,cai_line_graph_out,cpai_heatmap_out,dinucleotide_heatmap_out,enc_heatmap_out,enc_line_graph_out,gc_heatmap_out,junc_heatmap_out,codon_heatmap_out,codon_line_graph_out,free_energy_plot,mutations_plot,secondary_structures_plot,mutation_counts_file,structure_table_file,Table_S1_codon_usage,Table_S2_dn_usage,Table_S3_jdn_usage,Table_S4_gc_usage,Table_S5_cai,Table_S6_enc,Table_S7_cpai = gather_outputs(outdir, seq_type)
            
            
            filenames_list = [codon_data,
                            dn_data,
                            jdn_data,
                            enc_cai_data,
                            mut_data,
                            cds_data,
                            ss_data,
                            syn_removed_ss_data]
            #load data
            codon_gc_df, dn_df,jdn_df,enc_cai_df,mutation_df,cds_df, ss_df, syn_removed_ss_df = run_multiprocessing(load_tsv, filenames_list, 9)
            
            variant_proportion_df = create_variant_over_time_df(counts_metadata_combined)
            vait_pie_input_df = ss_df.copy()
            
            #result dicts and lists to hold variant proportion weighted data per date
            dates = sorted(list(variant_proportion_df.columns))

            wt_fes = []
            syn_removed_fes = []
            syn_mut_counts = []
            nonsyn_mut_counts = []
            gc_counts_dict = {'GC1%': [], 'GC2%': [], 'GC3%': [], 'GC%': []}

            codon_counts_dict = {'TTT': [], 'TTC': [], 'TTA': [], 'TTG': [], 'CTT': [], 'CTC': [], 'CTA': [], 'CTG': [], 'ATT': [], 'ATC': [], 'ATA': [], 'ATG': [], 'GTT': [], 'GTC': [], 'GTA': [], 'GTG': [], 'TAT': [], 'TAC': [], 'TAA': [], 'TAG': [], 'CAT': [], 'CAC': [], 'CAA': [], 'CAG': [], 'AAT': [], 'AAC': [], 'AAA': [], 'AAG': [], 'GAT': [], 'GAC': [], 'GAA': [], 'GAG': [], 'TCT': [], 'TCC': [], 'TCA': [], 'TCG': [], 'CCT': [], 'CCC': [], 'CCA': [], 'CCG': [], 'ACT': [], 'ACC': [], 'ACA': [], 'ACG': [], 'GCT': [], 'GCC': [], 'GCA': [], 'GCG': [], 'TGT': [], 'TGC': [], 'TGA': [], 'TGG': [], 'CGT': [], 'CGC': [], 'CGA': [], 'CGG': [], 'AGT': [], 'AGC': [], 'AGA': [], 'AGG': [], 'GGT': [], 'GGC': [], 'GGA': [], 'GGG': []}
            
            dinucleotide_counts_dict = {'TT': [], 'TC': [], 'TA': [], 'TG': [], 'CT': [], 'CC': [], 'CA': [], 'CG': [], 'AT': [], 'AC': [], 'AA': [], 'AG': [], 'GT': [], 'GC': [], 'GA': [], 'GG': []}
            
            jdinucleotide_counts_dict = {'TT': [], 'TC': [], 'TA': [], 'TG': [], 'CT': [], 'CC': [], 'CA': [], 'CG': [], 'AT': [], 'AC': [], 'AA': [], 'AG': [], 'GT': [], 'GC': [], 'GA': [], 'GG': []}
            
            enc_counts_dict = {'ENC': [], 'ENC_GC_corrected': [], 'ENC_CP': [], 'ENC_CP_GC_corrected': []}
            
            struct_counts_dict = {'f_count': [], 's_count': [], 'i_count': [], 'm_count': [], 'h_count': [], 't_count': []}
                
            cai_dict = {'whole cds CAI': [], 'e cds CAI': [], 'm cds CAI': [], 'n cds CAI': [], 'orf10 cds CAI': [], 'orf1ab cds CAI': [], 'orf3a cds CAI': [], 'orf3b cds CAI': [], 'orf6 cds CAI': [], 'orf7a cds CAI': [], 'orf7b cds CAI': [], 'orf8 cds CAI': [], 'spike cds CAI': []}
            
            cpai_dict = {'whole cds CPAI': [], 'e cds CPAI': [], 'm cds CPAI': [], 'n cds CPAI': [], 'orf10 cds CPAI': [], 'orf1ab cds CPAI': [], 'orf3a cds CPAI': [], 'orf3b cds CPAI': [], 'orf6 cds CPAI': [], 'orf7a cds CPAI': [], 'orf7b cds CPAI': [], 'orf8 cds CPAI': [], 'spike cds CPAI': []}
            
            vait_dict = {'stem-loop hits': []}
            
            #trim/reformat dataframes
            #codon/dinucleotide reformatting
            codon_gc_df = codon_gc_df[codon_gc_df['gene'] == 'whole_cds'].set_index('terminal_lineage')[list(gc_counts_dict.keys())+list(codon_counts_dict.keys())]
            codon_gc_df[list(codon_counts_dict.keys())] = codon_gc_df[list(codon_counts_dict.keys())].divide(codon_gc_df[list(codon_counts_dict.keys())].sum(axis = 1), axis = 0) * 100
            codon_gc_df[list(codon_counts_dict.keys())] = codon_gc_df[list(codon_counts_dict.keys())].fillna(0)

            dn_df = dn_df[dn_df['gene'] == 'whole_cds'].set_index('terminal_lineage')[list(dinucleotide_counts_dict.keys())]
            dn_df[list(dinucleotide_counts_dict.keys())] = dn_df[list(dinucleotide_counts_dict.keys())].divide(dn_df[list(dinucleotide_counts_dict.keys())].sum(axis = 1), axis = 0) * 100
            dn_df = dn_df.fillna(0)
            jdn_df = jdn_df[jdn_df['gene'] == 'whole_cds'].set_index('terminal_lineage')[list(jdinucleotide_counts_dict.keys())]
            jdn_df[list(jdinucleotide_counts_dict.keys())] = jdn_df[list(jdinucleotide_counts_dict.keys())].divide(jdn_df[list(jdinucleotide_counts_dict.keys())].sum(axis = 1), axis = 0) * 100
            jdn_df.columns = [i + '_jdn' for i in jdn_df.columns]
            jdn_df = jdn_df.fillna(0)
            #stats
            cai_df = enc_cai_df.pivot_table(index = 'terminal_lineage', columns = 'gene', values = ['CAI','CPAI'])
            cai_df.columns =  [' cds '.join(col[::-1]) for col in cai_df.columns]
            cai_df.rename({'whole_cds cds CPAI': 'whole cds CPAI', 'whole_cds cds CAI': 'whole cds CAI'},inplace = True, axis = 1)
            enc_df = enc_cai_df[enc_cai_df['gene'] == 'whole_cds'].set_index('terminal_lineage')[enc_counts_dict.keys()]
            #mutation reformatting
            mutation_df = mutation_df.set_index('terminal_lineage')[[i for i in mutation_df.columns if '_cds' in i]]
            syn_df = mutation_df.copy()
            for column in syn_df.columns:
                syn_df[column] = syn_df[column].apply(calc_type_muts, args=['SYNONYMOUS'])
            syn_df['total_syn'] = syn_df.sum(axis = 1)
            syn_df = syn_df[['total_syn']]
            syn_df = syn_df.fillna(0)

            nonsyn_df = mutation_df.copy()
            for column in nonsyn_df.columns:
                nonsyn_df[column] = nonsyn_df[column].apply(calc_type_muts, args=['NONSYNONYMOUS'])
            nonsyn_df['total_nonsyn'] = nonsyn_df.sum(axis = 1)
            nonsyn_df = nonsyn_df[['total_nonsyn']]
            nonsyn_df = nonsyn_df.fillna(0)

            #secondary structure reformatting
            ss_df = ss_df.set_index('terminal_lineage')
            ss_df['stem-loop hits'] = ss_df['forgi_output'].apply(find_vaits)
            ss_df.rename({'f':'f_count', 's':'s_count', 'i':'i_count', 'm':'m_count', 'h':'h_count', 't':'t_count'},inplace = True, axis = 1)
            ss_df = ss_df[['stem-loop hits'] + list(struct_counts_dict.keys()) + ['ensemble_free_energy']]
            ss_df[['stem-loop hits'] + list(struct_counts_dict.keys())] = ss_df[['stem-loop hits'] + list(struct_counts_dict.keys())].fillna(0)
            syn_removed_ss_df = syn_removed_ss_df.set_index('terminal_lineage')[['ensemble_free_energy']].rename({'ensemble_free_energy':'syn_removed_ensemble_free_energy'}, axis = 1)
            
            #subtract wt
            cai_df = cai_df.subtract(cai_df.loc['referenceB']) 
            codon_gc_df = codon_gc_df.subtract(codon_gc_df.loc['referenceB']) 
            codon_gc_df.to_csv('codon_subtract.csv')
            enc_df = enc_df.subtract(enc_df.loc['referenceB']) 
            ss_df = ss_df.subtract(ss_df.loc['referenceB'])
            syn_removed_ss_df = syn_removed_ss_df.subtract(syn_removed_ss_df.loc['referenceB'])
            dn_df = dn_df.subtract(dn_df.loc['referenceB']) 
            jdn_df = jdn_df.subtract(jdn_df.loc['referenceB'])
            
            #combine all the dataframes
            merge_df = pd.merge(variant_proportion_df, codon_gc_df, left_index=True, right_index=True, how = 'left')
            merge_df = pd.merge(merge_df, dn_df,left_index=True, right_index=True, how = 'left')
            merge_df = pd.merge(merge_df, jdn_df,left_index=True, right_index=True, how = 'left')
            merge_df = pd.merge(merge_df, cai_df,left_index=True, right_index=True, how = 'left')
            merge_df = pd.merge(merge_df, enc_df,left_index=True, right_index=True, how = 'left')
            merge_df = pd.merge(merge_df, syn_df,left_index=True, right_index=True, how = 'left')
            merge_df = pd.merge(merge_df, nonsyn_df,left_index=True, right_index=True, how = 'left')
            merge_df = pd.merge(merge_df, ss_df,left_index=True, right_index=True, how = 'left')
            merge_df = pd.merge(merge_df, syn_removed_ss_df,left_index=True, right_index=True, how = 'left')
            merge_df.dropna(inplace = True)

            for i in dates:
                try:

                    variant_proportions = merge_df[i] / merge_df[i].sum()
                except:
                    print(input_dir)
                    print(seq_type)
                    print(i)
                    sys.exit()

                wt_fe = variant_proportions * merge_df['ensemble_free_energy'].dropna()
                wt_fe = wt_fe.sum()

                syn_removed_fe = variant_proportions * merge_df['syn_removed_ensemble_free_energy']
                syn_removed_fe = syn_removed_fe.sum()

                nonsyn_muts = variant_proportions * merge_df['total_nonsyn']
                nonsyn_muts = nonsyn_muts.sum()

                syn_muts = variant_proportions * merge_df['total_syn']
                syn_muts = syn_muts.sum()

                for j in gc_counts_dict.keys():
                    gc_count = variant_proportions * merge_df[j].dropna()
                    gc_count = gc_count.sum()
                    gc_counts_dict[j].append(gc_count)

                for j in codon_counts_dict.keys():
                    codon_counts = variant_proportions * merge_df[j].dropna()
                    codon_counts = float(codon_counts.sum())
                    codon_counts_dict[j].append(codon_counts)  

                for j in enc_counts_dict.keys():
                    enc_counts = variant_proportions * merge_df[j].dropna()
                    enc_counts = float(enc_counts.sum())
                    enc_counts_dict[j].append(enc_counts) 

                for j in struct_counts_dict.keys():
                    struct_counts = variant_proportions * merge_df[j].dropna()
                    struct_counts = float(struct_counts.sum())
                    struct_counts_dict[j].append(struct_counts)

                for j in dinucleotide_counts_dict.keys():
                    dn_counts = variant_proportions * merge_df[j].dropna()
                    dn_counts = float(dn_counts.sum())
                    dinucleotide_counts_dict[j].append(dn_counts)

                for j in jdinucleotide_counts_dict.keys():
                    j += '_jdn'
                    jdn_counts = variant_proportions * merge_df[j].dropna()
                    jdn_counts = float(jdn_counts.sum())
                    jdinucleotide_counts_dict[j.replace('_jdn','')].append(jdn_counts)

                for j in cai_dict.keys():
                    cai_counts = variant_proportions * merge_df[j].dropna()
                    cai_counts = float(cai_counts.sum())
                    cai_dict[j].append(cai_counts)

                for j in cpai_dict.keys():
                    cpai_counts = variant_proportions * merge_df[j].dropna()
                    cpai_counts = float(cpai_counts.sum())
                    cpai_dict[j].append(cpai_counts) 

                for j in vait_dict.keys():
                    vait_counts = variant_proportions * merge_df[j].dropna()
                    vait_counts = float(vait_counts.sum())
                    vait_dict[j].append(vait_counts)  


                #uncomment below line after ss calcs
                wt_fes.append(wt_fe)
                syn_removed_fes.append(syn_removed_fe)
                syn_mut_counts.append(syn_muts)
                nonsyn_mut_counts.append(nonsyn_muts)
            
            # #now start plotting
            string_dates = ['-'.join(str(i).split(' ')[0].split('-')[0:2]) for i in dates]
            print(string_dates)
            print(wt_fes)

            # Calculate trendlines
            wt_fes_coeffs = np.polyfit(range(len(dates)), wt_fes, 1)  # Linear fit for WT data
            wt_fes_trendline = np.polyval(wt_fes_coeffs, range(len(dates)))  # Trendline values

            syn_removed_fes_coeffs = np.polyfit(range(len(dates)), syn_removed_fes, 1)  # Linear fit for Synonymous Removed data
            syn_removed_fes_trendline = np.polyval(syn_removed_fes_coeffs, range(len(dates)))  # Trendline values

            # Plot data and trendlines
            plt.plot(dates, wt_fes, label='Ensemble Free Energy - Difference From WT', marker='o')
            plt.plot(dates, wt_fes_trendline, label=f'Trendline WT (m={wt_fes_coeffs[0]:.2f})', linestyle='--')

            plt.plot(dates, syn_removed_fes, label='Synonymous Removed Ensemble Free Energy', marker='o')
            plt.plot(dates, syn_removed_fes_trendline, label=f'Trendline Syn Removed (m={syn_removed_fes_coeffs[0]:.2f})', linestyle='--')

            # Add labels, legend, and other details
            plt.legend()
            plt.ylabel('Ensemble Free Energy difference from WT (kcal/mol)', fontsize = 16)
            plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval = 4))
            plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
            plt.gca().set_xticklabels(plt.gca().xaxis.get_majorticklabels(), rotation = 45, ha='right', rotation_mode = 'anchor', fontsize = 14)
            plt.gca().set_yticklabels(plt.gca().yaxis.get_majorticklabels(), fontsize = 14)
            # plt.xticks(rotation = 90)
            plt.title('Ensemble Free Energy over Time (Difference from WT)', fontsize = 20)
            plt.tight_layout()


            # Save the plot
            plt.savefig(free_energy_plot)
            plt.clf()

            plot_acf(wt_fes, lags=3)
            plt.ylabel('Autocorrelation Coefficient', fontsize = 12)
            plt.xlabel('# Lags', fontsize = 12)
            plt.title('Autocorrelation for Ensemble Free Energy over Time', fontsize = 12)
            ax = plt.gca()
            ax.set_ylim([-1.1, 1.1])
            ax.tick_params(axis = 'both', labelsize = 12)

            plt.tight_layout()
            plt.savefig(autocorrelate_out)
            plt.clf()
            # adftest = adfuller(wt_fes, autolag=None, maxlag=5, regression = 'ct')
            adftest = adfuller(wt_fes, autolag='AIC',regression = 'ct')
            print("ADF Test Results")
            print("Null Hypothesis: The series has a unit root (non-stationary)")
            print("ADF-Statistic:", adftest[0])
            print("P-Value:", adftest[1])
            print("Number of lags:", adftest[2])
            print("Number of observations:", adftest[3])
            print("Critical Values:", adftest[4])
            print("Note: If P-Value is smaller than 0.05, we reject the null hypothesis and the series is stationary")
                
            #make mutation counts plot
            df = pd.DataFrame({'nonsyn muts':nonsyn_mut_counts, 'syn muts': syn_mut_counts}, index = string_dates)
            df.to_csv(mutation_counts_file)
            # Calculate trendlines
            nonsyn_mut_counts_coeffs = np.polyfit(range(len(dates)), nonsyn_mut_counts, 1)  # Linear fit for Nonsynonymous Mutations
            nonsyn_mut_counts_trendline = np.polyval(nonsyn_mut_counts_coeffs, range(len(dates)))  # Trendline values

            syn_mut_counts_coeffs = np.polyfit(range(len(dates)), syn_mut_counts, 1)  # Linear fit for Synonymous Mutations
            syn_mut_counts_trendline = np.polyval(syn_mut_counts_coeffs, range(len(dates)))  # Trendline values

            # Plot data and trendlines
            plt.plot(dates, nonsyn_mut_counts, label='Nonsynonymous Mutations', color='green', marker='o')
            plt.plot(dates, nonsyn_mut_counts_trendline, label=f'Trendline Nonsyn (m={nonsyn_mut_counts_coeffs[0]:.2f})', linestyle='--', color='green')

            plt.plot(dates, syn_mut_counts, label='Synonymous Mutations', color='purple', marker='o')
            plt.plot(dates, syn_mut_counts_trendline, label=f'Trendline Syn (m={syn_mut_counts_coeffs[0]:.2f})', linestyle='--', color='purple')

            # Add labels, legend, and other details
            plt.legend()
            plt.ylabel('Mutation Count', fontsize = 16)

            # Configure x-axis ticks
            plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=4))
            plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
            plt.gca().set_xticklabels(plt.gca().xaxis.get_majorticklabels(), rotation=45, ha='right', rotation_mode='anchor', fontsize = 14)
            plt.gca().set_yticklabels(plt.gca().yaxis.get_majorticklabels(), fontsize = 14)
            plt.title('Mutations over Time', fontsize = 20)
            plt.tight_layout()

            # Save the plot
            plt.savefig(mutations_plot)
            plt.clf()
            # Make a vait plot with trendlines
            for i in vait_dict.keys():
                # Calculate the trendline for the current series
                trendline_coeffs = np.polyfit(range(len(dates)), vait_dict[i], 1)  # Linear fit
                trendline_values = np.polyval(trendline_coeffs, range(len(dates)))  # Trendline values

                # Plot the data and trendline
                plt.plot(dates, vait_dict[i], label=f'{i}', marker='o')
                plt.plot(dates, trendline_values, label=f'Trendline {i} (m={trendline_coeffs[0]:.2f})', linestyle='--')

            # Add labels, legend, and other details
            plt.legend()
            plt.title('Stem Loop Hits over Time (compared to WT)')
            plt.xticks(rotation=45, ha='right', rotation_mode='anchor')
            plt.ylabel('# Stem Loops added/lost compared to WT')
            plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=4))
            plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
            plt.tight_layout()

            # Save the vait plot
            plt.savefig(vait_plot)
            plt.clf()

            #make ENC plot
            for i in enc_counts_dict.keys():
                plt.plot(dates, enc_counts_dict[i], label = i)
            plt.legend()
            plt.xticks(rotation = 45)
            plt.tight_layout()
            plt.savefig(enc_line_graph_out)
            plt.clf()

            #make CAI plot
            for i in cai_dict.keys():
                plt.plot(dates, cai_dict[i], label = i)
            plt.legend()
            plt.xticks(rotation = 45)
            plt.tight_layout()
            plt.savefig(cai_line_graph_out)
            plt.clf()

            #make codons plot
            for i in codon_counts_dict.keys():
                # print(dates)
                # print(codon_counts_dict[i])
                plt.plot(dates, codon_counts_dict[i], label = i)
                # break
            plt.legend()
            plt.xticks(rotation = 45)
            plt.tight_layout()
            plt.savefig(codon_line_graph_out)
            plt.clf()

            #make structure plot

            replace_dict = {'f_count': "5' unpaired nucleotides", 't_count':"3' unpaired nucleotides", 's_count':"Stem nucleotides", 'i_count':"Interior Loop nucleotides", 'm_count':"Multiloop nucleotides", "h_count":"Hairpin Loop nucleotides"}
            colors = sns.color_palette('tab10', len(struct_counts_dict))
            plt.figure(figsize=(10, 5))
            df = pd.DataFrame.from_dict(struct_counts_dict, orient = 'index', columns = string_dates)
            df.to_csv(structure_table_file)
            # Iterate over the keys in the dictionary
            fig, axs = plt.subplots(2, 3, figsize=(15, 10), sharex=True, sharey=True)  # Adjust rows/columns as needed
            for ax, (i, color) in zip(axs.flat, zip(struct_counts_dict.keys(), colors)):
                numerical_dates = range(len(dates))  # Or use numerical conversion if needed
                coeffs = np.polyfit(numerical_dates, struct_counts_dict[i], 1)
                trendline = np.polyval(coeffs, numerical_dates)
                
                # Plot the data and trendline
                ax.plot(dates, struct_counts_dict[i], label=f'{replace_dict[i]}', marker='o', color=color)
                ax.plot(dates, trendline, label=f'Trendline (m={coeffs[0]:.2f})', linestyle='--', color=color)
                
                # Set title and legend for each subplot
                ax.set_title(replace_dict[i])
                ax.legend()

            fig.tight_layout()
            plt.savefig(secondary_structures_plot)
            plt.clf()


            #CODON HEATMAP

            df = pd.DataFrame.from_dict(codon_counts_dict, orient='index',
                                columns=string_dates)
            
            df.index = [cb.utils.translate(i, return_str= True) + '-'+i for i in df.index]
            df.sort_index(inplace=True)
            amino_acids = [i.split('-')[0] for i in df.index]
            df.index = [i.split('-')[-1] for i in df.index]

            amino_acid_dict = {
                'A': 'Ala',  # Alanine
                'R': 'Arg',  # Arginine
                'N': 'Asn',  # Asparagine
                'D': 'Asp',  # Aspartic Acid
                'C': 'Cys',  # Cysteine
                'E': 'Glu',  # Glutamic Acid
                'Q': 'Gln',  # Glutamine
                'G': 'Gly',  # Glycine
                'H': 'His',  # Histidine
                'I': 'Ile',  # Isoleucine
                'L': 'Leu',  # Leucine
                'K': 'Lys',  # Lysine
                'M': 'Met',  # Methionine
                'F': 'Phe',  # Phenylalanine
                'P': 'Pro',  # Proline
                'S': 'Ser',  # Serine
                'T': 'Thr',  # Threonine
                'W': 'Trp',  # Tryptophan
                'Y': 'Tyr',  # Tyrosine
                'V': 'Val',  # Valine
                '*': 'STOP'
            }
            amino_acids = [amino_acid_dict[i] for i in amino_acids][::-1]

            df.to_csv(Table_S1_codon_usage)

            plt.figure(figsize=(18, 15))
            ax = sns.heatmap(df, cmap='viridis', fmt=".2f", cbar_kws={'label': 'Codon Abundance % Changes'}, linewidths=0.5, linecolor='white')
            ax.set_xticks(range(len(string_dates)))
            ax.set_xticklabels(string_dates, ha="center")
            ax.set_yticks(range(len(df.index)))
            ax.set_yticklabels(df.index, va='center')
            pos, textvals = plt.yticks()
            plt.yticks(np.arange(len(pos)) + 0.5,textvals, 
                rotation=0, fontsize="14", va="center")
            pos, textvals = plt.xticks()
            plt.xticks(np.arange(len(pos)) + 0.5,textvals, 
                rotation=90, fontsize="16", ha="center")
            

            ax.collections[0].colorbar.ax.tick_params(labelsize=15)
            ax.figure.axes[-1].yaxis.label.set_size(20)


            # plt.xticks(rotation=45)

            #add groupings
            ax2 = ax.twinx()
            ax2.spines['left'].set_position(('axes', -0.05))
            ax2.tick_params('both', length = 0, which = 'minor', labelsize = 16)
            ax2.tick_params('both', direction = 'in', which = 'major', length = 12, labelsize = 16)
            ax2.yaxis.set_ticks_position('left')
            ax2.yaxis.set_label_position('left')

            amino_acids_final = []
            for i in amino_acids:
                if i not in amino_acids_final:
                    amino_acids_final.append(i)

            ticks = [0]
            locators = []

            for i in amino_acids_final:
                aa_num = amino_acids.count(i)
                pretick = aa_num/len(amino_acids)
                tick = ticks[-1] + pretick
                ticks.append(tick)
                locators.append(tick - (pretick / 2))

            ax2.set_yticks(ticks)
            ax2.yaxis.set_major_formatter(ticker.NullFormatter())
            ax2.yaxis.set_minor_locator(ticker.FixedLocator(locators))
            ax2.yaxis.set_minor_formatter(ticker.FixedFormatter(amino_acids_final))


            # Customize the plot
            plt.title('Codon Abundance % Changes Heatmap', fontsize=30)
            plt.xlabel('')
            plt.ylabel('')

            # Show the plot
            plt.tight_layout()
            plt.savefig(codon_heatmap_out)
            plt.clf()

            #GC HEATMAP
            df = pd.DataFrame.from_dict(gc_counts_dict, orient='index',
                                columns=string_dates)
            
            df.to_csv(Table_S4_gc_usage)
            # Set up the Seaborn heatmap
            plt.figure(figsize=(12, 8))
            ax = sns.heatmap(df, cmap='viridis', fmt=".2f", cbar_kws={'label': 'GC% Changes'}, linewidths=0.5, linecolor='white')
            ax.set_xticks(range(len(string_dates)))
            ax.set_xticklabels(string_dates, ha="center")

            pos, textvals = plt.yticks()
            ax.set_yticks(np.arange(len(df.index)) + 0.5)
            ax.set_yticklabels(df.index,  va="center")
            pos, textvals = plt.yticks()
            plt.yticks(np.arange(len(pos)) + 0.5,textvals, 
                rotation=45, fontsize="16", va="center")
            pos, textvals = plt.xticks()
            plt.xticks(np.arange(len(pos)) + 0.5,textvals, 
                rotation=90, fontsize="12", ha="center")
            

            ax.collections[0].colorbar.ax.tick_params(labelsize=15)
            ax.figure.axes[-1].yaxis.label.set_size(20)

            plt.title('GC Abundance % Changes Heatmap', fontsize = 30)
            plt.xlabel('')
            plt.ylabel('')

            # Show the plot
            plt.tight_layout()
            plt.savefig(gc_heatmap_out)
            plt.clf()

            #repeat for ENC
            df = pd.DataFrame.from_dict(enc_counts_dict, orient='index',
                                columns=string_dates)
            
            df.to_csv(Table_S6_enc)
            # Set up the Seaborn heatmap
            plt.figure(figsize=(12, 8))
            ax = sns.heatmap(df, cmap='viridis', fmt=".2f", cbar_kws={'label': 'ENC Changes'}, linewidths=0.5, linecolor='white')
            ax.set_xticks(range(len(string_dates)))
            ax.set_xticklabels(string_dates, ha="center", fontsize = 12)

            pos, textvals = plt.yticks()
            ax.set_yticks(np.arange(len(df.index)) + 0.5)
            ax.set_yticklabels(df.index,  va="center", fontsize = 14)
            pos, textvals = plt.yticks()
            plt.yticks(np.arange(len(pos)) + 0.5,textvals, 
                rotation=45, fontsize="12", va="center")
            pos, textvals = plt.xticks()
            plt.xticks(np.arange(len(pos)) + 0.5,textvals, 
                rotation=90, fontsize="12", ha="center")
            

            ax.collections[0].colorbar.ax.tick_params(labelsize=15)
            ax.figure.axes[-1].yaxis.label.set_size(20)


            # Customize the plot
            plt.title('ENC Abundance Heatmap',fontsize=20)
            plt.xlabel('Dates')
            plt.ylabel('')

            # Show the plot
            plt.tight_layout()
            plt.savefig(enc_heatmap_out)
            plt.clf()

            #repeat for dinucleotides
            df = pd.DataFrame.from_dict(dinucleotide_counts_dict, orient='index',columns=string_dates)
            
            df.to_csv(Table_S2_dn_usage)
            plt.figure(figsize=(12, 8))
            ax = sns.heatmap(df, cmap='viridis', fmt=".0f", cbar_kws={'label': 'Dinucleotide % Changes'}, linewidths=0.5, linecolor='white')
            ax.set_xticks(range(len(string_dates)))
            ax.set_xticklabels(string_dates, ha="center")

            pos, textvals = plt.yticks()
            ax.set_yticks(np.arange(len(df.index)) + 0.5)
            ax.set_yticklabels(df.index,  va="center")
            pos, textvals = plt.yticks()
            ax.set_yticks(np.arange(len(df.index)) + 0.5)
            ax.set_yticklabels(df.index,  va="center")
            pos, textvals = plt.yticks()
            plt.yticks(np.arange(len(pos)) + 0.5,textvals, 
                rotation=45, fontsize="16", va="center")
            pos, textvals = plt.xticks()
            plt.xticks(np.arange(len(pos)) + 0.5,textvals, 
                rotation=90, fontsize="12", ha="center")
            

            ax.collections[0].colorbar.ax.tick_params(labelsize=15)
            ax.figure.axes[-1].yaxis.label.set_size(20)

            # Customize the plot
            plt.title('Dinucleotide % Changes Heatmap', fontsize = 30)
            plt.xlabel('Dates')
            plt.ylabel('')

            # Show the plot
            plt.tight_layout()
            plt.savefig(dinucleotide_heatmap_out)
            plt.clf()

            #repeat for junction dinucleotides
            df = pd.DataFrame.from_dict(jdinucleotide_counts_dict, orient='index',columns=string_dates)
            
            df.to_csv(Table_S3_jdn_usage)
            plt.figure(figsize=(12, 8))
            ax = sns.heatmap(df, cmap='viridis', fmt=".0f", cbar_kws={'label': 'Junction Dinucleotide % Changes'}, linewidths=0.5, linecolor='white')
            ax.set_xticks(range(len(string_dates)))
            ax.set_xticklabels(string_dates, ha="center")

            pos, textvals = plt.yticks()
            ax.set_yticks(np.arange(len(df.index)) + 0.5)
            ax.set_yticklabels(df.index,  va="center")
            pos, textvals = plt.yticks()
            ax.set_yticks(np.arange(len(df.index)) + 0.5)
            ax.set_yticklabels(df.index,  va="center")
            pos, textvals = plt.yticks()
            plt.yticks(np.arange(len(pos)) + 0.5,textvals, 
                rotation=45, fontsize="16", va="center")
            pos, textvals = plt.xticks()
            plt.xticks(np.arange(len(pos)) + 0.5,textvals, 
                rotation=90, fontsize="12", ha="center")
            

            ax.collections[0].colorbar.ax.tick_params(labelsize=15)
            ax.figure.axes[-1].yaxis.label.set_size(20)

            # Customize the plot
            plt.title('Junction Dinucleotide % Changes Heatmap', fontsize = 30)
            plt.xlabel('Dates')
            plt.ylabel('')

            # Show the plot
            plt.tight_layout()
            plt.savefig(junc_heatmap_out)
            plt.clf()

            #CAI HEATMAP
            cmap = plt.cm.viridis
            cmap.set_over('red')
            df = pd.DataFrame.from_dict(cai_dict, orient='index',columns=string_dates)
            
            df.to_csv(Table_S5_cai)
            plt.figure(figsize=(12, 8))
            ax = sns.heatmap(df, fmt=".0f", cbar_kws={'label': 'CAI Changes'}, linewidths=0.5, linecolor='white', cmap = cmap, vmax = 0.01 )
            
            ax.set_xticks(range(len(string_dates)))
            ax.set_xticklabels(string_dates, ha="center", fontsize = 12)

            pos, textvals = plt.yticks()
            ax.set_yticks(np.arange(len(df.index)) + 0.5)
            ax.set_yticklabels(df.index,  va="center", fontsize = 14)
            pos, textvals = plt.yticks()
            plt.yticks(np.arange(len(pos)) + 0.5,textvals, 
                rotation=0, fontsize="12", va="center")
            pos, textvals = plt.xticks()
            plt.xticks(np.arange(len(pos)) + 0.5,textvals, 
                rotation=90, fontsize="12", ha="center")
            

            ax.collections[0].colorbar.ax.tick_params(labelsize=15)
            ax.figure.axes[-1].yaxis.label.set_size(20)

            # Customize the plot
            plt.title('CAI Changes Heatmap', fontsize = 20)
            plt.xlabel('Dates')
            plt.ylabel('')

            # Show the plot
            plt.tight_layout()
            plt.savefig(cai_heatmap_out)
            plt.clf()

            #CPAI HEATMAP
            df = pd.DataFrame.from_dict(cpai_dict, orient='index',columns=string_dates)
            df.to_csv(Table_S7_cpai)
            plt.figure(figsize=(12, 8))
            ax = sns.heatmap(df, cmap=cmap, fmt=".0f", cbar_kws={'label': 'CPAI Changes'}, linewidths=0.5, linecolor='white', vmax = 0.01)
            ax.set_xticks(range(len(string_dates)))
            ax.set_xticklabels(string_dates, ha="center", fontsize = 12)

            ax.set_yticklabels(df.index,  va="center")
            pos, textvals = plt.yticks()
            ax.set_yticks(np.arange(len(df.index)) + 0.5)
            ax.set_yticklabels(df.index,  va="center", fontsize = 14)
            pos, textvals = plt.yticks()
            plt.yticks(np.arange(len(pos)) + 0.5,textvals, 
                rotation=0, fontsize="12", va="center")
            pos, textvals = plt.xticks()
            plt.xticks(np.arange(len(pos)) + 0.5,textvals, 
                rotation=90, fontsize="12", ha="center")
            

            ax.collections[0].colorbar.ax.tick_params(labelsize=15)
            ax.figure.axes[-1].yaxis.label.set_size(20)

            # Customize the plot
            plt.title('CPAI Changes Heatmap', fontsize = 20)
            plt.xlabel('Dates')
            plt.ylabel('')

            # Show the plot
            plt.tight_layout()
            plt.savefig(cpai_heatmap_out)
            plt.clf()        
                
            #make vait pie chart
            vait_pie_charts(vait_pie_input_df, vait_pie_chart_path)
        
            