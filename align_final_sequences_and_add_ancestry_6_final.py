import pandas as pd
import subprocess
from Bio import SeqIO
import os
import sys
from pango_aliasor.aliasor import Aliasor
import numpy as np

def align_final_seqs(input_tsv):
    df = pd.read_csv(input_tsv,sep = '\t')
    seqs = list(df['sequence'].str.replace('-',''))
    variant_names = list(df['terminal_lineage'])
    #make temp fasta
    with open (input_tsv + '.tmp.fasta', 'w') as tmp_fasta_file:
        for index, variant in enumerate(variant_names):
            tmp_fasta_file.write('>'+variant + '\n' + seqs[index] + '\n')
    subprocess.run(f"halign -t 150 -Xmx600g {input_tsv + '.tmp.fasta'}", shell = True, capture_output=True)
    aligned_seqs = [record.seq for record in SeqIO.parse(input_tsv + '.tmp.fasta.aligned', "fasta")]

    df['sequence'] = aligned_seqs

    dir_name = os.path.dirname(input_tsv)
    base_name = os.path.basename(input_tsv)
    aligned_name = f"aligned_{base_name}"
    aligned_filename = os.path.join(dir_name, aligned_name)

    df.to_csv(aligned_filename, sep = '\t', index = False)
    os.remove(input_tsv + '.tmp.fasta')
    os.remove(input_tsv + '.tmp.fasta.aligned')
    return aligned_filename

def add_ancestry(input_csv_name):
    aliasor = Aliasor()
    df = pd.read_csv(input_csv_name, sep = '\t')

    ref_df = pd.DataFrame()
    ref_df['parent_lineages'] = df[df['terminal_lineage'].str.contains('reference')]['terminal_lineage'].str.split('_').str[0]
    ref_df['terminal_lineage'] = df[df['terminal_lineage'].str.contains('reference')]['terminal_lineage']
    ref_df['max_count'] = np.nan
    ref_df['total_count'] = np.nan
    ref_df['proportion'] = np.nan
    ref_df['sequence'] = df[df['terminal_lineage'].str.contains('reference')]['sequence']
    ref_df['ID'] = df[df['terminal_lineage'].str.contains('reference')][f'ID']
    

    df = df[~df['terminal_lineage'].str.contains('reference')]

    df.sort_values(by=[df.columns[0]], inplace=True)
    lineages =df[df.columns[0]].tolist()
    all_parents = []
    all_parents_lengths = []
    for index, lineage in enumerate(lineages):
        stripped_lineage = lineage.split('_')[0]

        uncompress_lineage = aliasor.uncompress(stripped_lineage)
        uncompress_lineage_split_list = uncompress_lineage.split('.')

        parents = []
        for i in range(len(uncompress_lineage_split_list)):
            compress_parent_lineage = aliasor.compress('.'.join(uncompress_lineage_split_list[0:i+1]))
            parents.append(compress_parent_lineage)
        parents.append(lineage)
        all_parents.append(parents)
        all_parents_lengths.append(len(parents))

    max_length = max(all_parents_lengths)
    all_parents_equal = []
    for index, i in enumerate(all_parents):
        parent_length = len(i)
        difference = max_length - parent_length
        diff_list = [''] * difference
        diff_list.extend(i)
        all_parents_equal.append(diff_list)

    df_parents = pd.DataFrame(all_parents_equal)
    df_parents['parent_lineages'] = df_parents[df_parents.columns[:-1]].astype(str).agg('/'.join, axis=1).str.strip('/')
    # print(df_parents)
    df_parents['terminal_lineage'] = df_parents[df_parents.columns[-2]]
    df_parents = df_parents[['parent_lineages', 'terminal_lineage']]

    df = df.rename({df.columns[0]:'terminal_lineage'}, axis=1)


    merge_df = pd.merge(df_parents, df, on='terminal_lineage', how='left')
    merge_df = merge_df[~merge_df['terminal_lineage'].str.contains('no_lineage')]
    merge_df = merge_df[~merge_df['terminal_lineage'].str.contains('unclassifiable')]
    merge_df = pd.concat([merge_df, ref_df])
    merge_df.insert(2, 'gene','whole_seq')

    dir_name = os.path.dirname(input_csv_name)
    base_name = os.path.basename(input_csv_name)
    out_base_name = f"all_lineages_added_{base_name}"
    out_filename = os.path.join(dir_name, out_base_name)

    merge_df.to_csv(out_filename, index=False, sep = '\t')

if __name__=='__main__':
    #this script takes the output from the consensus/ancestral/representative sequence building script 
    #and adds variant lineage information, as well as aligning all of the sequences. 
    input_tsvs = [
        '../gisaidData/final_files_USA_only/ancestral_seqs.tsv',
        '../gisaidData/final_files_USA_only/consensus_seqs.tsv',
        '../gisaidData/final_files_USA_only/representative_seqs.tsv',
        '../ncbiData/final_files_USA_only/ancestral_seqs.tsv',
        '../ncbiData/final_files_USA_only/consensus_seqs.tsv',
        '../ncbiData/final_files_USA_only/representative_seqs.tsv',
        '../comparison_data/final_files_USA_only/ancestral_seqs.tsv',
        '../comparison_data/final_files_USA_only/consensus_seqs.tsv',
        '../comparison_data/final_files_USA_only/representative_seqs.tsv'        
    ]
    aligned_tsvs = [align_final_seqs(i) for i in input_tsvs]
    
    for i in aligned_tsvs:
        add_ancestry(i)
