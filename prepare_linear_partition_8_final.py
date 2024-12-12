import xmltodict
import pandas as pd
import subprocess
import time
from collections import Counter
import numpy as np
import sys
import os

#the primary function of this script is to generate the "modified SHAPE" file,
#which is an optional input for linearPartition. We modify the SHAPE data from
#manfredonia et al so that all mutated positions in the alignment are redacted (along with their paired sites)
#this enables us to preserve experimentally informed secondary structure predictions at highly conserved regions
#while still allowing linearpartition software to make predictions at more variable regions
#the same shape file is used for each variant, modified only slightly to account for differences in length for particular
#variants

#inputs/outputs
#inputs
input_dirs = [
    'gisaidData',
    'comparison_data',
    'ncbiData'
]
input_seq_types = [
    'representative',
    'ancestral',
    'consensus'
]
input_syn_or_normals =[
    'syn_removed',
    'all_lineages_added_aligned'
]
for input_dir in input_dirs:
    for input_seq_type in input_seq_types:
        for input_syn_or_normal in input_syn_or_normals:

            paper_shape_data = '../other_data/SHAPE_invitro.xml'
            sequence_tsv = f'../{input_dir}/final_files_USA_only/{input_syn_or_normal}_{input_seq_type}_seqs.tsv'

            #outputs
            linear_partition_input_shape_file = '../other_data/reference_shape_file.txt'
            reference_sequence_tmp_file = '../other_data/wuhan1_reference_sequence.txt'
            reference_linear_partition_output_bpseq_style = '../other_data/wuhan1_ref_bpseq_struct_predict.txt'
            clustalo_input = '../other_data/clustalo_input.fa'
            clustalo_output = '../other_data/clustalo_output_refB.fa'
            final_shape_file = f'../other_data/USA_only_modified_shape_{input_dir}_{input_seq_type}_{input_syn_or_normal}.txt'

            #main code block below
            #parse shape data from paper
            with open(paper_shape_data, 'r', encoding='utf-8') as file:
                my_xml = file.read()
            my_dict = xmltodict.parse(my_xml)
            reference_sequence = [i for i in ''.join(my_dict['data']['transcript']['sequence'].split('\t')).replace('\n', '')]
            reactivity = ''.join(my_dict['data']['transcript']['reactivity'].split('\t')).replace('\n', '').split(',')
            shape_df = pd.DataFrame()
            shape_df['reactivity'] = reactivity
            shape_df['index'] = shape_df.index + 1


            #make reference shape csv for reading into linearpartition
            shape_df[['index','reactivity']].fillna('NA').to_csv(linear_partition_input_shape_file, sep='\t',header = None, index = False)


            #make tmpfile that contains reference sequence for reading into linearpartition
            with open(reference_sequence_tmp_file, 'w') as tmpfile:
                tmpfile.write(''.join(reference_sequence))

            #run linearPartition on reference using shape file
            p = subprocess.run(f'cat {reference_sequence_tmp_file} | linearpartition -V -M --bpseq --shape {linear_partition_input_shape_file}', capture_output=True, text = True, shell = True)

            #make a bpseq style structure file
            with open(reference_linear_partition_output_bpseq_style, 'w') as reffile:
                reffile.write(p.stdout)

            #read in your sequences
            sequence_df = pd.read_csv(sequence_tsv, sep = '\t')

            #figure out mutated positions
            wt_index = sequence_df[sequence_df['parent_lineages'] == 'referenceB'].index[0]
            position_df = sequence_df.sequence.str.split('', expand=True).iloc[:, 1:-1]
            all_same = position_df.nunique() == 1
            position_df = position_df.T
            position_df['all_same'] = all_same
            position_df = position_df[[wt_index, position_df.columns[-1]]]
            index = 1
            indices = []
            for i in position_df[position_df.columns[0]]:
                if i != '-':
                    indices.append(index)
                    index +=1
                else:
                    indices.append(np.nan)
            position_df['bp_only_index'] = indices
            position_df['bp_only_index'] = position_df['bp_only_index']

            # # resolve length discrepancy between gisaid wt and shape wt
            #This step is useful if your WT does not match the SHAPE data WT - previously we had used the GISAID reference sequence which did not match the NCBI
            #reference exactly (had a few extra "A"s at the polyA tail at the end, otherwise identical) - now for consistency we've switched to the NCBI reference
            #so this step isn't 100% necessary, but I've left it in incase we decide to switch references again in the future. 
            our_data_ref_seq = [i for i in list(sequence_df[sequence_df['parent_lineages'] == 'referenceB']['sequence'])[0].replace('-','')]
            with open(clustalo_input, 'w') as clustalofile:
                clustalofile.write(f">wuhan1_ref\n{''.join(reference_sequence)}\n")
                clustalofile.write(f">our_data_wuhan1_seq\n{''.join(our_data_ref_seq)}")
            p = subprocess.run(f'clustalo -i {clustalo_input} -o {clustalo_output} --force', shell = True)

            with open(clustalo_output, 'r') as alignmentfile:
                alignstring = alignmentfile.read()
            our_data_ref_seq = [i for i in alignstring.split('>')[2].split('\n',1)[-1].replace('\n', '').strip(' ')]

            #read in structure prediction of WT based on shape data
            struct_df = pd.read_csv(reference_linear_partition_output_bpseq_style, sep =' ', skiprows=1, header=None)
            shape_df['index'] = shape_df['index'].astype(float)
            shape_df['wt_aligned'] = our_data_ref_seq
            shape_df['paired_index'] = struct_df[2]
            shape_df = shape_df[shape_df['wt_aligned'] != '-']

            #merge everything together
            merge_df = position_df.merge(shape_df, left_on='bp_only_index', right_on='index',how='left')

            #replace both index and paired index with np.nan if the position is not 100% conserved
            for index, i in enumerate(merge_df['all_same']):
                if i == False and np.isnan(merge_df.loc[merge_df.index[index], 'index']) == False:
                    merge_df.loc[merge_df.index[index], 'reactivity'] = np.nan
                    paired_index = merge_df.loc[merge_df.index[index], 'paired_index'] -1
                    if paired_index != -1:
                        merge_df.loc[merge_df.index[paired_index], 'reactivity'] = np.nan
            merge_df['final_index'] = merge_df.index + 1
            merge_df[['final_index','reactivity']].fillna('NA').to_csv(final_shape_file,sep = '\t', header = False, index = False)


                