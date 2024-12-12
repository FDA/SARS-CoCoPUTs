import pandas as pd
import threading
from concurrent.futures import ThreadPoolExecutor
import time
import concurrent.futures
import subprocess
import sys
import pandas as pd
from collections import Counter
import re
import os

def timeit(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()  # Record the start time
        result = func(*args, **kwargs)  # Execute the function
        end_time = time.time()  # Record the end time
        print(f"Function '{func.__name__}' executed in {end_time - start_time:.4f} seconds")
        return result
    return wrapper
@timeit
def calculate_mea_struct_and_forgi(seq_and_name):
    search_out = 'nothing'
    try:
        name = seq_and_name[0]
        term_lineage = name.split('/')[-1]
        print('Processing: ', name)
        seq = seq_and_name[-2]
        #check if already complete
        p = subprocess.run(f'grep -P "\t{term_lineage}\t" {outfile_tmp}', capture_output=True, text = True, shell = True)
        search_out = p.stdout
        if search_out == '':
            new_shape = modified_shape.copy()
            seq_list = list(seq)
            new_shape['testseq'] = seq_list
            new_shape = new_shape[new_shape['testseq'] != '-']
            new_shape.reset_index(inplace = True)
            new_shape['index'] = new_shape.index + 1
            shapefilepath = f'{variant_shape_file_dir}{name.replace("/","_")}_shape.txt'
            tmpfilename = 'tmp_' + name.replace("/","_")
            with open(tmpfilename, 'w') as tmpfile:
                tmpfile.write(''.join(list(new_shape['testseq'])))
            new_shape[['index',1]].fillna('NA').to_csv(shapefilepath,sep = '\t', header = None, index = False)
            p = subprocess.run(f'cat {tmpfilename} | linearpartition -V -M --shape {shapefilepath}', capture_output=True, text = True, shell = True)

            structure = p.stdout.split('\n')[2]
            free_energy = p.stderr.split(' ')[4]
            with open(tmpfilename, 'w') as tmpfile:
                tmpfile.write(structure)
            p = subprocess.run(f'conda run -n forgi python {forgi_rnaConvert_path} {tmpfilename} -T element_string', 
                                capture_output = True, 
                                text = True, 
                                shell = True)
            forgi_output = p.stdout.split('\n')[1]

            subprocess.run(f'rm {tmpfilename}', shell = True)
            outputstring = '\t'.join([str(i) for i in seq_and_name]) + '\t' + structure + '\t'+free_energy+'\t'+forgi_output+'\n'
    
    except:
        print('ERROR ENCOUNTERED: ', seq_and_name[0])
        outputstring = '\t'.join([str(i) for i in seq_and_name]) + '\t' + 'Error' +'\t' +'Error' +'\t' +'Error' + '\n'
    if search_out == '':
        Lock.acquire()
        with open(outfile_tmp, 'a') as finaloutputfile:
            finaloutputfile.write(outputstring)
        Lock.release()
        return f'Completed {name}'
    else:
        return (f'{name} already processed')



if __name__ == '__main__':
    #please note - as is, this script takes several days to run on the ~18 different datasets here (9x9x2)
    #this script generates both the dot-bracket secondary structures as well as the forgi strings for our
    #representative/consensus/ancestral sequences. All of this is saved to output tsvs.
    #resume or start fresh
    resume_or_start = 'start'

    input_dirs = [
    'comparison_data',
    'gisaidData',
    'ncbiData'
    ]
    input_seq_types = [
        'ancestral',
        'consensus',
        'representative'
    ]
    input_syn_or_normals =[
        'all_lineages_added_aligned',
        'syn_removed'
        
    ]
    
    for input_seq_type in input_seq_types:
        for input_dir in input_dirs:
            for input_syn_or_normal in input_syn_or_normals:
                #inputs/outputs
                #inputs
                modified_shape_file = f'../other_data/USA_only_modified_shape_{input_dir}_{input_seq_type}_{input_syn_or_normal}.txt'
                sequence_file = f'../{input_dir}/final_files_USA_only/{input_syn_or_normal}_{input_seq_type}_seqs.tsv'
                forgi_rnaConvert_path = '../forgi-master/examples/rnaConvert.py'

                #outputs
                variant_shape_file_dir = f'../{input_dir}/final_files_USA_only/variant_shapes_{input_seq_type}_{input_syn_or_normal}/'
                outfile_tmp = f'../{input_dir}/final_files_USA_only/{input_syn_or_normal}_{input_seq_type}_seqs_secondary_structure_file.tsv'
                outfile_final = f'../{input_dir}/final_files_USA_only/{input_syn_or_normal}_{input_seq_type}_seqs_final_secondary_structure_file.tsv'
                
                if os.path.exists(outfile_final):
                    print(f'{outfile_final} already generated, skipping')
                    pass
                else:
                    print(f'{outfile_final} not yet generated')

                    if not os.path.exists(variant_shape_file_dir):
                        os.makedirs(variant_shape_file_dir)
                    print(f'Now starting generation of {outfile_final} - will take several hours')
                    #delete tmp files if they exist already
                    subprocess.run('rm tmp_*', shell = True)
                    
                    #read in modified shape data prepared by prepare_linear_partition_5.py
                    modified_shape = pd.read_csv(modified_shape_file, sep = '\t', header = None)
                    seqfile= pd.read_csv(sequence_file,sep='\t')
                    names_and_genomes = list(zip(*[seqfile[seqfile.columns[index]] for index, i in enumerate(seqfile.columns)]))
                    
                    #create output file
                    if resume_or_start == 'start':
                        with open(outfile_tmp, 'w') as finaloutputfile:
                            finaloutputfile.write('\t'.join(seqfile.columns) + '\t' + 'linearpartition_predict' + '\t'+'ensemble_free_energy'+'\t'+'forgi_output\n')
                    
                    # run linear partition with multithreading
                    Lock = threading.Lock()
                    workers = 75
                    total = 0
                    with ThreadPoolExecutor(workers) as executor:
                        # print([i for i in range(len(transcript_set))])

                        results = executor.map(calculate_mea_struct_and_forgi, names_and_genomes)
                        for i in results:
                            print(i)

                    #add some calculations - # of nucleotides per secondary structure type
                    struct_df = pd.read_csv(outfile_tmp, sep = '\t')
                    def count_nucleotides(seq):
                        return dict(Counter(seq))

                    struct_df['count_nt'] = struct_df['forgi_output'].apply(count_nucleotides)

                    struct_df = pd.concat([struct_df, struct_df['count_nt'].apply(pd.Series)],axis = 1)

                    # struct_df['terminal_lineage'] = struct_df['parent_lineages'].str.split('/').str[-1]

                    forgi_strings = list(struct_df['forgi_output'])
                    linearpredict_strings = list(struct_df['linearpartition_predict'])
                    msa_strings = list(struct_df['sequence'])
                    new_forgi_strings = []

                    for index, msa_string in enumerate(msa_strings):
                        new_forgi_string = ''
                        forgi_index = 0
                        print('completed: ', index + 1)
                        for index2, i in enumerate(msa_string):
                            if i != '-':
                                new_forgi_string += forgi_strings[index][forgi_index]
                                forgi_index +=1
                            else:
                                new_forgi_string += '-'
                        new_forgi_strings.append(new_forgi_string)

                    struct_df['msa_forgi_output'] = new_forgi_strings
                    new_linearpredict_strings = []
                    for index, msa_string in enumerate(msa_strings):
                        new_forgi_string = ''
                        forgi_index = 0
                        print('completed: ', index + 1)
                        for index2, i in enumerate(msa_string):
                            if i != '-':
                                new_forgi_string += linearpredict_strings[index][forgi_index]
                                forgi_index +=1
                            else:
                                new_forgi_string += '-'
                        new_linearpredict_strings.append(new_forgi_string)

                    struct_df['msa_linearpartition_predict'] = new_linearpredict_strings
                    struct_df.drop('count_nt', inplace = True, axis = 1)

                    struct_df.to_csv(outfile_final, sep = '\t', index = False)