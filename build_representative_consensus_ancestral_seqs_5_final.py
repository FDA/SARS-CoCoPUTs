import os
from collections import defaultdict, Counter
import subprocess
import pandas as pd
import Levenshtein as lev
import numpy as np
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from Bio import SeqIO
import time

def timeit(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()  # Record the start time
        result = func(*args, **kwargs)  # Execute the function
        end_time = time.time()  # Record the end time
        print(f"Function '{func.__name__}' executed in {end_time - start_time:.4f} seconds")
        return result
    return wrapper
def get_wt_sequence():
    out = subprocess.run(f"grep -A1 'referenceB' {variant_file_dir}/B.fasta.aligned", shell = True, capture_output = True)
    out = out.stdout.decode('utf-8').split('\n')
    return out[1]

def get_consensus_sequence(fasta_file):
    # Initialize a defaultdict to keep track of nucleotide counts at each position
    start_time = time.time()
    print(f'Building consensus for {fasta_file}')
    position_nucleotide_counts = defaultdict(Counter)
    sequence_length = None

    # Parse the FASTA file one sequence at a time
    for index, record in enumerate(SeqIO.parse(fasta_file, "fasta")):
        #skip the first row because its always the wild type reference
        if index == 0:
            pass
        else:
            multiplier = int(record.description.split('|')[-1])
            seq = str(record.seq)
            if sequence_length is None:
                sequence_length = len(seq)  # Set the length of sequences based on the first sequence
            elif len(seq) != sequence_length:
                raise ValueError("All sequences must be of the same length.")
            # Update the nucleotide count for each position in the sequence
            for i, nucleotide in enumerate(seq):
                position_nucleotide_counts[i][nucleotide] += multiplier
        if (index +1) % 100000 == 0:
            print(f'Processed {index +1} for {fasta_file} in {time.time() - start_time} seconds')


    # Build the consensus sequence
    consensus_seq = []
    for i in range(sequence_length):
        # Get the most common nucleotide for each position
        consensus_nucleotide = position_nucleotide_counts[i].most_common(1)[0][0]
        consensus_seq.append(consensus_nucleotide)
    print(f'Finished consensus for {fasta_file}')

    # Join the list of nucleotides into a string and return
    return ''.join(consensus_seq).replace('-','')


def get_representative_seq(fasta_file):
    max_count = 0
    potential_rep_seqs = []
    total_count = 0
    for index, record in enumerate(SeqIO.parse(fasta_file, "fasta")):
        if index == 0:
            wt_seq = str(record.seq)
        elif index == 1:
            max_count = int(record.description.split('|')[-1])
            potential_rep_seqs.append([record.description, str(record.seq)])
            total_count += int(max_count)
        else:
            count = int(record.description.split('|')[-1])
            if count < max_count:
                total_count += count
                
            else:
                total_count += count
                potential_rep_seqs.append([record.description, str(record.seq)])
    if len(potential_rep_seqs)==1:
        max_count_seq = potential_rep_seqs[0]
        max_count_seq[0] = f'{max_count_seq[0]}-{str(total_count)}-{str(max_count / total_count)}'
        return max_count_seq
    else:
        lev_distances = [lev.distance(wt_seq, i[0]) for i in potential_rep_seqs]
        max_count_seq = potential_rep_seqs[np.argmin(lev_distances)]
        max_count_seq[0] = f'{max_count_seq[0]}-{str(total_count)}-{str(max_count / total_count)}'
        return max_count_seq

def identify_percent_mut_found_in_all_seqs(test_muts,all_muts_dict, num_seqs):  
    test_mut_representativeness = np.mean([all_muts_dict[i] / num_seqs for i in test_muts.split(',') if i != ''])
    return test_mut_representativeness
def adjust_missing_col(missing_vals):
    if missing_vals == '':
        return ''
    else:
        return ','.join(['missing_'+ i for i in missing_vals.split(',')])
def convert_deletion_to_individual_muts(deletion_val):
    if deletion_val != '':
        deletions = deletion_val.split(',')
        
        all_deleted_nts = []
        for deletion in deletions:
            if '-' in deletion:
                start = int(deletion.split('-')[0])
                end = int(deletion.split('-')[1])
                deleted_nts = [str(i) +'_del' for i in range(start,end + 1)]
                all_deleted_nts += deleted_nts
            else:
                all_deleted_nts.append(deletion + '_del')
                

        return ','.join(all_deleted_nts)
    else:
        return ''

def add_deletions_due_to_missing_coverage(start_end):
    alignment_start = int(start_end.split('_')[0])
    alignment_end = int(start_end.split('_')[1])
    sars_cds_start = 266
    sars_cds_end = 29674
    missing_from_start = alignment_start - sars_cds_start
    missing_from_end = sars_cds_end - alignment_end
    deletions = []
    if missing_from_start > 0:
        for i in range(sars_cds_start, alignment_start + 1):
            deletions.append(str(i) + '_del')
    if missing_from_end > 0:
        for i in range(alignment_end, sars_cds_end + 1):
            deletions.append(str(i) + '_del')
    return ','.join(deletions)

def get_ancestral_seq(nextclade_out_tsv_path, fasta_file_path):
    df = pd.read_csv(nextclade_out_tsv_path, sep = '\t', dtype=str)
    df = df[~df['seqName'].str.contains('referenceB')]
    df['mut_counts'] =  df[['totalSubstitutions',
                        'totalDeletions',	
                        'totalInsertions',	
                        'totalFrameShifts',	
                        'totalMissing',	
                        'totalNonACGTNs']].astype(int).sum(axis=1)    

    df = df.sort_values(by=['mut_counts', 'coverage'], ascending=[True, False])
    df['missing'].fillna('',inplace=True)
    df['missing'] = df['missing'].apply(adjust_missing_col)
    df['deletions'] = df['deletions'].fillna('').apply(convert_deletion_to_individual_muts)
    df['start_end'] = df['alignmentStart'] + '_' + df['alignmentEnd']
    df['start_end'] = df['start_end'].apply(add_deletions_due_to_missing_coverage)
    

    df['muts'] = df['substitutions'].fillna('') + ','+df['deletions']+ ','+df['insertions'].fillna('')+ ','+df['missing'].fillna('') + ','+df['start_end'].fillna('')

    num_seqs = len(df.index)
    all_muts_dict = Counter([i for i in ','.join(list(df['muts'])).split(',') if i != ''])
    df['representativeness'] = df['muts'].apply(identify_percent_mut_found_in_all_seqs, args = (all_muts_dict, num_seqs))
    df['representativeness'].fillna(1,inplace = True)
    best_seq = df.loc[df['representativeness'].idxmax()]['seqName']
    p = subprocess.run(f'grep -A 1 "{best_seq}" {fasta_file_path}', capture_output=True, shell = True)
    sequence = p.stdout.decode('utf-8').strip('\n').split('\n')[-1]
    return [best_seq, sequence]

def write_to_tsv(seqfilepath):
    repSeqID, repSeq = get_representative_seq(seqfilepath)
    conseq = get_consensus_sequence(seqfilepath)
    anseqID, anseq = get_ancestral_seq(seqfilepath.replace('.aligned','.tsv'), seqfilepath)
    with lock:
        with open(out_representative_seq_tsv_file,'a') as rep_tsv_file:
            split_ID = repSeqID.split('|')
            rep_tsv_file.write(f"{split_ID[1]}\t{split_ID[-1].split('-')[0]}\t{split_ID[-1].split('-')[1]}\t{split_ID[-1].split('-')[2]}\t{repSeq}\t{split_ID[0].strip('>')}\n")
        with open(out_consensus_seq_tsv_file, 'a') as cons_tsv_file:
            cons_tsv_file.write(f"{seqfilepath.replace('.fasta.aligned','').split('/')[-1]}\t{'NaN'}\t{'NaN'}\t{'NaN'}\t{conseq}\t{'NaN'}\n")
        with open(out_ancestral_seq_tsv_file, 'a') as ans_tsv_file:
            split_ID = anseqID.split('|')
            print(f'Writing ancestral sequence for {seqfilepath}')
            ans_tsv_file.write(f"{split_ID[1]}\t{split_ID[-1]}\t{'NaN'}\t{'NaN'}\t{anseq}\t{split_ID[0].strip('>')}\n")           
def run_multiprocessing(function, input_list, workers):
    workers = workers
    with ProcessPoolExecutor(workers) as executor:
        results = executor.map(function, input_list)
        for i in results:
            pass

if __name__ == '__main__':
    #This script builds/collects three types of sequences for each pango lineage, which are used for downstream analyses
    #the reasoning/methodology is explained within the text of our manuscript
    #filepaths
    #inputs
    variant_file_dir = '../ncbiData/aligned_unique_files_USA_only/'
    #outputs
    out_representative_seq_tsv_file = '../ncbiData/final_files_USA_only/representative_seqs.tsv'
    out_consensus_seq_tsv_file = '../ncbiData/final_files_USA_only/consensus_seqs.tsv'
    out_ancestral_seq_tsv_file = '../ncbiData/final_files_USA_only/ancestral_seqs.tsv'

    if not os.path.exists(os.path.dirname(out_representative_seq_tsv_file)):
        os.makedirs(os.path.dirname(out_representative_seq_tsv_file))
    if not os.path.exists(os.path.dirname(out_consensus_seq_tsv_file)):
        os.makedirs(os.path.dirname(out_consensus_seq_tsv_file))
    if not os.path.exists(os.path.dirname(out_ancestral_seq_tsv_file)):
        os.makedirs(os.path.dirname(out_ancestral_seq_tsv_file))
    
    seq_filepaths = [variant_file_dir + i for i in os.listdir(variant_file_dir) if i.endswith('aligned')]
    wt_seq = get_wt_sequence()
    lock = multiprocessing.Lock()

    with open(out_representative_seq_tsv_file, 'w') as outtsvfile:
        outtsvfile.write('terminal_lineage\tmax_count\ttotal_count\tproportion\tsequence\tID\n')
    with open(out_consensus_seq_tsv_file, 'w') as outtsvfile:
        outtsvfile.write('terminal_lineage\tmax_count\ttotal_count\tproportion\tsequence\tID\n')
    with open(out_ancestral_seq_tsv_file, 'w') as outtsvfile:
        outtsvfile.write('terminal_lineage\tmax_count\ttotal_count\tproportion\tsequence\tID\n')

    run_multiprocessing(write_to_tsv, seq_filepaths,100)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
    
    with open(out_representative_seq_tsv_file,'a') as outtsvfile:
        outtsvfile.write(f"referenceB\tNaN\tNaN\tNaN\t{wt_seq}\tMN908947\n")
    with open(out_consensus_seq_tsv_file,'a') as outtsvfile:
        outtsvfile.write(f"referenceB\tNaN\tNaN\tNaN\t{wt_seq}\tMN908947\n")
    with open(out_ancestral_seq_tsv_file,'a') as outtsvfile:
        outtsvfile.write(f"referenceB\tNaN\tNaN\tNaN\t{wt_seq}\tMN908947\n")