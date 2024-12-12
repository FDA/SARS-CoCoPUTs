from Bio import SeqIO
from collections import defaultdict
import subprocess
import os
import time
import sys

def write_unique_sequences(reference_fasta, input_fasta, output_fasta):
    # Dictionary to store unique sequences and their corresponding headers and counts
    with open(halign_nextclade_log, 'a') as halign_file:
        start_time = time.time()
        halign_file.write(f'Starting unique sequence collection for {input_fasta}\n')
        unique_sequences = defaultdict(lambda: {"header": None, "count": 0})

        # Read the input fasta file one record at a time
        with open(input_fasta, "r") as input_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                sequence_str = str(record.seq)
                location = record.description.split('|')[3]
                if 'USA' in location:

                    # If sequence is seen for the first time, store its header and set count to 1
                    if sequence_str not in unique_sequences:
                        unique_sequences[sequence_str]["header"] = record.description
                        unique_sequences[sequence_str]["count"] = 1
                    else:
                        # Increment count for duplicate sequences
                        unique_sequences[sequence_str]["count"] += 1

        sorted_sequences = sorted(unique_sequences.items(), key=lambda x: x[1]["count"], reverse=True)
        if len(sorted_sequences) != 0:
            finished_unique_seqs_time = time.time()
            halign_file.write(f'Finished unique sequence collection in {finished_unique_seqs_time - start_time} seconds, now writing to output file\n')
            with open(output_fasta, "w") as output_handle:
                with open(reference_fasta, 'r') as ref_handle:
                    for record in SeqIO.parse(ref_handle, 'fasta'):
                        ref_header = 'MN908947|referenceB|2019-12-01|Asia / China / Wuhan|Wuhan-Hu-1 '
                        ref_seq = str(record.seq)
                        output_handle.write(f">{ref_header}\n{ref_seq}\n")
            # Write unique sequences to the output fasta file
                for sequence, info in sorted_sequences:
                    # Create a modified header with the count added
                    new_header = f"{info['header']}|{info['count']}".replace(")","").replace("()","").replace(";","")
                    # Write the unique sequence to the output file in fasta format
                    output_handle.write(f">{new_header}\n{sequence}\n")
            final_finish_time = time.time()
            halign_file.write(f'Finished outfile writing in {final_finish_time - start_time} seconds\n')
        else:
            finished_unique_seqs_time = time.time()
            halign_file.write(f'No USA sequences found. No outfile created. Finished in {finished_unique_seqs_time - start_time}\n')
def run_halign(input_fasta_halign):
    if os.path.exists(input_fasta_halign):
        p = subprocess.run(f"halign -t 150 -Xmx600g {input_fasta_halign}", shell = True, capture_output = True)
        with open(halign_nextclade_log, 'a') as halign_file:
            halign_file.write(p.stdout.decode('UTF-8') + '\n')
    else:
        with open(halign_nextclade_log, 'a') as halign_file:
            halign_file.write(f'No input detected for running halign on {input_fasta_halign}, skipping.\n')
def run_nextclade(input_fasta_nextclade):
    if os.path.exists(input_fasta_nextclade):
        p = subprocess.run(f"nextclade run --input-dataset {nextclade_dataset_dir} --output-tsv {input_fasta_nextclade + '.tsv'} {input_fasta_nextclade}", shell = True, capture_output = True)
        with open(halign_nextclade_log, 'a') as halign_nextclade_file:
            halign_nextclade_file.write(p.stderr.decode('UTF-8') + '\n')
        os.remove(input_fasta_nextclade)
    else:
        with open(halign_nextclade_log, 'a') as halign_nextclade_file:
            halign_nextclade_file.write(f'No input detected for running nextclade on {input_fasta_nextclade}, skipping.\n')
    
if __name__ == '__main__':
    #First, this script reduces a given input fasta file down to only unique sequences (preserving the count of each 
    #sequence) to allow for faster computation
    #Then, this script generates multiple sequence alignments (for generation of our consensus sequences) and individual 
    #nextclade QC files for each fasta file, which is useful for our ancestral sequence generation (we use the mutation
    #data from these files, specifically)
    #inputs
    nextclade_dataset_dir = '../other_data/wuhan_dataset_nextclade/'
    reference_fasta = '../other_data/wuhan_dataset_nextclade/reference.fasta'
    input_dirs = [
        '../ncbiData/variant_split_files5',
        '../gisaidData/variant_split_files4',
        '../comparison_data/variant_split_files4'
    ]
    #outputs
    logs = [
        'ncbiData_unique_align_nextclade_USA.log',
        'gisaidData_unique_align_nextclade_USA.log',
        'comparison_data_unique_align_nextclade_USA.log'

    ]
    output_dirs = [
        '../ncbiData/aligned_unique_files_USA_only',
        '../gisaidData/aligned_unique_files_USA_only',
        '../comparison_data/aligned_unique_files_USA_only'
    ]
    for index, input_dir in enumerate(input_dirs):
        output_dir = output_dirs[index]
        halign_nextclade_log = logs[index]
        if os.path.exists(halign_nextclade_log):
            os.remove(halign_nextclade_log)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        for index, file in enumerate(os.listdir(input_dir)):
            if file.endswith('.fasta'):
                input_fasta = input_dir + '/' + file
                output_fasta = output_dir + '/' + file
                write_unique_sequences(reference_fasta, input_fasta, output_fasta)
                run_halign(output_fasta)
                run_nextclade(output_fasta)
