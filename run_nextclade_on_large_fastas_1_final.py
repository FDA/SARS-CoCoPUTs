import os
from Bio import SeqIO
import os
from Bio import SeqIO
import subprocess

def split_fasta_streaming(input_fasta, output_dir, output_prefix="split", sequences_per_file=1000000):
    """
    Splits a large FASTA file into smaller FASTA files, each containing a specified number of sequences.
    Writes each sequence to the output file immediately as it is parsed.

    Parameters:
    - input_fasta (str): Path to the input FASTA file.
    - output_dir (str): Directory where the split FASTA files will be saved.
    - output_prefix (str): Prefix for the output FASTA files. Defaults to "split".
    - sequences_per_file (int): Number of sequences per output file. Defaults to 1,000,000.

    Raises:
    - FileNotFoundError: If the input FASTA file does not exist.
    - IOError: If there is an error writing to the output files.
    """
    if not os.path.isfile(input_fasta):
        raise FileNotFoundError(f"The file '{input_fasta}' does not exist.")

    os.makedirs(output_dir, exist_ok=True)

    file_index = 1
    sequence_count = 0
    total_sequences = 0
    basename = os.path.basename(input_fasta)
    output_file_path = os.path.join(output_dir, f"{output_prefix}_{file_index}_{basename}")
    try:
        out_handle = open(output_file_path, 'w')
    except IOError as e:
        raise IOError(f"Error opening output file '{output_file_path}': {e}")

    print(f"Creating {output_file_path}...")

    try:
        with open(input_fasta, 'r') as infile:
            fasta_parser = SeqIO.parse(infile, 'fasta')
            for record in fasta_parser:
                try:
                    SeqIO.write(record, out_handle, 'fasta')
                except IOError as e:
                    out_handle.close()
                    raise IOError(f"Error writing to output file '{output_file_path}': {e}")

                sequence_count += 1
                total_sequences += 1

                if sequence_count >= sequences_per_file:
                    out_handle.close()
                    print(f"Created {output_file_path} with {sequence_count} sequences.")

                    file_index += 1
                    sequence_count = 0
                    output_file_path = os.path.join(output_dir, f"{output_prefix}_{file_index}_{basename}")
                    try:
                        out_handle = open(output_file_path, 'w')
                        print(f"Creating {output_file_path}...")
                    except IOError as e:
                        raise IOError(f"Error opening output file '{output_file_path}': {e}")

        # After the loop, close the last output file if it's still open
        if not out_handle.closed:
            out_handle.close()
            print(f"Created {output_file_path} with {sequence_count} sequences.")

        print(f"Finished splitting. Total sequences processed: {total_sequences}. Total files created: {file_index}")

    except Exception as e:
        if not out_handle.closed:
            out_handle.close()
        raise e
        

def run_nextclade(input_fasta, nextclade_dataset):
    outfile = input_fasta.replace('.fasta','').replace('.fa','') + '_nextclade_qc.tsv'
    print('Generating: ',outfile)
    subprocess.run(f'nextclade run --input-dataset {nextclade_dataset} --output-tsv {outfile} {input_fasta}', shell = True)
    
if __name__ == '__main__':
    #this script runs nextclade on the initial raw data from both gisaid and ncbi
    #it will generate an output tsv from nextclade that contains QC information for all sequences
    #this will be used in subsequent scripts
    ncbi_raw_fasta = '../ncbiData/rawData/sars_cov_2_ncbi_sequences.fa'
    ncbi_outdir = os.path.dirname(ncbi_raw_fasta)
    gisaid_raw_fasta = '../gisaidData/rawData/sars_cov_2_gisaid_sequences.fasta'
    gisaid_outdir = os.path.dirname(gisaid_raw_fasta)
    split_file_prefix = "split"
    nextclade_dataset = ' ../other_data/wuhan_dataset_nextclade/'
    
    #first we must split the large fastas into smaller ones because nextclade crashes when the input fasta is too large
    split_fasta_streaming(ncbi_raw_fasta, ncbi_outdir, output_prefix="split", sequences_per_file=1000000)
    split_fasta_streaming(gisaid_raw_fasta, gisaid_outdir, output_prefix="split", sequences_per_file=1000000)
    
    #run nextclade on the splitfiles
    # ncbi
    with os.scandir(ncbi_outdir) as entries:
        split_files = [
            os.path.join(ncbi_outdir,entry.name) for entry in entries
            if entry.is_file() and entry.name.startswith(split_file_prefix) and bool(entry.name.endswith(".fa") or entry.name.endswith(".fasta"))
        ]
        
    for split_fasta in split_files:
        run_nextclade(split_fasta, nextclade_dataset)

    #gisaid
    with os.scandir(gisaid_outdir) as entries:
        split_files = [
            os.path.join(gisaid_outdir,entry.name) for entry in entries
            if entry.is_file() and entry.name.startswith(split_file_prefix) and bool(entry.name.endswith(".fa") or entry.name.endswith(".fasta"))
        ]
        
    for split_fasta in split_files:
        run_nextclade(split_fasta, nextclade_dataset)
    
    #combine split files
    subprocess.run(f"awk 'FNR==1 && NR!=1 {{ next }} {{ print }}' {ncbi_outdir}/{split_file_prefix}_*.tsv > {ncbi_outdir}/sars_cov_2_ncbi_sequences_nextclade_qc.tsv", shell = True)
    subprocess.run(f"awk 'FNR==1 && NR!=1 {{ next }} {{ print }}' {gisaid_outdir}/{split_file_prefix}_*.tsv > {gisaid_outdir}/sars_cov_2_gisaid_sequences_nextclade_qc.tsv", shell = True)