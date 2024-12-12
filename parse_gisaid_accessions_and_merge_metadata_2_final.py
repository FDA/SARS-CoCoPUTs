import dask.dataframe as dd
from dask.distributed import Client,  wait
import sqlite3
from Bio import SeqIO
import time
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import subprocess
import re
import hashlib
from filelock import FileLock

def timeit(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()  # Record the start time
        result = func(*args, **kwargs)  # Execute the function
        end_time = time.time()  # Record the end time
        print(f"Function '{func.__name__}' executed in {end_time - start_time:.4f} seconds")
        return result
    return wrapper
def fasta_to_dict(file_path):
    with open(file_path, 'r') as file:
        fasta_content = file.read()
    split_file = [i.split('\n', 1) for i in fasta_content.split('>')][1:]
    headers =  [i[0] for i in split_file]
    seqs = [i[1].replace('\n','') for i in split_file]
    return dict(zip(headers, seqs))
def cut_string(input_string):
    pattern = r'^(?:[^/]*/){3}[^|]*\|'
    match = re.search(pattern, input_string)
    if match:
        return input_string[:match.end()]
    return input_string
@timeit
def create_db_from_fasta_with_hash(db_name, table_name, fasta_file):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    c.execute(f'''
        CREATE TABLE IF NOT EXISTS {table_name} (
            header TEXT,
            sequence TEXT,
            sequence_hash TEXT
        )
    ''')
    conn.commit()
    index = 0
    with open('invalid_id_file_gisaid.txt', 'w') as invalid_id_file:
        with open(fasta_file, "r") as fasta:
            for record in SeqIO.parse(fasta, "fasta"):
                try:
                    header = cut_string(record.description).rstrip('|')
                    sequence = str(record.seq).replace('-', '')

                    # Create a SHA-256 hash of the sequence
                    sequence_hash = hashlib.sha256(sequence.encode('utf-8')).hexdigest()

                    c.execute('INSERT INTO sequences (header, sequence, sequence_hash) VALUES (?, ?, ?)', (header, sequence, sequence_hash))
                except:
                    print(f"Error for header {header} skipped.")
                    invalid_id_file.write('Invalid ID: ' + str(record.description) + '\n')

                index += 1

    print('Finished adding sequences to db, now making indices')
    c.execute("PRAGMA cache_size = -2000000")  # Set cache size to 2GB (in pages, assuming page size is 1024 bytes)
    c.execute("PRAGMA journal_mode = WAL")     # Use Write-Ahead Logging for better concurrency
    c.execute('CREATE INDEX IF NOT EXISTS idx_header ON sequences (header)')
    
    print('accession index finished, now making index on sequence hash column')
    c.execute('CREATE INDEX IF NOT EXISTS idx_sequence_hash ON sequences (sequence_hash)')
    conn.commit()
    conn.close()    
def query_sequence(db_name, header):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    c.execute('SELECT sequence FROM sequences WHERE header = ?', (header,))
    result = c.fetchone()
    conn.close()
    if result:
        return result[0]
    else:
        return "Header not found"
def query_nextclade(db_name, virus_name):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    c.execute('SELECT * FROM nextclade_qc WHERE "Virus name" = ?', (virus_name,))
    result = c.fetchone()
    conn.close()
    if result:
        return result
    else:
        return "Header not found"
def get_all_headers(db_name):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    
    # SQL query to select all headers
    c.execute('''SELECT header FROM sequences''')
    
    # Fetch all results
    headers = c.fetchall()
    
    # Close the connection
    conn.close()
    
    # Convert list of tuples to a list of strings
    header_list = [header[0] for header in headers]
    return header_list


def is_valid_sequence_regex(sequence):
    return bool(re.fullmatch("[ACTG-]+", sequence))
def count_non_actg_chars(sequence):
    """
    Count the number of characters in a string that are not 'A', 'C', 'T', or 'G'.

    :param sequence: A string representing the sequence.
    :return: The count of non-ACTG characters.
    """
    valid_chars = {'A', 'C', 'T', 'G'}
    return sum(1 for char in sequence if char.upper() not in valid_chars)
def pair_gisaid_accession_with_variant(virus_name):
    try:
        nextclade_result = query_nextclade(gisaid_nextclade_qc_sqldb, virus_name)
        if nextclade_result[10] == 'good':
            query_result = query_metadata_db(virus_name)
            accession, completeness, lineage, ncontent,collectionDate,location = [query_result[i] for i in [0,19,13,22,5,6]]

            try:
                ncontent = float(ncontent)
            except:
                ncontent = 0
            completeness = bool(completeness)
            if completeness and ncontent < 0.01:
                sequence = query_sequence(gisaid_fasta_sqldb, virus_name)

                new_header = f'>{accession}|{lineage}|{collectionDate}|{location}|{virus_name}'
                # with lock:
                lock = FileLock(f'{variant_split_files_path}/{lineage}.fasta.lock')
                with lock:
                    with open(f'{variant_split_files_path}/{lineage}.fasta', 'a') as variant_file:
                        variant_file.write(f'{new_header}\n{sequence}\n')
                lock = FileLock(qc_success_path + '.lock')
    except:
        lock = FileLock(invalid_headers_path + '.lock')
        with lock:
            invalid_id_file.write(virus_name + '\n')


def process_chunk(chunk, func, num_threads=None):
    """
    This function runs in a separate process and uses threads to process a chunk of data.

    Parameters:
    - chunk: The subset of data to process.
    - func: The function to apply to each item in the chunk.
    - num_threads: Optional; The number of threads to use in this process. Defaults to the number of CPU cores.

    Returns:
    - A list of results obtained by applying 'func' to each item in the chunk.
    """
    if num_threads is None:
        num_threads = multiprocessing.cpu_count()
    elif num_threads < 1:
        raise ValueError("num_threads must be at least 1")

    results = []
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(func, item) for item in chunk]
        for future in as_completed(futures):
            try:
                results.append(future.result())
            except Exception as e:
                print(f"An error occurred in thread: {e}")
    return results
@timeit
def run_function_on_inputs(func, i,num_processes=None, num_threads=None):
    """
    Applies the given function 'func' to each item in the list 'i' using multiprocessing and multithreading.

    Parameters:
    - i: List of inputs to process.
    - func: The function to apply to each item in 'i'.
    - num_processes: Optional; The number of processes to use. Defaults to the number of CPU cores.
    - num_threads: Optional; The number of threads to use in each process. Defaults to the number of CPU cores.

    Returns:
    - A list of results obtained by applying 'func' to each item in 'i'.
    """
    if num_processes is None:
        num_processes = multiprocessing.cpu_count()
    elif num_processes < 1:
        raise ValueError("num_processes must be at least 1")

    chunk_size = max(1, len(i) // num_processes)
    chunks = [i[x:x + chunk_size] for x in range(0, len(i), chunk_size)]
    all_results = []

    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        futures = [executor.submit(process_chunk, chunk, func, num_threads) for chunk in chunks]
        for future in as_completed(futures):
            try:
                all_results.extend(future.result())
            except Exception as e:
                print(f"An error occurred in process: {e}")
    return all_results
   
def convert_tsv_to_sqllite(tsv_file, db_name, table_name):
    with Client(n_workers = 40) as client:
        df = dd.read_csv(tsv_file, delimiter='\t', dtype=str)
        df = df.persist()
        df = df.set_index('Accession ID').persist()
        df['submitter_sample_id'] = df['Virus name'].str.split('/').str[-2].str.lower()
        wait(df)
        df.to_sql(table_name, f'sqlite:///{db_name}', if_exists='replace', index = True)
        wait(df)
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_accession_id ON metadata ("Accession ID")')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_virus_name ON metadata ("Virus name")')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_submitter_sample_id ON metadata ("submitter_sample_id")')
    conn.commit()
    conn.close()
def convert_nextclade_tsv_to_sqllite(tsv_file, db_name, table_name):
    with Client(n_workers = 40) as client:
        df = dd.read_csv(tsv_file, delimiter='\t', dtype=str)
        df['Virus name'] = df['seqName'].apply(cut_string, meta=('seqName', 'object')).str.rstrip('|')
        columns = ['Virus name'] + [col for col in df.columns if col != 'Virus name']
        df = df[columns]
        df = df.persist()
        wait(df)
        df.to_sql(table_name, f'sqlite:///{db_name}', if_exists='replace', index = False)
        wait(df)
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    cursor.execute(f'CREATE INDEX IF NOT EXISTS idx_virus_name ON {table_name} ("Virus name")')
    conn.commit()
    conn.close()
# @timeit                    
def query_metadata_db(accession):
    conn = sqlite3.connect(gisaid_metadata_sqldb)
    cursor = conn.cursor()
    query = f"""SELECT * FROM metadata
            WHERE "Virus name" = ?"""
    cursor.execute(query, (accession,))
    results = cursor.fetchone()
    conn.close()
    return results

def test_query_metadata(db_name, lineage):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    c.execute('SELECT * FROM metadata WHERE `Pango lineage` = ?', (lineage,))
    result = c.fetchall()
    conn.close()
    if result:
        return result
    else:
        return "Lineage not found"
def split_list_into_chunks(header_list, chunk_size):
    """Split the list into chunks of specified size."""
    return [header_list[i:i + chunk_size] for i in range(0, len(header_list), chunk_size)]    
    
if __name__ == "__main__":
    
    #because the gisaid sequence data and metadata are separated, its important (for easier downstream analyses) to pair them. 
    #we also must filter out any sequences that are not classified by nextclade as "good"
    #this script does both of the above tasks, creating a new directory that has one fasta file per variant (ex: B.1.160.fasta would
    #contain all sequences that are classified as "B.1.160.fasta" and passed the above QC filter).
    #The intermediate step in this script involves the creation of a sqllite database, for easy and fast searching of the sequences (which are then
    # paired with the metadata, which is also stored in a sqllite database). 
    # The metadata for each sequence is stored in the headers within the fasta files. 

    #filepaths
    #inputs
    gisaid_metadata_filepath = '../gisaidData/rawData/sars_cov_2_gisaid_metadata.tsv'
    gisaid_fasta_filepath = '../gisaidData/rawData/sars_cov_2_gisaid_sequences.fasta'
    gisaid_nextclade_qc_filepath = '../gisaidData/rawData/sars_cov_2_gisaid_sequences_nextclade_qc.tsv'
    
    #outputs
    gisaid_metadata_sqldb = '../gisaidData/sql_dbs/sars_cov_2_gisaid_metadata.db'
    gisaid_fasta_sqldb = '../gisaidData/sql_dbs/sars_cov_2_gisaid_sequences.db'
    gisaid_nextclade_qc_sqldb = '../gisaidData/sql_dbs/sars_cov_2_gisaid_sequences_nextclade_qc.db'
    variant_split_files_path = '../gisaidData/variant_split_files4'
    variant_split_files_path = variant_split_files_path.rstrip('/')
    invalid_headers_path = 'invalid_id_file_gisaid_variant_splitting4.txt'
    qc_success_path = 'gisaid_qc_success.txt'
    qc_fail_path = 'gisaid_qc_fail.txt'

    #make output dir if does not exist
    subprocess.run(f'mkdir -p {variant_split_files_path}', shell = True)
    
    # DONT FORGET TO COMMENT OUT THIS LINE IF YOU DON'T WANT TO DELETE PREV ITERATIONS
    subprocess.run(f'rm -rf {variant_split_files_path}/*', shell = True)
    
    # building sequence and metadata sqllite dbs
    convert_nextclade_tsv_to_sqllite(gisaid_nextclade_qc_filepath, gisaid_nextclade_qc_sqldb, 'nextclade_qc')
    convert_tsv_to_sqllite(gisaid_metadata_filepath,gisaid_metadata_sqldb , 'metadata')
    create_db_from_fasta_with_hash(gisaid_fasta_sqldb, 'sequences',gisaid_fasta_filepath)

   
    total = 0
    with open(invalid_headers_path,'w') as invalid_id_file, open(qc_success_path, 'w') as qc_success_file, open(qc_fail_path, 'w') as qc_fail_file:
        header_list = get_all_headers(gisaid_fasta_sqldb)
        header_list = split_list_into_chunks(header_list, 1000000)
        for index, i in enumerate(header_list):
            run_function_on_inputs(pair_gisaid_accession_with_variant, i, 10, 120)
            total += len(i)
            print(f'Finished {str(total)}')
