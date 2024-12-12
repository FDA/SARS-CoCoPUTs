import dask.dataframe as dd
from dask.distributed import Client, wait
import sqlite3
from Bio import SeqIO
import time
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import subprocess
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
def query_metadata_db(accession):
    conn = sqlite3.connect(output_meta_sqldb)
    cursor = conn.cursor()
    query = f"""SELECT * FROM metadata
            WHERE "Accession" = ?"""
    cursor.execute(query, (accession,))
    results = cursor.fetchone()
    conn.close()
    return results
def query_nextclade(db_name, accession):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    c.execute('SELECT * FROM nextclade_qc WHERE "Accession" = ?', (accession,))
    result = c.fetchone()
    conn.close()
    if result:
        return result
    else:
        return "Header not found"
def query_seq_db(accession):
    conn = sqlite3.connect(output_fasta_sqldb)
    cursor = conn.cursor()
    query = f"""SELECT sequence FROM sequences
            WHERE "header" = ?"""
    cursor.execute(query, (accession,))
    results = cursor.fetchone()
    conn.close()
    return results    
@timeit
def create_db_from_fasta_with_hash(db_name, table_name, fasta_file):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    c.execute(f'''
        CREATE TABLE IF NOT EXISTS {table_name} (
            header TEXT PRIMARY KEY,
            sequence TEXT,
            sequence_hash TEXT
        )
    ''')
    conn.commit()
    index = 0
    with open('invalid_id_file_ncbi.txt', 'w') as invalid_id_file:
        with open(fasta_file, "r") as fasta:
            for record in SeqIO.parse(fasta, "fasta"):
                try:
                    header = str(record.description.split(' ')[0])
                    
                    sequence = str(record.seq).replace('-', '')
                    
                    # Create a SHA-256 hash of the sequence
                    sequence_hash = hashlib.sha256(sequence.encode('utf-8')).hexdigest()
                    

                    try:
                        c.execute('INSERT INTO sequences (header, sequence, sequence_hash) VALUES (?, ?, ?)', 
                                  (header, sequence, sequence_hash))
                    except sqlite3.IntegrityError:  # This handles the case where the same header is encountered twice
                        print(f"Error for {header}, skipped.")
                except:
                    invalid_id_file.write('Invalid ID: ' + record.id + '\n')

                index += 1
                # if index == 10000:
                #     break
    print('Finished adding sequences to db, now making indices')
    c.execute("PRAGMA cache_size = -2000000")  # Set cache size to 2GB (in pages, assuming page size is 1024 bytes)
    c.execute("PRAGMA journal_mode = WAL")     # Use Write-Ahead Logging for better concurrency
    c.execute('CREATE INDEX IF NOT EXISTS idx_header ON sequences (header)')
    
    print('accession index finished, now making index on sequence hash column')
    c.execute('CREATE INDEX IF NOT EXISTS idx_sequence_hash ON sequences (sequence_hash)')
    conn.commit()
    conn.close() 
def convert_tsv_to_sqllite(tsv_file, db_name, table_name):
    with Client(n_workers = 40) as client:
        df = dd.read_csv(tsv_file, dtype=str, sep = '\t')
        df = df.persist()
        df = df.set_index('Accession').persist()
        df['submitter_sample_id'] = df['Isolate Lineage'].str.split('/').str[-2].str.lower()
        wait(df)
        df.to_sql(table_name, f'sqlite:///{db_name}', if_exists='replace', index = True)
        wait(df)
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    cursor.execute(f'CREATE INDEX IF NOT EXISTS idx_accession_id ON {table_name} ("Accession")')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_submitter_sample_id ON metadata ("submitter_sample_id")')
    conn.commit()
    conn.close()
def convert_nextclade_tsv_to_sqllite(tsv_file, db_name, table_name):
    with Client(n_workers = 40) as client:
        df = dd.read_csv(tsv_file, delimiter='\t', dtype=str)
        df['Accession'] = df['seqName'].str.split(' ').str[0]
        columns = ['Accession'] + [col for col in df.columns if col != 'Accession']
        df = df[columns]
        df = df.persist()
        wait(df)
        df.to_sql(table_name, f'sqlite:///{db_name}', if_exists='replace', index = False)
        wait(df)
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    cursor.execute(f'CREATE INDEX IF NOT EXISTS idx_accession_id ON {table_name} ("Accession")')
    conn.commit()
    conn.close()
def pair_ncbi_accession_with_variant(accession):
    try:
        nextclade_result = query_nextclade(ncbi_nextclade_qc_sqldb, accession)
        if nextclade_result[10] == 'good':
            query_result = query_metadata_db(accession)
            lineage,collectionDate,location,virus_name = [query_result[i] for i in [46,21,5,22]]

            sequence = query_seq_db(accession)[0]
            seq_len = len(sequence)
            n_content = bool((sequence.count('N') / seq_len) < 0.01)

            if seq_len > 29000 and n_content:
                new_header = f'>{accession}|{lineage}|{collectionDate}|{location}|{virus_name}'
                lock = FileLock(f'{variant_split_files_path}/{lineage}.fasta.lock')
                with lock:
                    with open(f'{variant_split_files_path}/{lineage}.fasta', 'a') as variant_file:
                        variant_file.write(f'{new_header}\n{sequence}\n')
    except:
        lock = FileLock(invalid_headers_path + '.lock')
        with lock:
            invalid_id_file.write(accession + '\n')
    
    
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
def split_list_into_chunks(header_list, chunk_size):
    """Split the list into chunks of specified size."""
    return [header_list[i:i + chunk_size] for i in range(0, len(header_list), chunk_size)]
if __name__ == '__main__':
    #This script functions similarly to the gisaid version, except with minor adjustments to account for differences
    #in formatting between raw GISAID data and raw NCBI data
    #inputs
    input_meta_tsv = '../ncbiData/rawData/sars_cov_2_ncbi_metadata.tsv'
    input_raw_fasta = '../ncbiData/rawData/sars_cov_2_ncbi_sequences.fa'
    ncbi_nextclade_qc_filepath = '../ncbiData/rawData/sars_cov_2_ncbi_sequences_nextclade_qc.tsv'
    #outputs
    variant_split_files_path = '../ncbiData/variant_split_files5'
    output_meta_sqldb = '../ncbiData/sql_dbs/sars_cov_2_ncbi_metadata.db'
    output_fasta_sqldb = '../ncbiData/sql_dbs/sars_cov_2_ncbi_sequences.db'
    combined_db = '../ncbiData/sql_dbs/sars_cov_2_ncbi_all_data.db'
    ncbi_nextclade_qc_sqldb = '../ncbiData/sql_dbs/sars_cov_2_ncbi_sequences_nextclade_qc.db'
    invalid_headers_path = 'invalid_id_file_ncbi_variant_splitting5.txt'
    
    #This will convert everything to sqllite databases, for fast/easy searching and matching of metadata with sequence data
    convert_tsv_to_sqllite(input_meta_tsv, output_meta_sqldb, 'metadata')
    convert_nextclade_tsv_to_sqllite(ncbi_nextclade_qc_filepath, ncbi_nextclade_qc_sqldb, 'nextclade_qc')
    create_db_from_fasta_with_hash(output_fasta_sqldb, 'sequences', input_raw_fasta)
    
    # DONT FORGET TO COMMENT OUT THIS LINE IF YOU DON'T WANT TO DELETE PREV ITERATIONS
    subprocess.run(f'rm -rf {variant_split_files_path}/*', shell = True)
    
    total = 0
    with open(invalid_headers_path,'w') as invalid_id_file:
        header_list = get_all_headers(output_fasta_sqldb)
        header_list = split_list_into_chunks(header_list, 100000)
        for index, i in enumerate(header_list):
            run_function_on_inputs(pair_ncbi_accession_with_variant, i,num_processes=10, num_threads=120)
            total += len(i)
            print(f'Finished {str(total)}')