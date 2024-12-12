import dask.dataframe as dd
from dask.distributed import Client, wait
import sqlite3
import time
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import subprocess
from datetime import datetime
from dateutil import parser
from filelock import FileLock

def timeit(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()  # Record the start time
        result = func(*args, **kwargs)  # Execute the function
        end_time = time.time()  # Record the end time
        print(f"Function '{func.__name__}' executed in {end_time - start_time:.4f} seconds")
        return result
    return wrapper

def get_all_headers(db_name):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    
    # SQL query to select all headers
    c.execute('SELECT header FROM sequences')
    
    # Fetch all results
    headers = c.fetchall()
    
    # Close the connection
    conn.close()
    
    # Convert list of tuples to a list of strings
    header_list = [header[0] for header in headers]
    return header_list
# @timeit
def get_ncbi_sequence_hash(header):
    conn = sqlite3.connect(ncbi_seq_db)
    c = conn.cursor()
    c.execute('SELECT sequence_hash FROM sequences WHERE header = ?', (header,))
    result = c.fetchone()
    conn.close()
    if result:
        return result[0]
    else:
        return "Header not found"
def get_gisaid_sequence_hash(header):
    conn = sqlite3.connect(gisaid_seq_db)
    c = conn.cursor()
    c.execute('SELECT sequence_hash FROM sequences WHERE header = ?', (header,))
    result = c.fetchone()
    conn.close()
    if result:
        return result[0]
    else:
        return "Header not found"
# @timeit
def find_ncbi_seq_in_gisaid_db(ncbi_sequence_hash):
    conn = sqlite3.connect(gisaid_seq_db)
    c = conn.cursor()
    c.execute('SELECT header FROM sequences WHERE sequence_hash = ?', (ncbi_sequence_hash,))
    result = c.fetchall()
    if result:
        return result
    else:
        return 'seqhash not found'

def get_gisaid_sequence(header):
    conn = sqlite3.connect(gisaid_seq_db)
    c = conn.cursor()
    c.execute('SELECT sequence FROM sequences WHERE header = ?', (header,))
    result = c.fetchone()
    conn.close()
    if result:
        return result[0]
    else:
        return "Header not found"
def get_ncbi_sequence(header):
    conn = sqlite3.connect(ncbi_seq_db)
    c = conn.cursor()
    c.execute('SELECT sequence FROM sequences WHERE header = ?', (header,))
    result = c.fetchone()
    conn.close()
    if result:
        return result[0]
    else:
        return "Header not found"
    
def query_gisaid_metadata(accession):
    conn = sqlite3.connect(gisaid_meta_db)
    c = conn.cursor()
    query = f"""SELECT * FROM metadata
            WHERE "Virus name" = ?"""
    c.execute(query, (accession,))
    results = c.fetchone()
    conn.close()
    return results
def query_gisaid_metadata_submitter_id(submitter_id):
    conn = sqlite3.connect(gisaid_meta_db)
    c = conn.cursor()
    query = f"""SELECT * FROM metadata
            WHERE "submitter_sample_id" = ?"""
    c.execute(query, (submitter_id,))
    results = c.fetchone()
    return results
def query_ncbi_metadata(accession):
    conn = sqlite3.connect(ncbi_meta_db)
    c = conn.cursor()
    query = f"""SELECT * FROM metadata
            WHERE "Accession" = ?"""
    c.execute(query, (accession,))
    results = c.fetchone()
    return results
def numOfDays(date1, date2):
  #check which date is greater to avoid days output in -ve number
    if date2 > date1:   
        return (date2-date1).days
    else:
        return (date1-date2).days
def pair_ncbi_seqhash_with_gisaid_seqhash(accession):
    metadata_match = 'Undecided'
    ncbi_seq_hash = get_ncbi_sequence_hash(accession)
    ncbi_metadata = query_ncbi_metadata(accession)
    # print(ncbi_metadata)
    if ncbi_metadata[-1] != '':
        try:
            gisaid_metadata = query_gisaid_metadata_submitter_id(ncbi_metadata[-1])
            # print(gisaid_metadata)
            gisaid_hash = get_gisaid_sequence_hash(gisaid_metadata[1])
            # print(gisaid_hash)
            # sys.exit()
            if gisaid_hash == ncbi_seq_hash:
                metadata_match = 'Succeeded'
                # print('metadata match succeeded for: ', accession)
                final_result = ['' if value is None else value for value in ncbi_metadata + gisaid_metadata + ('seq_and_meta_match',)]
                return final_result
            else:
                metadata_match = 'Failed due to mismatched seqhashes despite metadata match'
                # print('Sequence hashes do not match for: ', accession)
                # sys.exit()
        except:
            metadata_match = 'Failed because ncbi metadata not found in gisaid metadata'
            # print('metadata match fail for: ', accession)
            

    gisaid_match_headers = find_ncbi_seq_in_gisaid_db(ncbi_seq_hash)
    if gisaid_match_headers == 'seqhash not found':
        if metadata_match == 'Failed due to mismatched seqhashes despite metadata match':
            final_result = ['' if value is None else value for value in ncbi_metadata + tuple(gisaid_metadata) + ('seq_mismatch_but_meta_match',)]
            return final_result
        elif metadata_match == 'Failed because ncbi metadata not found in gisaid metadata':
            final_result = ['' if value is None else value for value in ncbi_metadata + tuple('' for i in gisaid_meta_columns) + ('seq_mismatch_and_meta_mismatch',)]
            return final_result
    else:
        gisaid_match_headers = [i[0] for i in gisaid_match_headers]
        gisaid_match_metadata = [query_gisaid_metadata(i) for i in gisaid_match_headers]
        try:
            ncbi_loc = ncbi_metadata[5].split(':')[0].lower()
        except:
            ncbi_loc = ''
        # ncbi_date = datetime.strptime(ncbi_metadata[21], "%Y-%m-%d").date()
        try:
            ncbi_date = parser.parse(ncbi_metadata[21], default=datetime(1, 1, 1)).date()
        except:
            ncbi_date = parser.parse('2019-09-01', default=datetime(1, 1, 1)).date()
        # print(ncbi_date)
        gisaid_locs = []
        for i in gisaid_match_metadata:
            try:
                gisaid_locs.append(i[6].split(' / ')[1].lower())
            except:
                try:
                    gisaid_locs.append(i.strip(' /'))
                except:
                    gisaid_locs.append('')
        gisaid_dates = []
        for i in gisaid_match_metadata:
            try:
                date_string = i[5]
                gisaid_date_object = parser.parse(date_string, default=datetime(1, 1, 1)).date()
                difference = numOfDays(ncbi_date,gisaid_date_object)
                gisaid_dates.append(difference)
            except:
                gisaid_date_object = parser.parse('2019-09-01', default=datetime(1, 1, 1)).date()
                difference = numOfDays(ncbi_date,gisaid_date_object)
                gisaid_dates.append(difference)
        for index, date in enumerate(gisaid_dates):
            if date == 0 and gisaid_locs[index] == ncbi_loc:
                return ['' if value is None else value for value in ncbi_metadata + tuple(gisaid_match_metadata[index]) + ('seq_match_and_likely_meta_match',)]
            
        sorted_pairs = sorted(zip(gisaid_dates, gisaid_match_metadata, gisaid_locs))
        sorted_gisaid_dates, sorted_gisaid_metadata, gisaid_locs = zip(*sorted_pairs)
        # print(ncbi_metadata + tuple(sorted_gisaid_metadata[0]))
        
        return ['' if value is None else value for value in ncbi_metadata + tuple(sorted_gisaid_metadata[0]) + ('seq_match_and_attempted_closest_meta_date_match',)]
        
def starts_or_ends_with_invalid_char(sequence):
    valid_chars = "ATCG"
    return sequence[0] not in valid_chars or sequence[-1] not in valid_chars
@timeit            
def run_multiprocessing(function, input_list, workers):
    workers = workers
    with ProcessPoolExecutor(workers) as executor:
        results = list(executor.map(function, input_list))
    return results
@timeit            
def run_multithreading(function, input_list, workers):
    workers = workers
    with ThreadPoolExecutor(workers) as executor:
        results = list(executor.map(function, input_list))
    return results 
def get_table_headers():
    gisaid_conn = sqlite3.connect(gisaid_meta_db)
    ncbi_conn = sqlite3.connect(ncbi_meta_db)
    gisaid_cursor = gisaid_conn.cursor()
    ncbi_cursor = ncbi_conn.cursor()
    #get metadata table headers
    gisaid_cursor.execute("PRAGMA table_info(metadata)")
    gisaid_meta_columns = [i[1] + '_GISAID' for i in gisaid_cursor.fetchall()]
    ncbi_cursor.execute("PRAGMA table_info(metadata)")
    ncbi_meta_columns = [i[1] + '_NCBI' for i in ncbi_cursor.fetchall()]
    both_columns = ncbi_meta_columns + gisaid_meta_columns + ['match_type']   
    gisaid_conn.close()
    ncbi_conn.close()
    return gisaid_meta_columns, ncbi_meta_columns, both_columns
    
    
def split_list_into_chunks(header_list, chunk_size):
    """Split the list into chunks of specified size."""
    return [header_list[i:i + chunk_size] for i in range(0, len(header_list), chunk_size)]
def add_gisaid_unmatched_to_tsv(initial_combined_tsv, gisaid_metadata, out_tsv):
    with Client(n_workers = 40) as client:   
        init_combined_df = dd.read_csv(initial_combined_tsv, delimiter='\t', dtype=str)
        gisaid_meta_df = dd.read_csv(gisaid_metadata, delimiter='\t', dtype=str)
        exclude_match_types = ['seq_and_meta_match','seq_match_and_likely_meta_match', 'seq_mismatch_but_meta_match']
        both_columns = init_combined_df.columns
        print(both_columns)
        print(gisaid_meta_df.columns)
        exclude_accession_set = set(init_combined_df[init_combined_df['match_type'].isin(exclude_match_types)]['Accession ID_GISAID'].compute())
        keep_indices = [i for i in list(gisaid_meta_df['Accession ID'].compute()) if i not in exclude_accession_set]
        gisaid_meta_df = gisaid_meta_df.set_index('Accession ID')
        gisaid_meta_df = gisaid_meta_df.loc[keep_indices].persist()
        gisaid_meta_df = gisaid_meta_df.reset_index()
        gisaid_meta_df['submitter_sample_id'] = gisaid_meta_df['Virus name'].str.split('/').str[-2].str.lower()
        gisaid_meta_df.columns = gisaid_meta_df.columns + '_GISAID'
        gisaid_meta_df['match_type'] = 'seq_mismatch_and_meta_mismatch'
        gisaid_meta_df = gisaid_meta_df.persist()
        wait(gisaid_meta_df)
        for i in both_columns:
            if '_NCBI' in i:
                gisaid_meta_df[i] = ''
                gisaid_meta_df = gisaid_meta_df.persist()
        gisaid_meta_df = gisaid_meta_df[both_columns].persist()
        wait(gisaid_meta_df)
        final_df = dd.concat([gisaid_meta_df, init_combined_df]).persist()
        wait(final_df)
        # print(gisaid_meta_df.columns)
        final_df.to_csv(out_tsv, sep = '\t', index = False, single_file=True)
        wait(final_df)

def get_ncbi_headers_to_add_to_combined(ncbi_accessions_in_gisaid_tsv_out_final):
    with Client(n_workers = 40) as client:
        combined_df = dd.read_csv(ncbi_accessions_in_gisaid_tsv_out_final, delimiter='\t', dtype=str)
        categories_to_add = ['seq_match_and_attempted_closest_meta_date_match','seq_mismatch_and_meta_mismatch']
        accessions_to_add = list(combined_df[(combined_df['match_type'].isin(categories_to_add))].dropna(subset='Accession_NCBI')['Accession_NCBI'].compute())
    return accessions_to_add
        
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
                lock = FileLock(f'{combined_variant_split_files}/{lineage}.fasta.lock')
                with lock:
                    with open(f'{combined_variant_split_files}/{lineage}.fasta', 'a') as variant_file:
                        variant_file.write(f'{new_header}\n{sequence}\n')
    except:
        lock = FileLock(invalid_headers_path + '.lock')
        with lock:
            invalid_id_file.write(accession + '\n')        
        
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
def query_metadata_db(accession):
    conn = sqlite3.connect(ncbi_meta_db)
    cursor = conn.cursor()
    query = f"""SELECT * FROM metadata
            WHERE "Accession" = ?"""
    cursor.execute(query, (accession,))
    results = cursor.fetchone()
    conn.close()
    return results
def query_seq_db(accession):
    conn = sqlite3.connect(ncbi_seq_db)
    cursor = conn.cursor()
    query = f"""SELECT sequence FROM sequences
            WHERE "header" = ?"""
    cursor.execute(query, (accession,))
    results = cursor.fetchone()
    conn.close()
    return results

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


if __name__ =='__main__':
    #this data combines NCBI and GISAID metadata and sequence data for a combined dataset. 
    #input filepaths
    gisaid_meta_tsv = '../gisaidData/rawData/sars_cov_2_gisaid_metadata.tsv'
    gisaid_seq_db = '../gisaidData/sql_dbs/sequences.db'
    gisaid_meta_db = '../gisaidData/sql_dbs/metadata.db'
    ncbi_seq_db = '../ncbiData/sql_dbs/sars_cov_2_ncbi_sequences.db'
    ncbi_meta_db = '../ncbiData/sql_dbs/sars_cov_2_ncbi_metadata.db'
    ncbi_nextclade_qc_sqldb = '../ncbiData/sql_dbs/sars_cov_2_ncbi_sequences_nextclade_qc.db'
    gisaid_variant_split_files = '../gisaidData/variant_split_files4'
    
    #output filepaths
    ncbi_accessions_in_gisaid_tsv_out = '../comparison_data/check_ncbi_in_gisaid.tsv'
    ncbi_accessions_in_gisaid_tsv_out_final = '../comparison_data/check_ncbi_in_gisaid_final.tsv'
    combined_variant_split_files = '../comparison_data/variant_split_files5'
    invalid_headers_path = 'combined_invalid_ncbi_ids.txt'

    
    #metadata table headers
    gisaid_meta_columns, ncbi_meta_columns, both_columns = get_table_headers()
    #get all ncbi headers
    ncbi_headers = get_all_headers(ncbi_seq_db)
    
    #MAIN CODE BLOCK
    ncbi_headers = split_list_into_chunks(ncbi_headers, 10000)
    finished = 0
    with open(ncbi_accessions_in_gisaid_tsv_out, 'w') as outfile:
        outfile.write('\t'.join(both_columns) + '\n')
        for ncbi_header_subset in ncbi_headers:
            results = run_multithreading(pair_ncbi_seqhash_with_gisaid_seqhash, ncbi_header_subset, 120)
            for result in results:
                outfile.write('\t'.join(result) + '\n')
            finished += len(ncbi_header_subset)
            print('Finished: ', finished)
    
    # create full combined metadata
    add_gisaid_unmatched_to_tsv(ncbi_accessions_in_gisaid_tsv_out, gisaid_meta_tsv, ncbi_accessions_in_gisaid_tsv_out_final)
    
    #copy GISAID to combined variant split directory
    subprocess.run(f'cp -r {gisaid_variant_split_files} {combined_variant_split_files}', shell = True)
    
    #get list of ncbi headers to add to the combined variant split direcotry
    ncbi_headers_to_add = get_ncbi_headers_to_add_to_combined(ncbi_accessions_in_gisaid_tsv_out_final)
    total = 0
    with open(invalid_headers_path,'w') as invalid_id_file:
        header_list = split_list_into_chunks(ncbi_headers_to_add, 100000)
        for index, i in enumerate(header_list):
           
            run_function_on_inputs(pair_ncbi_accession_with_variant, i,num_processes=30, num_threads=2)
            total += len(i)
            print(f'Finished {str(total)}')
    