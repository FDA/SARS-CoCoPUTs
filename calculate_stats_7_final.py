import pandas as pd
import codonbias as cb
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import  Manager
import warnings
import re
import os
import numpy as np

def split_representative_seqs_into_genes(rep_seq_file):
    #load the tsv that was generated from intial 1)sequence parsing, 2)representative sequence parsing and 3) lineage ancestry addition
    #this function will create a new tsv that has all the genes split from the initial whole genome sequences
    #of note, orf1ab will be adjusted for the frameshift identified in the literature - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7310627/
    #this is done so that the coding sequence of ORF1ab will reflect the actual translated protein
    #this must be done after whole genome secondary structure prediction, since preserving actual natural sequence is important for this step
    wholeseqdf = pd.read_csv(rep_seq_file, sep = '\t')
    
    #adjust whole genome sequence for orf1ab frameshift
    #all gene indices within this code are from database/literature sources
    reference_sequence = wholeseqdf[wholeseqdf['terminal_lineage'] == 'referenceB']['sequence'].values[0]
    orf1ab_cds = reference_sequence.replace('-','')[265:21555]
    orf1ab_cds_preframeshift = orf1ab_cds[0:13203]
    #the regex search process below is basically just to get the end index of pre-frameshift orf1ab within the MSA
    regex_string = ''
    for i in orf1ab_cds_preframeshift:
        regex_string += i
        regex_string += '-*'
    #this just cuts off the last '-*'
    regex_string = regex_string[0:-2]
    #now search for the string, get the indices of orf1ab pre-frameshift (only the end index is needed below)
    y = re.search(regex_string, reference_sequence)
    #the goal is to repeat one nucleotide (in this case a "C" to mimic the coding sequence that matches the actual translated sequence
    #this loop will increase the subtract index until the first value of the post-frameshift sequence is not a '-' 
    #since we want to continue moving backwards until we are repeating an actual nucleotide (in this case a 'C'), not a '-'
    #realistically this is just to account for an edge case if the MSA changes in the future and a GAP is added at this location
    subtract_index = 1
    continue_subtracting = True
    while continue_subtracting:
        post_frameshift_sequence_msa = reference_sequence[y.end()-subtract_index:]
        if post_frameshift_sequence_msa[0] != '-':
            break
        else:
            subtract_index +=1
    #this line adjusts the entire MSA for this adjustment
    wholeseqdf['sequence'] = wholeseqdf['sequence'].str[:y.end()] + wholeseqdf['sequence'].str[y.end() -subtract_index:]
    
    #now that the frameshift is accounted for, we can begin splitting sequences into individual genes
    #redefine reference sequence after adjustment
    reference_sequence = wholeseqdf[wholeseqdf['terminal_lineage'] == 'referenceB']['sequence'].values[0]
    #all of these indices are from gisaid (https://gisaid.org/wiv04/), but adjusted for python indexing which starts at 0
    #all except the first start index have a +1 added due to the frameshift adjustment
    #the below block will get just wild type coding sequences
    cds_seq_names =['orf1ab_cds',
                    'spike_cds',
                    'orf3a_cds',
                    'orf3b_cds',
                    'e_cds',
                    'm_cds',
                    'orf6_cds',
                    'orf7a_cds',
                    'orf7b_cds',
                    'orf8_cds',
                    'n_cds',
                    'orf10_cds']
    start_indices = [265,21563,25393,25765,26245,26523,27202,27394,27756,27894,28274,29558]
    end_indices = [21556,25385, 26221,26221,26473,27192,27388,27760,27888,28260,29534,29675]
    cds_seqs = [reference_sequence.replace('-','')[i:end_indices[index]] for index, i in enumerate(start_indices)]
    #now we will get indices for the whole MSA for each coding sequence
    start_indices_msa = []
    end_indices_msa = []
    msa_cds_seqs = []
    for cds_seq in cds_seqs:
        regex_string = ''
        for i in cds_seq:
            regex_string += i
            regex_string += '-*'
        regex_string = regex_string[0:-2]
        x = re.search(regex_string, reference_sequence)
        start_indices_msa.append(x.start())
        end_indices_msa.append(x.end())
        msa_cds_seqs.append(reference_sequence[x.start():x.end()])
    
    for index, i in enumerate(cds_seq_names):
        wholeseqdf[i] = wholeseqdf['sequence'].str[start_indices_msa[index]:end_indices_msa[index]]
    
    #now repeat for noncoding sequences
    noncoding_starts = [0] + end_indices_msa
    noncoding_ends = start_indices_msa + [len(reference_sequence)]
    noncoding_seq_names =['pre_orf1ab_sequence',
                          'pre_spike_sequence',
                          'pre_orf3a_sequence',
                          'pre_orf3b_sequence',
                          'pre_e_sequence',
                          'pre_m_sequence',
                          'pre_orf6_sequence',
                          'pre_orf7a_sequence',
                          'pre_orf7b_sequence',
                          'pre_orf8_sequence',
                          'pre_n_sequence',
                          'pre_orf10_sequence',
                          'post_orf10_sequence']
    #this cleaning is necessary because some coding sequences overlap with each other and so it is necessary to remove their indices 
    #since there are no noncoding sequences between them
    clean_noncoding_seq_names = [i for index, i in enumerate(noncoding_seq_names) if noncoding_starts[index] < noncoding_ends[index]]
    clean_noncoding_starts = [i for index, i in enumerate(noncoding_starts) if noncoding_starts[index] < noncoding_ends[index]]
    clean_noncoding_ends = [i for index, i in enumerate(noncoding_ends) if noncoding_starts[index] < noncoding_ends[index]]
    for index, i in enumerate(clean_noncoding_seq_names):
        wholeseqdf[i] = wholeseqdf['sequence'].str[clean_noncoding_starts[index]:clean_noncoding_ends[index]]
    return wholeseqdf, cds_seq_names, noncoding_seq_names

def get_mutations(query, reference):
    #this function compares a query to a reference (both aligned) and determines the synonymous, nonsynonymous, insertion, and deletion mutations
    #create a dictionary of amino acids to codons
    amino_acid_to_codons = {
    "A": ["GCT", "GCC", "GCA", "GCG"],
    "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "N": ["AAT", "AAC"],
    "D": ["GAT", "GAC"],
    "C": ["TGT", "TGC"],
    "Q": ["CAA", "CAG"],
    "E": ["GAA", "GAG"],
    "G": ["GGT", "GGC", "GGA", "GGG"],
    "H": ["CAT", "CAC"],
    "I": ["ATT", "ATC", "ATA"],
    "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "K": ["AAA", "AAG"],
    "M": ["ATG"],
    "F": ["TTT", "TTC"],
    "P": ["CCT", "CCC", "CCA", "CCG"],
    "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "T": ["ACT", "ACC", "ACA", "ACG"],
    "W": ["TGG"],
    "Y": ["TAT", "TAC"],
    "V": ["GTT", "GTC", "GTA", "GTG"],
    "*": ["TAA", "TAG", "TGA"]}
    codon_to_amino_acid = {
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "AAT": "N", "AAC": "N",
    "GAT": "D", "GAC": "D",
    "TGT": "C", "TGC": "C",
    "CAA": "Q", "CAG": "Q",
    "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATT": "I", "ATC": "I", "ATA": "I",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "AAA": "K", "AAG": "K",
    "ATG": "M",
    "TTT": "F", "TTC": "F",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TAA": "*", "TAG": "*", "TGA": "*"}

    #remove where both are dashes
    clean_reference = ''.join([i for index, i in enumerate(reference) if i != '-' or query[index] != '-'])
    clean_query = ''.join([i for index, i in enumerate(query) if i != '-' or reference[index] != '-'])
    #get number of codons in reference
    codons = len(clean_reference.replace('-',''))/3
    #check if codonis whole number
    if codons % 1 != 0:
        raise ValueError("Codons in reference are not whole numbers")
    codons_ref_msa =[]
    ref_index = []
    codons_query_msa = []
    codon_ref = ''
    codon_query = ''
    codon_ref_index = ''
    codon_count = 0
    for index, i in enumerate(clean_reference):
        if i != '-':
            codon_ref += i
            codon_query += clean_query[index]
            codon_ref_index += str(index +1)+' '
            codon_count += 1
        else:
            codon_ref += i
            codon_query += clean_query[index]
            codon_ref_index += '-' + ' '
        if codon_count == 3:
            codons_ref_msa.append(codon_ref)
            codons_query_msa.append(codon_query)
            ref_index.append(codon_ref_index.strip(' '))
            codon_ref = ''
            codon_query = ''
            codon_ref_index = ''
            codon_count = 0
    mutations = []
    syn_removed_query_seq = ''
    for index, ref_codon in enumerate(codons_ref_msa):
        query_codon = codons_query_msa[index]
        if ref_codon != query_codon and all(char in "ATCG-" for char in query_codon):
            
            position = [i for i in ref_index[index].split(' ') if i != '-'][0]
            
            if ref_codon.find('-') > -1:
                codon_mutation = f'{position}_{ref_codon}>{query_codon}_INSERTION'
                syn_removed_query_seq += query_codon                
            elif query_codon.find('-') > -1:
                codon_mutation = f'{position}_{ref_codon}>{query_codon}_DELETION'
                syn_removed_query_seq += query_codon
            elif codon_to_amino_acid[ref_codon] ==  codon_to_amino_acid[query_codon]:
                codon_mutation = f'{position}_{ref_codon}>{query_codon}_SYNONYMOUS'
                syn_removed_query_seq += ref_codon
            elif codon_to_amino_acid[ref_codon] !=  codon_to_amino_acid[query_codon]:
                codon_mutation = f'{position}_{ref_codon}>{query_codon}_NONSYNONYMOUS'
                syn_removed_query_seq += query_codon
            mutations.append(codon_mutation)
        else:
            syn_removed_query_seq += query_codon
    syn_removed_query_seq = syn_removed_query_seq.replace('-','')
    final_syn_removed_query_seq = ''
    count = 0
    for index, nt in enumerate(query):
        if nt == '-':
            final_syn_removed_query_seq += nt
        else:
            final_syn_removed_query_seq += syn_removed_query_seq[count]
            count +=1
    if len(final_syn_removed_query_seq) != len(query):
        raise Exception('error')
    return [mutations, final_syn_removed_query_seq, query]
def create_mutation_df_and_syn_removed_df(split_gene_df, cds_seq_names):
    # create syn removed cds, basically just applies the get_mutations function
    mutations_df = split_gene_df[split_gene_df.columns[:8]]
    syn_removed_cds_df = split_gene_df[split_gene_df.columns[:8]]
    for index, i in enumerate(cds_seq_names[::-1]):
        reference_sequence = split_gene_df[split_gene_df['terminal_lineage'] == 'referenceB'][i].values[0]
        tmp_df = pd.DataFrame(split_gene_df[i].apply(get_mutations, reference=reference_sequence).to_list(), columns=[f'mutations_{i}', f'syn_removed_{i}', i])
        mutations_df = pd.concat([mutations_df, tmp_df[f'mutations_{i}']], axis = 1)
        syn_removed_cds_df = pd.concat([syn_removed_cds_df, tmp_df[f'syn_removed_{i}']], axis = 1)
    return mutations_df, syn_removed_cds_df
def finalize_syn_removed_df(syn_removed_cds_df,wt_split_gene_df):
    orf7b_start_index = 4
    continue_adding = True
    while continue_adding:
        orf7b_cds_reference = wt_split_gene_df[wt_split_gene_df['terminal_lineage'] == 'referenceB']['orf7b_cds'].values[0]
        overlap = orf7b_cds_reference[:orf7b_start_index]
        if '-' not in overlap:
            break
        else:
            orf7b_start_index +=1
    #we have to undo the added nucleotide to orf1ab since the main point of doing this synonymous mutation removal is to 
    #observe its impact on secondary structure, and therefore we want the sequences to be otherwise as close to WT as possible
    reference_sequence = wt_split_gene_df[wt_split_gene_df['terminal_lineage'] == 'referenceB']['sequence'].values[0]
    orf1ab_cds = reference_sequence.replace('-','')[265:21555]
    orf1ab_cds_preframeshift = orf1ab_cds[0:13203]
    regex_string = ''
    for i in orf1ab_cds_preframeshift:
        regex_string += i
        regex_string += '-*'
    #this just cuts off the last '-*'
    regex_string = regex_string[0:-2]
    #now search for the string, get the indices of orf1ab pre-frameshift (only the end index is needed below)
    wt_syn_removed_orf1ab_cds = syn_removed_cds_df[syn_removed_cds_df['terminal_lineage'] == 'referenceB']['syn_removed_orf1ab_cds'].values[0]
    y = re.search(regex_string, wt_syn_removed_orf1ab_cds)
    orf1ab_add_index = 1
    continue_adding = True
    while continue_adding:
        post_frameshift_sequence_msa = wt_syn_removed_orf1ab_cds[y.end()+orf1ab_add_index:]
        if post_frameshift_sequence_msa[0] != '-':
            break
        else:
            orf1ab_add_index +=1
    #this line adjusts the entire MSA for this adjustment
    syn_removed_cds_df['syn_removed_orf1ab_cds'] = syn_removed_cds_df['syn_removed_orf1ab_cds'].str[:y.end()] + syn_removed_cds_df['syn_removed_orf1ab_cds'].str[y.end() +orf1ab_add_index:]
    
    
    
    final_syn_removed_df = syn_removed_df[syn_removed_df.columns[:8]]
    
    insert_column = wt_split_gene_df['pre_orf1ab_sequence'] + syn_removed_cds_df['syn_removed_orf1ab_cds'] + wt_split_gene_df['pre_spike_sequence'] + syn_removed_cds_df['syn_removed_spike_cds'] + wt_split_gene_df['pre_orf3a_sequence'] + syn_removed_cds_df['syn_removed_orf3a_cds'] + wt_split_gene_df['pre_e_sequence'] + syn_removed_cds_df['syn_removed_e_cds'] + wt_split_gene_df['pre_m_sequence'] + syn_removed_cds_df['syn_removed_m_cds'] + wt_split_gene_df['pre_orf6_sequence'] + syn_removed_cds_df['syn_removed_orf6_cds'] + wt_split_gene_df['pre_orf7a_sequence'] +syn_removed_cds_df['syn_removed_orf7a_cds'] +syn_removed_cds_df['syn_removed_orf7b_cds'].str[orf7b_start_index:] + wt_split_gene_df['pre_orf8_sequence'] + syn_removed_cds_df['syn_removed_orf8_cds'] + wt_split_gene_df['pre_n_sequence'] + syn_removed_cds_df['syn_removed_n_cds'] + wt_split_gene_df['pre_orf10_sequence'] + syn_removed_cds_df['syn_removed_orf10_cds'] + wt_split_gene_df['post_orf10_sequence']
    
    final_syn_removed_df.drop('sequence',inplace = True, axis = 1)
    final_syn_removed_df.insert(6, 'sequence',insert_column)
    return final_syn_removed_df

def run_multiprocessing(function, input_list, workers):
    workers = workers
    with ProcessPoolExecutor(workers) as executor:
        results = executor.map(function, input_list)
        for i in results:
            pass
def get_substrings_with_indices(sequence):
    # Regex to match groups of "-" or groups of non "-" characters
    pattern = re.compile(r'-+|[^-]+')
    
    substrings_with_indices = []
    
    # Use finditer to get match objects, including start and end indices
    for match in pattern.finditer(sequence):
        start = match.start()
        end = match.end()
        substrings_with_indices.append([match.group(), start, end])
        
    return substrings_with_indices

def split_into_codons_and_get_indices(sequence):
    codon_list = []
    codon = ''
    for index, i in enumerate(sequence):
        codon += i


def adjust_seqs (reference_seq, query_seq):
    cleaned_ref_tmp = ''
    cleaned_query_tmp = ''
    for index, i in enumerate(reference_seq):
        if i != '-' or query_seq[index] != '-':
            cleaned_ref_tmp += i
            if query_seq[index] in '-ATCG':

                cleaned_query_tmp += query_seq[index]
            else:
                cleaned_query_tmp += i
    #do it again in case replacing N's created gaps
    cleaned_ref = ''
    cleaned_query = ''
    for index, i in enumerate(cleaned_ref_tmp):
        if i != '-' or cleaned_query_tmp[index] != '-':
            cleaned_ref += i
            cleaned_query += cleaned_query_tmp[index]
    print(cleaned_ref)
    print(cleaned_query)

    
    





def get_orf(query_seq, column_name):
    query_seq = query_seq.strip('\n')

    reference_seq = sars_referenceB[column_name].strip('\n')
    query_seq = ''.join([i if i in '-ATCG' else reference_seq[index] for index, i in enumerate(query_seq)])
    query_seq = query_seq.replace('-','')
    possible_orfs = []
    for reading_frame_start in range(3):
        seq = cb.utils.translate(query_seq[reading_frame_start:], return_str = True, genetic_code = 1).split('*')[0] + '*'
        possible_orfs.append(len(seq))
    start_index = np.argmax(possible_orfs)
    end_index = start_index + (possible_orfs[start_index] * 3)
    return query_seq[start_index:end_index]

def calculate_codon_usage(sequence):
    seq = sequence.replace('-','')
    codon_stats = cb.stats.CodonCounter(seqs = seq, ignore_stop = False)
    codon_usage = dict(codon_stats.get_codon_table())
    gc_content = gc_codon(seq)
    seqlength = {'adjusted_cds_length':len(seq)}
    gc_content.update(seqlength)
    gc_content.update(codon_usage)
    return gc_content
def calculate_bicodon_usage(sequence):
    seq = sequence.replace('-','')
    bicodon_stats = cb.stats.CodonCounter(seqs = seq, k_mer=2, ignore_stop = False)
    bicodon_stats = dict(bicodon_stats.get_codon_table())
    seqlength = {'adjusted_cds_length':len(seq)}
    seqlength.update(bicodon_stats)
    return seqlength
def calculate_dinucleotides(sequence):
    seq = sequence.replace('-','')
    nt_stats = cb.stats.BaseCounter(seqs=seq, k_mer=2, step=1)
    dinucleotide_stats = dict(nt_stats.get_table())
    gc_content = gc_dinucleotide(seq)
    seqlength = {'adjusted_cds_length':len(seq)}
    gc_content.update(seqlength)
    gc_content.update(dinucleotide_stats)
    return gc_content
def calculate_junction_dinucleotides(sequence):
    seq = sequence.replace('-','')
    nt_stats = cb.stats.BaseCounter(seqs=seq, k_mer=2, step=3, frame = 3)
    junction_dnt = dict(nt_stats.get_table())
    gc_content = gc_junction_dinucleotide(seq)
    seqlength = {'adjusted_cds_length':len(seq)}
    gc_content.update(seqlength)
    gc_content.update(junction_dnt)
    return gc_content
def calc_enc(sequence, bg_correction, k_mer):
    seq = sequence.replace('-','')
    enc = cb.scores.EffectiveNumberOfCodons(k_mer=k_mer, bg_correction=bg_correction)
    gene_scores = enc.get_score(seq)
    return gene_scores
def gc_codon(sequence):
    seq = sequence.replace('-','')
    nt_stats = cb.stats.BaseCounter(seqs=seq, k_mer=1, step=1, frame = 1)
    wholeseq_stats = nt_stats.get_table()
    gc = (wholeseq_stats['G'] + wholeseq_stats['C']) / wholeseq_stats.sum()
    nt_stats = cb.stats.BaseCounter(seqs=seq, k_mer=1, step=3, frame = 1)
    wholeseq_stats = nt_stats.get_table()
    gc1 = (wholeseq_stats['G'] + wholeseq_stats['C']) / wholeseq_stats.sum()
    nt_stats = cb.stats.BaseCounter(seqs=seq, k_mer=1, step=3, frame = 2)
    wholeseq_stats = nt_stats.get_table()
    gc2 = (wholeseq_stats['G'] + wholeseq_stats['C']) / wholeseq_stats.sum()
    nt_stats = cb.stats.BaseCounter(seqs=seq, k_mer=1, step=3, frame = 3)
    wholeseq_stats = nt_stats.get_table()
    gc3 = (wholeseq_stats['G'] + wholeseq_stats['C']) / wholeseq_stats.sum()
    return {'GC1%':gc1,'GC2%':gc2,'GC3%':gc3, 'GC%':gc}
def gc_dinucleotide(sequence):
    seq = sequence.replace('-','')
    nt_stats_gc1 = cb.stats.BaseCounter(seqs=seq[:-1], k_mer=1, step=1, frame = 1)
    wholeseq_stats_gc1 = nt_stats_gc1.get_table()
    gc1 = (wholeseq_stats_gc1['G'] + wholeseq_stats_gc1['C']) / wholeseq_stats_gc1.sum()
    nt_stats_gc2 = cb.stats.BaseCounter(seqs=seq[1:], k_mer=1, step=1, frame = 1)
    wholeseq_stats_gc2 = nt_stats_gc2.get_table()
    gc2 = (wholeseq_stats_gc2['G'] + wholeseq_stats_gc2['C']) / wholeseq_stats_gc2.sum()
    gc = (gc1 + gc2) / 2
    return {'GC1%':gc1,'GC2%':gc2, 'GC%':gc}
def gc_junction_dinucleotide(sequence):
    seq = sequence.replace('-','')
    final_index = ((len(seq) // 3) * 3) - 1
    seq = seq[:final_index]
    nt_stats_gc1 = cb.stats.BaseCounter(seqs=seq, k_mer=1, step=3, frame = 3)
    wholeseq_stats_gc1 = nt_stats_gc1.get_table()
    gc1 = (wholeseq_stats_gc1['G'] + wholeseq_stats_gc1['C']) / wholeseq_stats_gc1.sum()
    nt_stats_gc2 = cb.stats.BaseCounter(seqs=seq, k_mer=1, step=3, frame = 4)
    wholeseq_stats_gc2 = nt_stats_gc2.get_table()
    gc2 = (wholeseq_stats_gc2['G'] + wholeseq_stats_gc2['C']) / wholeseq_stats_gc2.sum()
    wholeseq_stats = wholeseq_stats_gc1 + wholeseq_stats_gc2
    gc = (wholeseq_stats['G'] + wholeseq_stats['C']) / wholeseq_stats.sum()
    return {'GC1%':gc1,'GC2%':gc2, 'GC%':gc}
def read_fasta(file_path):
    sequences = []
    with open(file_path, 'r') as file:
        sequence = ''
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if sequence:
                    sequences.append(sequence)
                    sequence = ''
            else:
                sequence += line
        if sequence:
            sequences.append(sequence)
    return sequences

def calc_cai(sequence):
    return cai_object.get_score(sequence)
def calc_cpai(sequence):
    return cpai_object.get_score(sequence)

def run_calcs_on_columns(column_name):
    i = column_name
    tmp_df = split_gene_df[['parent_lineages','terminal_lineage', i]]
    tmp_df[i] = tmp_df[i].apply(get_orf, args = (i,))
    print(f'Starting codon Usage calc for {i}')
    codon_usage_df = tmp_df.copy()
    codon_usage_df['gene'] = i.replace('_cds','')
    codon_usage_df[i] = codon_usage_df[i].apply(calculate_codon_usage)
    codon_usage_df.rename(columns={i:'codon_usage'}, inplace = True)
    codon_usage_df = pd.concat([codon_usage_df.drop(['codon_usage'], axis=1), codon_usage_df['codon_usage'].apply(pd.Series)], axis=1)
    final_codon_usage_df[i] = codon_usage_df

    #bicodon usage
    print(f'Starting bicodon Usage calc for {i}')
    bicodon_usage_df = tmp_df.copy()
    bicodon_usage_df['gene'] = i.replace('_cds','')
    bicodon_usage_df[i] = bicodon_usage_df[i].apply(calculate_bicodon_usage)
    bicodon_usage_df.rename(columns={i:'bicodon_usage'}, inplace = True)
    bicodon_usage_df = pd.concat([bicodon_usage_df.drop(['bicodon_usage'], axis=1), bicodon_usage_df['bicodon_usage'].apply(pd.Series)], axis=1)
    final_bicodon_usage_df[i] = bicodon_usage_df
    
    #dinucleotide usage
    print(f'Starting dinucleotide Usage calc for {i}')
    dinucleotide_usage_df = tmp_df.copy()
    dinucleotide_usage_df['gene'] = i.replace('_cds','')
    dinucleotide_usage_df[i] = dinucleotide_usage_df[i].apply(calculate_dinucleotides)
    dinucleotide_usage_df.rename(columns={i:'dinucleotide_usage'}, inplace = True)
    dinucleotide_usage_df = pd.concat([dinucleotide_usage_df.drop(['dinucleotide_usage'], axis=1), dinucleotide_usage_df['dinucleotide_usage'].apply(pd.Series)], axis=1)
    final_dinucleotide_df[i] = dinucleotide_usage_df
    
    #junction dinucleotide usage
    print(f'Starting junction dinucleotide Usage calc for {i}')
    junction_dinucleotide_usage_df = tmp_df.copy()
    junction_dinucleotide_usage_df['gene'] = i.replace('_cds','')
    junction_dinucleotide_usage_df[i] = junction_dinucleotide_usage_df[i].apply(calculate_junction_dinucleotides)
    junction_dinucleotide_usage_df.rename(columns={i:'junction_dinucleotide_usage'}, inplace = True)
    junction_dinucleotide_usage_df = pd.concat([junction_dinucleotide_usage_df.drop(['junction_dinucleotide_usage'], axis=1), junction_dinucleotide_usage_df['junction_dinucleotide_usage'].apply(pd.Series)], axis=1)
    final_junction_dinucleotide_df[i] = junction_dinucleotide_usage_df
    
    #enc and CAI
    enc_df = tmp_df.copy()
    enc_df['gene'] = i.replace('_cds','')
    enc_df['adjusted_cds_length'] = enc_df[i].apply(len)
    
    print(f'Starting uncorrected enc calc for {i}')
    enc_df['ENC'] = enc_df[i].apply(calc_enc, args=[False, 1])
    
    print(f'Starting corrected enc calc for {i}')   
    enc_df['ENC_GC_corrected'] = enc_df[i].apply(calc_enc, args=[True, 1])
    
    print(f'Starting codon pair uncorrected enc calc for {i}')
    enc_df['ENC_CP'] = enc_df[i].apply(calc_enc, args=[False, 2])
    
    print(f'Starting codon pair corrected enc calc for {i}')    
    enc_df['ENC_CP_GC_corrected'] = enc_df[i].apply(calc_enc, args=[True, 2])
    
    print(f'Starting CAI calculation for {i}')
    enc_df['CAI'] = enc_df[i].apply(calc_cai)
    
    print(f'Starting CPAI calculation for {i}')
    enc_df['CPAI'] = enc_df[i].apply(calc_cpai)
    final_enc_df[i] = enc_df.drop([i], axis = 1)
    
    return f'Finished {column_name}'
def add_whole_cds_calcs_codons(final_codon_usage_df, cocoputs_codon_order):
    #add whole cds calc
    tmp_final_codon_usage_df = final_codon_usage_df.copy()
    gc_columns = ['GC1%','GC2%','GC3%','GC%']
    for i in gc_columns:
        tmp_final_codon_usage_df[i] *= tmp_final_codon_usage_df['adjusted_cds_length']
    tmp_final_codon_usage_df = tmp_final_codon_usage_df.groupby(['parent_lineages','terminal_lineage']).sum()
    for i in gc_columns:
        tmp_final_codon_usage_df[i] /= tmp_final_codon_usage_df['adjusted_cds_length']
    tmp_final_codon_usage_df = tmp_final_codon_usage_df.reset_index()
    tmp_final_codon_usage_df['gene'] = 'whole_cds'
    final_codon_usage_df = pd.concat([final_codon_usage_df, tmp_final_codon_usage_df])
    #change order of codon columns
    other_cols = list(final_codon_usage_df.columns[:8])
    other_cols.extend(cocoputs_codon_order)
    final_codon_usage_df = final_codon_usage_df[other_cols] 
    return final_codon_usage_df
def add_whole_cds_calcs_bicodons(final_bicodon_usage_df, cocoputs_bicodons_order):
    #add whole cds calc
    tmp_final_bicodon_usage_df = final_bicodon_usage_df.copy()
    tmp_final_bicodon_usage_df = tmp_final_bicodon_usage_df.groupby(['parent_lineages','terminal_lineage']).sum()
    tmp_final_bicodon_usage_df = tmp_final_bicodon_usage_df.reset_index()
    tmp_final_bicodon_usage_df['gene'] = 'whole_cds'
    final_bicodon_usage_df = pd.concat([final_bicodon_usage_df, tmp_final_bicodon_usage_df])
    other_cols = list(final_bicodon_usage_df.columns[:4])
    other_cols.extend(cocoputs_bicodons_order)
    final_bicodon_usage_df = final_bicodon_usage_df[other_cols]
    return final_bicodon_usage_df
def add_whole_cds_calcs_dinucleotides(final_dinucleotide_df):
    tmp_final_dinucleotide_df = final_dinucleotide_df.copy()
    gc_columns = ['GC1%','GC2%','GC%']
    for i in gc_columns:
        tmp_final_dinucleotide_df[i] *= tmp_final_dinucleotide_df['adjusted_cds_length']
    tmp_final_dinucleotide_df = tmp_final_dinucleotide_df.groupby(['parent_lineages','terminal_lineage']).sum()
    for i in gc_columns:
        tmp_final_dinucleotide_df[i] /= tmp_final_dinucleotide_df['adjusted_cds_length']
    tmp_final_dinucleotide_df = tmp_final_dinucleotide_df.reset_index()
    tmp_final_dinucleotide_df['gene'] = 'whole_cds'
    final_dinucleotide_df = pd.concat([final_dinucleotide_df, tmp_final_dinucleotide_df])
    cocoputs_dinucleotide_order = ["TT", "TC", "TA", "TG", "CT", "CC", "CA", "CG", "AT", "AC", "AA", "AG", "GT", "GC", "GA", "GG"]
    other_cols = list(final_dinucleotide_df.columns[:7])
    other_cols.extend(cocoputs_dinucleotide_order)
    final_dinucleotide_df = final_dinucleotide_df[other_cols]
    return final_dinucleotide_df
def add_whole_cds_calcs_enc_cai(final_enc_df):
    tmp_final_enc_df = final_enc_df.copy()
    gc_columns = ['ENC','ENC_GC_corrected','ENC_CP','ENC_CP_GC_corrected','CAI','CPAI']
    for i in gc_columns:
        tmp_final_enc_df[i] *= tmp_final_enc_df['adjusted_cds_length']
    tmp_final_enc_df = tmp_final_enc_df.groupby(['parent_lineages','terminal_lineage']).sum()
    for i in gc_columns:
        tmp_final_enc_df[i] /= tmp_final_enc_df['adjusted_cds_length']
    tmp_final_enc_df = tmp_final_enc_df.reset_index()
    tmp_final_enc_df['gene'] = 'whole_cds'
    final_enc_df = pd.concat([final_enc_df, tmp_final_enc_df])
    return final_enc_df

def create_outfile_name(input_tsv, prefix, replace_string):
    dir_name = os.path.dirname(input_tsv)
    base_name = os.path.basename(input_tsv).replace(replace_string,"")
    out_name = f"{prefix}_{base_name}"
    aligned_filename = os.path.join(dir_name, out_name)
    print(aligned_filename)
    return aligned_filename
if __name__ == '__main__':
    #this script performs many of the sequence composition calculations on our representative, consensus and ancestral sequences,
    #and also generates our synonymous-removed datasets. 
    #there are many warnings, so this helps to simplify outputs
    warnings.filterwarnings("ignore")    
    #file paths
    #input
    input_tsvs =[
        '../gisaidData/final_files_USA_only/all_lineages_added_aligned_ancestral_seqs.tsv',
        '../gisaidData/final_files_USA_only/all_lineages_added_aligned_consensus_seqs.tsv',
        '../gisaidData/final_files_USA_only/all_lineages_added_aligned_representative_seqs.tsv',
        '../ncbiData/final_files_USA_only/all_lineages_added_aligned_ancestral_seqs.tsv',
        '../ncbiData/final_files_USA_only/all_lineages_added_aligned_consensus_seqs.tsv',
        '../ncbiData/final_files_USA_only/all_lineages_added_aligned_representative_seqs.tsv',
        '../comparison_data/final_files_USA_only/all_lineages_added_aligned_ancestral_seqs.tsv',
        '../comparison_data/final_files_USA_only/all_lineages_added_aligned_consensus_seqs.tsv',
        '../comparison_data/final_files_USA_only/all_lineages_added_aligned_representative_seqs.tsv'  
    ]
    codon_order = '../other_data/codon_order.txt'
    human_coding_sequences = '../other_data/coding_sequences_mane.fasta'

    for input_tsv in input_tsvs:
        #output
        split_gene_file = create_outfile_name(input_tsv, 'split_genes', 'all_lineages_added_aligned_')
        mutation_out_file = create_outfile_name(input_tsv, 'mutation_data', 'all_lineages_added_aligned_')
        syn_out_file = create_outfile_name(input_tsv, 'syn_removed', 'all_lineages_added_aligned_')
        codon_out_file = create_outfile_name(input_tsv, 'codon_usage', 'all_lineages_added_aligned_')
        bicodon_out_file = create_outfile_name(input_tsv, 'bicodon_usage', 'all_lineages_added_aligned_')
        dinucleotide_out_file = create_outfile_name(input_tsv, 'dinucleotide_usage', 'all_lineages_added_aligned_')
        jn_dinucleotide_out_file = create_outfile_name(input_tsv, 'jn_dinucleotide_usage', 'all_lineages_added_aligned_')
        enc_out_file = create_outfile_name(input_tsv, 'enc_cai', 'all_lineages_added_aligned_')

        #split the whole sequences into individual genes
        split_gene_df, cds_seq_names, noncoding_seq_names = split_representative_seqs_into_genes(input_tsv)
        split_gene_df.to_csv(split_gene_file, sep = '\t', index = False)
        
        #get mutations and synonymous mutation redacted sequences
        mutation_df, syn_removed_df = create_mutation_df_and_syn_removed_df(split_gene_df, cds_seq_names)
        #save the mutation file to a tsv
        mutation_df.to_csv(mutation_out_file, sep = '\t', index = False)
        #finalize the synonymous redacted csv
        final_syn_removed_dataframe = finalize_syn_removed_df(syn_removed_df, split_gene_df)
        final_syn_removed_dataframe.to_csv(syn_out_file,sep = '\t', index = False)

        #run codon,codon pair, ENC, etc calculations
        #before running CAI calculations below we must initialize cai/cpai objects
        ref_seq = read_fasta(human_coding_sequences)
        cai_object = cb.scores.CodonAdaptationIndex(ref_seq, k_mer=1, genetic_code=1, ignore_stop=False, pseudocount=1)
        cpai_object = cb.scores.CodonAdaptationIndex(ref_seq, k_mer=2, genetic_code=1, ignore_stop=False, pseudocount=1)

        #Get the correct order of the codon columns - this is only necessary for compatibility with CoCoPUTs
        with open(codon_order, 'r') as codon_order_file:
            string_file = codon_order_file.read().replace('"','')
            codon_order_list = string_file.split('\n')[0].split(',')
            bicodon_order_list = string_file.split('\n')[1].split(',')
        
        sars_referenceB = split_gene_df[split_gene_df['terminal_lineage'] == 'referenceB'][cds_seq_names]
        sars_referenceB = dict(zip(sars_referenceB.columns, sars_referenceB.values[0]))
        
        manager = Manager()
        final_codon_usage_df = manager.dict()
        final_bicodon_usage_df = manager.dict()
        final_dinucleotide_df = manager.dict()
        final_junction_dinucleotide_df = manager.dict()
        final_enc_df = manager.dict()
        
        run_multiprocessing(run_calcs_on_columns, cds_seq_names, 20)
        
        final_codon_usage_df = add_whole_cds_calcs_codons(pd.concat(list(final_codon_usage_df.values())), codon_order_list)
        gc_columns_codon = ['GC1%', 'GC2%', 'GC3%', 'GC%']
        final_codon_usage_df[gc_columns_codon] *= 100
        final_bicodon_usage_df = add_whole_cds_calcs_bicodons(pd.concat(list(final_bicodon_usage_df.values())), bicodon_order_list)
        final_dinucleotide_df = add_whole_cds_calcs_dinucleotides(pd.concat(list(final_dinucleotide_df.values())))
        gc_columns_dn = ['GC1%', 'GC2%', 'GC%']
        final_dinucleotide_df[gc_columns_dn] *= 100
        final_junction_dinucleotide_df = add_whole_cds_calcs_dinucleotides(pd.concat(list(final_junction_dinucleotide_df.values())))
        final_junction_dinucleotide_df[gc_columns_dn] *= 100
        final_enc_df = add_whole_cds_calcs_enc_cai(pd.concat(list(final_enc_df.values())))
        
        final_codon_usage_df.to_csv(codon_out_file, index = False, sep = '\t')
        final_bicodon_usage_df.to_csv(bicodon_out_file, index = False, sep = '\t')
        final_dinucleotide_df.to_csv(dinucleotide_out_file, index = False, sep = '\t')
        final_junction_dinucleotide_df.to_csv(jn_dinucleotide_out_file, index = False, sep = '\t')
        final_enc_df.to_csv(enc_out_file, index = False, sep = '\t')
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        # final_syn_removed_df = split_gene[0]
        # wholeseqdf[wholeseqdf['terminal_lineage'] == 'referenceB']['sequence'].values[0]
        # print(final_syn_removed_df['orf7b_cds'].str[4:])
        
        