#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 09:49:27 2019

@author: jwalla12
"""

import pandas as pd
import glob2, argparse, os
import numpy as np
from pybiomart import Dataset
import sys

#htseq_directory='/Users/jwalla12/Repositories/htseq_counts_for_ck/'

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument( '-d', '--htseq-directory', dest='htseq_directory', nargs=1)
    #sys.argv.extend(['-d', '/Users/jwalla12/Repositories/htseq_counts_for_ck/'])
    return parser.parse_args()

def compile_htseq_files(htseq_directory):
    data_frames = []
    for file in (glob2.glob(os.path.join(htseq_directory, '*.htseq.counts.gz'))): 
        data_frames.append(pd.read_table(file, names=['gene', (file.strip('\n').split('/')[-1])]))        
    counts_matrix = reduce(lambda x, y: pd.merge(x, y, on='gene'), data_frames)
    return counts_matrix

def normalize_htseq_files(counts_matrix):
    counts_matrix_cleaned = counts_matrix[~counts_matrix['gene'].str.contains("__")] # get rid of rows where additional htseq output is reported
    counts_matrix_cleaned_indexed=counts_matrix_cleaned.set_index('gene')
    counts_matrix_norm = (counts_matrix_cleaned_indexed.div((counts_matrix_cleaned_indexed.apply(np.sum, axis=0))/1000000)).reset_index() 
    return counts_matrix_norm

def query_pybiomart():
    dataset = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
    gene_id_table = dataset.query(attributes=['ensembl_gene_id', 'ensembl_gene_id_version', 'hgnc_id', 'hgnc_symbol'])
    return gene_id_table

def append_gene_ids(counts_matrix, gene_id_table):
    count_ids = counts_matrix.merge(gene_id_table, how='left', left_on='gene', right_on='Gene stable ID version')
    counts_matrix_ids = count_ids.drop(labels=['Gene stable ID version'], axis=1).fillna(0)
    return counts_matrix_ids



def insert_descriptions(counts_matrix, htseq_directory, output_suffix):
    d=['#Column descriptions:', '#gene: ENSEMBL gene ID and version (as reported in GDC outputs)','#ensembl_gene_ID: ENSEMBL gene ID without version information','#hcng_symbol: HUGO Gene Nomenclature Committee (HGNC) gene symbol','#gene_id: HGNC gene ID', '#Column E:AW: UUIDs as per the sample manifest file (see manifest tab)', '#The raw tab is the raw read counts from each htseq file', '#The normalized tab is the raw counts normalized by library size per UUID -- that is to say (the number of reads mapped to each gene for a given UUID)/(millions of reads mapped across all genes for a given UUID)', '#Excluded from this calculation is the following additional information reported in the gene column of the raw counts', '#See the HTSeq documentation for more information: https://htseq.readthedocs.io/en/release_0.11.1/count.html', '#__no_feature: reads (or read pairs) which could not be assigned to any feature (set S as described above was empty).', '#__ambiguous: reads (or read pairs) which could have been assigned to more than one feature and hence were not counted for any of these, unless the --nonunique all option was used (set S had more than one element).', '#__too_low_aQual: reads (or read pairs) which were skipped due to the -a option', '#__not_aligned: reads (or read pairs) in the SAM file without alignment', '#__alignment_not_unique: reads (or read pairs) with more than one reported alignment. These reads are recognized from the NH optional SAM field tag. (If the aligner does not set this field, multiply aligned reads will be counted multiple times, unless they get filtered out by the -a option.) Note that if the --nonuniqueall option was used, these reads (or read pairs) are still assigned to features.']
    header_df = pd.DataFrame(data=d)
    #matrix_with_header = pd.concat([header_df, counts_matrix], ignore_index=True)
   # matrix_with_header = pd.concat([header_df, counts_matrix])
    matrix_with_header = pd.concat([header_df, counts_matrix], ignore_index=True)
    
    
    
    Gene_stable_ID = matrix_with_header['Gene stable ID']
    HGNC_ID = matrix_with_header['HGNC ID']
    HGNC_symbol = matrix_with_header['HGNC symbol']
    gene = matrix_with_header['gene']

    matrix_with_header = matrix_with_header.drop(['Gene stable ID', 'HGNC ID', 'HGNC symbol', 'gene'], axis=1)
    
    matrix_with_header.insert(0, 'HGNC symbol', HGNC_symbol)
    matrix_with_header.insert(0, 'HGNC ID', HGNC_ID)
    matrix_with_header.insert(0, 'Gene stable ID', Gene_stable_ID)
    matrix_with_header.insert(0, 'gene', gene)   
    
    matrix_with_header = matrix_with_header.set_index('gene')
    
    #matrix_with_header_fixed = matrix_with_header.set_index('gene')
    matrix_with_header.to_csv((htseq_directory+output_suffix)) #headers false too
    return matrix_with_header

#def concat_prepended_text():
#    args = get_args()
#    gene_id_table = query_pybiomart()
#    counts_matrix = compile_htseq_files(args.htseq_directory[0])
#    counts_matrix_norm = normalize_htseq_files(counts_matrix)
#    counts_matrix_with_ids = append_gene_ids(counts_matrix)
#    counts_matrix_norm_with_ids = append_gene_ids(counts_matrix_norm)
#    header = insert_descriptions()
#    
#    
#    append_gene_ids_norm()
#    insert_descriptions()
    
if __name__ == '__main__':
    args = get_args()
    counts_matrix = compile_htseq_files(args.htseq_directory[0])
    counts_matrix_norm = normalize_htseq_files(counts_matrix)
    gene_id_table = query_pybiomart()
    counts_matrix_with_ids = append_gene_ids(counts_matrix, gene_id_table)
    counts_matrix_norm_with_ids = append_gene_ids(counts_matrix_norm, gene_id_table)
    counts_matrix_with_ids_and_headers = insert_descriptions(counts_matrix_with_ids, args.htseq_directory[0], 'count_matrix.txt')
    counts_matrix_norm_with_ids_and_headers = insert_descriptions(counts_matrix_norm_with_ids, args.htseq_directory[0], 'count_matrix_norm.txt')
    

