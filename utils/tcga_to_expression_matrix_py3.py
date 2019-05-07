#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  6 09:30:49 2019

@author: jwalla12
"""

import pandas as pd
import glob2, argparse, os
from functools import reduce
import numpy as np


#import sys

htseq_directory='/Users/jwalla12/Repositories/htseq_counts_for_ck/'
gene_id_table = pd.read_csv('/Users/jwalla12/Repositories/htseq_counts_for_ck/gene_ids.txt', sep='\t')

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument( '-d', '--htseq-directory', dest='htseq_directory', nargs=1)
    #sys.argv.extend(['-d', '/Users/jwalla12/Repositories/htseq_counts_for_ck/'])
    return parser.parse_args()

def compile_htseq_files(htseq_directory):
    data_frames = []
    for file in (glob2.glob(os.path.join(htseq_directory, '*.htseq.counts.gz'))): 
        data_frames.append(pd.read_csv(file, names=['gene', (file.strip('\n').split('/')[-1])], sep='\t'))        
    counts_matrix = reduce(lambda  left,right: pd.merge(left,right,on=['gene'],  how='outer'), data_frames)
    return counts_matrix

def normalize_htseq_files(counts_matrix):
    counts_matrix_cleaned = counts_matrix[~counts_matrix['gene'].str.contains("__")] # get rid of rows where additional htseq output is reported
    counts_matrix_cleaned_indexed=counts_matrix_cleaned.set_index('gene')
    #counts_matrix_cleaned_indexed_dropped = counts_matrix_cleaned_indexed.drop(labels=['base_gene_id','version_gene_id', 'Gene description', 'Gene stable ID','Gene stable ID version', 'Transcript stable ID','Transcript stable ID version', 'HGNC ID', 'HGNC symbol', 'HGNC transcript name ID', 'Gene name'], axis=1)
    counts_matrix_norm = (counts_matrix_cleaned_indexed.div((counts_matrix_cleaned_indexed.apply(np.sum, axis=0))/1000000)).reset_index()
    #counts_matrix_norm.to_csv('/Users/jwalla12/Repositories/htseq_counts_for_ck/norm_counts_2.txt', sep='\t') 
    return counts_matrix_norm

def append_gene_ids(counts_matrix, gene_id_table):
    #make a new column called base_gene_id and populate it with the 1st entry of the . splitted gene columns
    gene_ids = pd.read_csv(gene_id_table, sep='\t')
    counts_matrix[['base_gene_id', 'version_gene_id']] = counts_matrix['gene'].str.split('.', expand=True)
    count_ids = counts_matrix.merge(gene_ids, how='left', left_on='base_gene_id', right_on='ensembl_gene_id')
    return count_ids

def rename_headers(count_ids, htseq_directory):
    path = (glob2.glob(os.path.join(htseq_directory, 'manifest.txt')))[0]
    u=open(path, "r+")
    u_stripped = [l.rstrip() for l in u.readlines()]
    d={}
    for item in u_stripped:
        uuid, file, md5, size, status = item.split('\t')
        d[file]=(uuid)
    count_ids=count_ids.rename(columns=d)
    count_ids.set_index('gene')
    return(count_ids)


def insert_descriptions(count_ids, htseq_directory):
    file = (glob2.glob(os.path.join(htseq_directory, 'new_headers.txt')))[0]
    header_df = pd.read_csv(file, sep='\t')
    matrix_with_header = pd.concat([header_df, count_ids]).reset_index()
    
    matrix_with_header_ordered = matrix_with_header[['#Columns:', 'gene', 'ensembl_gene_id', 'ensembl_gene_id_version', 'description', 'hgnc_id', 'hgnc_symbol', 'hgnc_trans_name','d9249d25-d852-4008-a71e-dd9be6cb68e9', '03aee74e-4e37-4a58-a720-c90e807d2f40', '7845948f-701e-49c5-8b76-2f0e2f0d5a76', 'c32379fa-928f-4eec-b849-a5b715d1f9cb', '0d2c466e-d8b8-4b8f-9909-0be2175fa6a0', '3a39302b-b097-4371-a785-753c8dabb0a2', '4e337123-ea4c-42e5-acb9-37fce1fa7fd5', '7a5b1cfb-e4d1-4989-bfc0-6eb2afb31e2e', '7dd1c33e-35ef-406d-b516-5b3b2a486473', '533b56af-6953-4e39-b23e-ddf9a06499a5', '42994602-7f8d-4d05-a191-9bcec686aaf7', '494626b6-541d-417f-8ad0-1c4c44e25bef', '10d08172-48d2-49e7-b760-721163492cc1', '9796011a-fbab-409a-a850-bd890098c576', 'f18234f6-9933-4aa7-8c17-f41f64569cf6', 'b223180f-7a91-40dd-a3f1-ab423a820255', 'f5cf4d3b-50f0-4c39-a8c3-450436f5844a', '4963b9af-7d7b-42cd-b057-2f894639e59f', '8e06ea1c-dbca-4b01-8f00-e4259507a77a', 'efb2d3a5-b73d-4895-865c-4c27d4c142bc', '4fc729ce-f4f2-492f-b305-10f51ca60972', '361cb97e-20ee-4a03-8c77-574fe10bf0d5', 'bd31db54-5544-41e6-abab-86f2e997cfe5', '8499d831-d9d9-4279-8fbe-84d48a562ed0', '85bc7f81-51fb-4446-b12d-8741eef4acee', '19e8fd21-f6c8-49b0-aa76-109eef46c2e9', '42b8d463-6209-4ea0-bb01-8023a1302fa0', 'be2d52c8-eb14-4b45-a65b-c1bdd6f5c4c2', 'b74b20b7-2797-4266-8f56-166978275576', 'b48c636d-ac52-40cd-8279-0215259f8406', 'b886634e-e0a7-459a-a2c9-20fabb779dc6', '13921de8-8580-4f26-88bb-49f71ea6e103', 'b554b9e8-01c4-46ad-8909-5ada19fe5860', 'afecdda2-735c-4304-a087-ef917ad9cd5a', 'f38840ce-3838-40df-b5f2-ecd58bfe7c57', 'e4e4df95-84e9-4811-9ff9-d3bf4f3be116', 'fd0ea67b-5b75-471f-be3c-a92142b91cf3', '6e2031e9-df75-48df-b094-8dc6be89bf8b', '1ace0df3-9837-467e-85de-c938efda8fc8', 'e7a3c801-4e2f-45fd-948b-cef353403afb', '5121695c-6f2c-4fd5-b75e-25d75d1aa567', '50e857b6-77be-4dd6-ac94-4b2186129fcb', '8e31c3bb-40c8-48e2-9de1-1ee1ee31aac8', '7cde9495-e573-4b38-b89c-991076cf8cf8', '3f34af32-53c9-4b97-8886-a3057991c466']]
    
    return matrix_with_header_ordered

def clinical_data(htseq_directory, clinical_data, manifest, json):
    
    
    clinical_data_file = pd.read_csv((os.path.join(htseq_directory, clinical_data))), sep='\t')
    original_manifest = (os.path.join(htseq_directory, manifest))
    json_file = (os.path.join(htseq_directory, json))
    
    
    
    path = (glob2.glob(os.path.join(htseq_directory, 'manifest.txt')))[0]
    
    clinical_data_df = pd.read_csv(glob2.glob(os.path.join(htseq_directory, 'clinical_data.txt')))[0]
    
    
    clinical_data_df = pd.read_csv('/Users/jwalla12/Repositories/htseq_counts_for_ck/clinical_data/clinical_data.txt', sep='\t')


clinical_data = pd.read_csv('/Users/jwalla12/Repositories/htseq_counts_for_ck/clinical_data/clinical_data.txt', sep='\t')
original_manifest = pd.read_csv('/Users/jwalla12/Repositories/htseq_counts_for_ck/clinical_data/gdc_files_for_ck.txt', sep='\t')
cleaned_manifest_json = pd.read_csv('/Users/jwalla12/Repositories/htseq_counts_for_ck/clinical_data/manifest_df_3_clean.csv', sep=',')

#merge the clinical data to the json data
clinical_data_json = clinical_data.merge(cleaned_manifest_json, how='outer', on='case_id')

clinical_data_json_original = original_manifest.merge(clinical_data_json, how = 'outer', left_on='file_uuid', right_on='file_id')

clinical_data_json_original.to_csv('/Users/jwalla12/Repositories/htseq_counts_for_ck/clinical_data/ck_data_clinical.csv', sep=',')



def excel_writer(raw_counts, normalized_counts):
    writer = pd.ExcelWriter((os.path.join(htseq_directory, 'counts.xlsx')), engine='xlsxwriter')
    raw_counts.to_excel(writer, sheet_name='raw_counts')
    normalized_counts.to_excel(writer, sheet_name='normalized_counts')
    writer.save()
    
if __name__ == '__main__':
    args = get_args()    
    
    counts = compile_htseq_files(args.htseq_directory[0])
    normalized_counts = normalize_htseq_files(counts)
        
    counts_ids = append_gene_ids(counts, gene_id_table = '/Users/jwalla12/Repositories/htseq_counts_for_ck/gene_ids.txt')
    normalized_counts_ids = append_gene_ids(normalized_counts, gene_id_table = '/Users/jwalla12/Repositories/htseq_counts_for_ck/gene_ids.txt')
    
    counts_ids_renamed = rename_headers(counts_ids, args.htseq_directory[0])
    normalized_counts_ids_renamed = rename_headers(normalized_counts_ids, args.htseq_directory[0])
    
    counts_ids_renamed_headers = insert_descriptions(counts_ids_renamed, args.htseq_directory[0])
    normalized_counts_ids_renamed_headers = insert_descriptions(normalized_counts_ids_renamed, args.htseq_directory[0])

    



    counts = compile_htseq_files(htseq_directory)
    normalized_counts = normalize_htseq_files(counts)
        
    counts_ids = append_gene_ids(counts, gene_id_table = '/Users/jwalla12/Repositories/htseq_counts_for_ck/gene_ids.txt')
    normalized_counts_ids = append_gene_ids(normalized_counts, gene_id_table = '/Users/jwalla12/Repositories/htseq_counts_for_ck/gene_ids.txt')
    
    counts_ids_renamed = rename_headers(counts_ids, htseq_directory)
    normalized_counts_ids_renamed = rename_headers(normalized_counts_ids, htseq_directory)
    
    counts_ids_renamed_headers = insert_descriptions(counts_ids_renamed, htseq_directory)
    normalized_counts_ids_renamed_headers = insert_descriptions(normalized_counts_ids_renamed, htseq_directory)    
        
    excel_writer(counts_ids_renamed_headers, normalized_counts_ids_renamed_headers)
