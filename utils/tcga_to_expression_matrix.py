#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 09:49:27 2019

@author: jwalla12
"""

import pandas as pd
import glob2
import numpy as np
import os
import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument( '-d', '--htseq-directory')
    return parser.parse_args()

def expression_matrices(file_directory):    
    data_frames = []
    for file in (glob2.glob(os.path.join(file_directory, '*.htseq.counts.gz'))): 
        data_frames.append(pd.read_table(file, names=['gene', (file.strip('\n').split('/')[-1])]))        
    raw_counts_expression_matrix = reduce(lambda x, y: pd.merge(x, y, on='gene'), data_frames)    
    raw_counts_expression_matrix.to_csv((os.path.join(file_directory,'raw_counts_expression_matrix.txt')), sep=',', index=False)   
    raw_counts_expression_matrix_only = raw_counts_expression_matrix[~raw_counts_expression_matrix['gene'].str.contains("__")]   
    raw_counts_expression_matrix_only.set_index('gene', inplace=True) 
    normalized_counts_expression_matrix = raw_counts_expression_matrix_only.div((raw_counts_expression_matrix_only.apply(np.sum, axis=0))/1000000) 
    normalized_counts_expression_matrix.to_csv((os.path.join(file_directory,'normalized_counts_expression_matrix.txt')), sep=',', index=True)

if __name__ == '__main__':
    my_args = get_args()
    expression_matrices(my_args.htseq_directory)
    
#expression_matrices('/Users/jwalla12/Repositories/htseq_counts_for_ck/')    