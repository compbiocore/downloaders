#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 14:44:31 2019

@author: jwalla12
"""
#see https://docs.gdc.cancer.gov/API/Users_Guide/Python_Examples/#post-request-to-download-multiple-files

import requests
import re
import argparse
import time
start=time.time()

parser=argparse.ArgumentParser()
parser.add_argument( '-mf', '--manifest-file', dest='manifest_file', help="Give the full path to the manifest file")
parser.add_argument( '-ep', '--data-endpoint', dest='data_endpoint', help="Give the address to the API -- for example, 'https://api.gdc.cancer.gov/data'")
args=parser.parse_args()

#python tcga_simple.py -mf '/Users/jwalla12/Downloads/gdc_files_for_ck.txt' -ep 'https://api.gdc.cancer.gov/data/'

data_endpt = str(args.data_endpoint)
manifest = open((str(args.manifest_file)), 'r').readlines()

#manifest_file = '/Users/jwalla12/Downloads/gdc_files_for_ck.txt'
#manifest = open((str(manifest_file)), 'r').readlines()
#data_endpt = 'https://api.gdc.cancer.gov/data/'

uuids = []
for entry in manifest:
    uuids.append(entry.strip('\n').split('\t')[0])

for uuid in uuids:
    response=requests.get((data_endpt+uuid), headers={'Content-Type': 'application/json'})
    response_head_cd = response.headers['Content-Disposition']
    file_name = re.findall('filename=(.+)', response_head_cd)[0]
    with open(file_name, 'wb') as output_file:
        output_file.write(response.content)
        
end=time.time()
time = str(end-start)
print(time+' seconds to download files')

