import sys, subprocess, argparse, requests, re, tqdm
#from IPython.parallel import Client
from ipyparallel import Client

#setup(manifest_file='/Users/jwalla12/Downloads/gdc_files_for_ck.txt', md5_file='/Users/jwalla12/Downloads/checksum.md5')

#maybe eventually make file types arguments that are passed to com (so .bam, .htseq.counts.gz, etc.)
def setup(manifest_file, token_file=None, add_uuid=True, md5_file="checksum.md5"):
    # manifest_file = "/Users/aragaven/Downloads/gdc_manifest_20180202_005135.txt"
    # token_file = "/Users/aragaven/Downloads/gdc-user-token.2018-03-06T18_16_48.240Z.txt"
    # md5_file = "/Users/aragaven/scratch/test_checksum.md5"
    if token_file is not None:
        token_txt = open(token_file, 'r').readlines()[0].strip('\n')
    manifest = open(manifest_file, 'r').readlines()
    md5_file = open(md5_file, 'w')
    coms_to_submit = []
    
    for ln in manifest[1:]:
        uuid, fname, md5, size, state = ln.strip('\n').split()
        # com = "trap 'exit 42' SIGINT SIGQUIT SIGTERM ; "
        com = ""
        com = com + "export ec=56; counter=1; while [[ $ec -eq 56 && $counter -lt 500 ]]; do "
        com = com + "curl --no-buffer --show-error -O -C - --max-time 3600 "
        if token_file is not None:
            com = com + " -H \"X-Auth-Token: " + token_txt +  "\ " 
        com = com + "'https://api.gdc.cancer.gov/data/" + uuid
        com = com + "' 2> logs/" + '_'.join(fname.split('.')[:-1]) + ".err "
        com = com + " 1>> logs/" + '_'.join(fname.split('.')[:-1]) + ".done;"
        com = com + " export ec=$?; let counter++;"
        com = com + " echo '" + fname + " Error Code:' $ec >> logs/" + '_'.join(fname.split('.')[:-1]) + ".curl_err;"
        com = com + " done; if [ $ec -ne 0 ] ;  then exit 1 ;"
        if add_uuid:
            fname = fname.replace(".htseq.counts.gz" , "_" + uuid + ".htseq.counts.gz")

        com = com + "else  mv " + uuid + " " + fname + " ; fi"
        md5_file.write(md5 + "  " + fname + "\n")
        coms_to_submit.append((com, fname))

    md5_file.close()
    return coms_to_submit



def run_curl(com):
    cmd, fname = com
    ret_val = subprocess.call(cmd, shell=True)
    value = "\n***** success: " + fname + " *****\n" + cmd + "\n"
    if ret_val != 0:
        value = "\n**** failed: " + fname + " *****\n" + cmd + "\n"
        #put the simple version of events right here...
        run_requests(my_args.manifest_file)
    f = open('download_status.log', 'a')
    f.write(value)
    f.close()
    return value

def run_requests(manifest_file):  
    uuids = []
    manifest = open(manifest_file, 'r').readlines()
    for entry in manifest:
        uuids.append(entry.strip('\n').split('\t')[0])
    
    for uuid in uuids:
        response=requests.get(('https://api.gdc.cancer.gov/data/'+uuid), headers={'Content-Type': 'application/json'})
        response_head_cd = response.headers['Content-Disposition']
        file_name = re.findall('filename=(.+)', response_head_cd)[0]
        with open(file_name, 'wb') as output_file:
            output_file.write(response.content)
            

def progress(ar):  # count, total, status=''
    total = len(ar.msg_ids)
    bar_len = 60
    status = "downloading"
    while not ar.ready():
        count = ar.progress
        filled_len = int(round(bar_len * count / float(total)))

        percents = round(100.0 * count / float(total), 1)
        bar = '=' * filled_len + '-' * (bar_len - filled_len)

        sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
        sys.stdout.flush()
        time.sleep(10)
    status = "complete"
    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
    sys.stdout.flush()
    return


def get_args():
    parser = argparse.ArgumentParser(description='download tcga data using curl and IPy Parallel')
    parser.add_argument( '-mf', '--manifest-file', help="Give the full path to the manifest file")
    parser.add_argument( '-tf', '--token-file', help="Give the full path to the token file", default=None)
    parser.add_argument('-u','--add-uuid-to-filename', choices=[True,False], default=True)
    parser.add_argument('-d','--dry-run', choices=['True','False'], default='False')
    parser.add_argument('-p','--parallel-run', choices=[True,False], default=False)
    return parser.parse_args()

if __name__ == '__main__':
    # m = "/Users/aragaven/Downloads/gdc_manifest_20180202_005135.txt"
    # t = "/Users/aragaven/Downloads/gdc-user-token.2018-03-06T18_16_48.240Z.txt"
    # md5 = "/Users/aragaven/scratch/test_checksum.md5"
    # coms = setup(m,t,md5)
    my_args = get_args()
    if not my_args.token_file:
        coms = setup(my_args.manifest_file, add_uuid=my_args.add_uuid_to_filename)
    else:
        coms = setup(my_args.manifest_file, my_args.token_file, add_uuid=my_args.add_uuid_to_filename)
   
    
    f = open('curl_cmds.log', 'w')
    for c in coms:
        f.write(c[1] + "\t" + c[0] + "\n")
    f.close()
    
    if my_args.dry_run == 'False':
        print "running job"
        if not my_args.parallel_run:
            for c in tqdm.tqdm(coms):
                #print "running:\n"
                #print c
                subprocess.check_output(c, shell=True)
        else:
            rc = Client()
            dview = rc[:]
            with dview.sync_imports():
                import time, subprocess
            lview = rc.load_balanced_view()
            jobs_run = lview.map_async(run_curl, coms)
            jobs_run.wait_interactive(interval=60)
    # progress(jobs_run)

