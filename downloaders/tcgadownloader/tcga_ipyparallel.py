import sys, subprocess, argparse
from IPython.parallel import Client


def setup(manifest_file, token_file, add_uuid=False, md5_file="checksum.md5"):
    # manifest_file = "/Users/aragaven/Downloads/gdc_manifest_20180202_005135.txt"
    # token_file = "/Users/aragaven/Downloads/gdc-user-token.2018-03-06T18_16_48.240Z.txt"
    # md5_file = "/Users/aragaven/scratch/test_checksum.md5"

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
        com = com + " -H \"X-Auth-Token: " + token_txt
        com = com + "\" 'https://api.gdc.cancer.gov/data/" + uuid
        com = com + "' 2> logs/" + '_'.join(fname.split('.')[:-1]) + ".err "
        com = com + " 1>> logs/" + '_'.join(fname.split('.')[:-1]) + ".done;"
        com = com + " export ec=$?; let counter++;"
        com = com + " echo '" + fname + " Error Code:' $ec >> logs/" + '_'.join(fname.split('.')[:-1]) + ".curl_err;"
        com = com + " done; if [ $ec -ne 0 ] ;  then exit 1 ;"
        if add_uuid:
            fname = fname.replace(".bam" , "_" + uuid + ".bam")

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
    f = open('donwload_status.log', 'a')
    f.write(value)
    f.close()
    return value



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
    parser.add_argument( '-tf', '--token-file', help="Give the full path to the token file")
    parser.add_argument('-u','--add-uuid-to-filename', choices=[True,False], default=False)
    return parser.parse_args()

if __name__ == '__main__':
    # m = "/Users/aragaven/Downloads/gdc_manifest_20180202_005135.txt"
    # t = "/Users/aragaven/Downloads/gdc-user-token.2018-03-06T18_16_48.240Z.txt"
    # md5 = "/Users/aragaven/scratch/test_checksum.md5"
    # coms = setup(m,t,md5)
    my_args =get_args()
    coms = setup(my_args.manifest_file, my_args.token_file, uuid=my_args.add_uuid_to_filename)
    f = open('curl_cmds.log', 'w')
    for c in coms:
        f.write(c[1] + "\t" + c[0] + "\n")
    f.close()

    rc = Client()
    dview = rc[:]
    with dview.sync_imports():
        import sys, time, subprocess
    lview = rc.load_balanced_view()
    jobs_run = lview.map_async(run_curl, coms)
    jobs_run.wait_interactive(interval=60)
    # progress(jobs_run)

