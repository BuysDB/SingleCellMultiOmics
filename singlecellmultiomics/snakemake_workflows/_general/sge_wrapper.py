#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import re
import subprocess
import uuid
from snakemake.utils import read_job_properties

# load
## loading job stdout & stderr path
log_path = sys.argv[-2]
## loading job script (provided by snakemake)
job_script = sys.argv[-1]
job_properties = read_job_properties(job_script)

# getting job parameters from snakemake-generated job script
try:
    threads = job_properties['threads']
except KeyError:
    threads = 1
n = threads

try:
    time = job_properties['cluster']['time'] # runtime is time in hours
except KeyError:
    try:
        time = job_properties['params']['runtime'].replace('h','') + ':00:00'
    except KeyError:
        time = '12:00:00'


try:
    mem = job_properties['cluster']['mem']
except KeyError:
    try:
        mem = int( int(job_properties['resources']['mem_mb'])/1000 )
    except KeyError:
        mem = 10


# removing 'special' characters in log paths (default for snakemake)
base_path = os.path.dirname(job_script)

try:
    os.makedirs(base_path+'/cluster_jobs')
except Exception as e :
    pass

job_id = uuid.uuid4()
std_out =  f'{base_path}/cluster_jobs/{job_id}.o'
std_err = f'{base_path}/cluster_jobs/{job_id}.e' 


# formatting qsub command
cmd = "qsub -V -pe threaded {n} -l h_vmem={mem}G -l h_rt={time} -o {std_out} -e {std_err} {job_script}"
cmd = cmd.format(n=n, mem=mem, time=time, std_out=std_out, std_err=std_err, job_script=job_script)

# subprocess job: qsub
try:
    res = subprocess.run(cmd, check=True, shell=True, stdout=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    raise e

# get qsub job ID
res = res.stdout.decode()
try:
    m = re.search("Your job (\d+)", res)
    jobid = m.group(1)
    print(jobid)
except Exception as e:
    print(e)
    raise
