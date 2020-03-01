#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from singlecellmultiomics.utils.submission import submit_job

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
    time = int( job_properties['resources']['time'])
except KeyError:
    try:
        time = int(job_properties['cluster']['time']) # runtime is time in hours
    except KeyError:
        try:
            time = job_properties['params']['runtime'].replace('h','')
        except KeyError:
            time = '12'


try:
    mem = job_properties['cluster']['mem']
except KeyError:
    try:
        mem = int( int(job_properties['resources']['mem_mb'])/1000 )
    except KeyError:
        mem = 10


# removing 'special' characters in log paths (default for snakemake)
base_path = os.path.dirname(job_script)

cluster_file_folder= base_path+'/cluster_jobs'

# formatting qsub command
job_id = submit_job(f'sh {job_script};', prefix='snake', target_directory=cluster_file_folder,  working_directory='.',
               threads_n=n, memory_gb=mem, time_h=time, scheduler='sge', copy_env=True,
               email=None, mail_when_finished=False, hold=None,submit=True)
print(job_id)
