#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import re
import itertools as it
import glob
import time
import datetime
import subprocess
import distutils.spawn
import uuid

def create_job_file_paths(target_directory,job_alias=None, prefix=None, job_file_name=None):

    if not os.path.exists(target_directory):
        os.makedirs(target_directory)

    if prefix is None:
        prefix = time.strftime("%d_%m_%Y_%H_%M_%S_") + str(uuid.uuid4())

    if job_file_name is None:
        job_file_name = '%s-%s' % (prefix, job_alias)

    jobfile = target_directory + '/%s.sh' % job_file_name
    stderr = target_directory + '/%s.stderr' % job_file_name
    stdout = target_directory + '/%s.stdout' % job_file_name

    if prefix is None:
        while os.path.exists(jobfile):
            job_file_name = '%s-%s' % (time.strftime("%d_%m_%Y_%H_%M_%S"), job_alias)
            jobfile = target_directory + '/%s.sh' % job_file_name
            stderr = target_directory + '/%s.stderr' % job_file_name
            stdout = target_directory + '/%s.stdout' % job_file_name
            time.sleep(1)
    else:
        if os.path.exists(jobfile):
            print(
                "Job %s already exists. Files might become corrupted if previous job is still running" %
                jobfile)

    return jobfile,stderr, stdout, job_file_name


def generate_job_script(scheduler, jobfile,stderr, stdout, job_name, memory_gb, working_directory, time_h, threads_n, email, mail_when_finished=False, copy_env=True, slurm_scratch_space_size=None ):
    if scheduler=='local':

        return [f'cd {working_directory}']

    if scheduler=='slurm':
        jobData = [
            '#!/bin/sh',
            '#SBATCH -J %s' % job_name, # Sets job name
            '#SBATCH -n %s' % threads_n,
            '#SBATCH -N 1', # Run on a single node
            '#SBATCH --time %s:00:00' % str(time_h).zfill(2),
            '#SBATCH --mem %sG' % memory_gb,
            '#SBATCH --chdir %s' % (working_directory),
            '#SBATCH -o %s' % stdout,
            '#SBATCH -e %s' % stderr
        ]

        if slurm_scratch_space_size is not None:
            jobData.append(f'#SBATCH --gres=tmpspace:{slurm_scratch_space_size}G')
        if email is not None:
            if mail_when_finished:
                raise NotImplementedError('email when finished is not implemented for slurm')

            jobData.append('#SBATCH --mail-type=FAIL')
            jobData.append('#SBATCH --mail-user=%s' % email)
    elif scheduler=='sge':
        jobData = [
            '#!/bin/sh',
            '#$ -S /bin/bash',
            '#$ -N %s' % job_name,
            '#$ -l h_rt=%s:00:00' % time_h,
            '#$ -l h_vmem=%sG' % memory_gb,
            # '#$ -l hostname=\'!n00[18]*\'',
            '#$ -wd %s' % (working_directory),
            '#$ -o %s' % stdout,
            '#$ -e %s' % stderr,
            '#$ -q all.q'
        ]

        if email is not None:
            jobData.append('#$ -M %s' % email)
            jobData.append('#$ -m %sas' % ('e' if mail_when_finished else ''))

        if copy_env:
            jobData.append('#$ -V')

        if threads_n > 1:
            jobData.append('#$ -pe threaded %s' % threads_n)

    # Make sure we land in the right directory
    if working_directory is not None:
        jobData.append(f'cd {working_directory}')
    return jobData

def write_cmd_to_submission_file(cmd, job_data, jobfile, scheduler='sge' ):
    if scheduler in ('slurm','sge','local'):
        job_data.append('%s' % cmd)
    else:
        raise NotImplementedError()

    with open(jobfile, 'w') as f:
        f.write('\n'.join(job_data) + '\n')

def generate_submission_command(jobfile, hold, scheduler='sge'):

    if scheduler=='slurm':
        if hold is not None and len(hold)>0 and hold[0]!='none':

            js = 'afterany:' + ':'.join( [f'{h.strip()}' for h in hold] )
            qs = f'sbatch --dependency={js} {jobfile}'

        else:
            qs = 'sbatch %s' % jobfile
    else:
        qs = 'qsub %s %s' % ((('-hold_jid %s' % ','.join(hold))
                              if (hold is not None and hold[0] != 'none') else ''), jobfile)
    return qs


def submit_job(command,  target_directory,  working_directory,
               threads_n=1, memory_gb=8, time_h=8, scheduler='sge', copy_env=True,
               email=None,job_alias=None, mail_when_finished=False,
               hold=None,submit=True, prefix=None, job_file_name=None, job_name=None, silent=False, slurm_scratch_space_size=None):
    """
    Submit a job

    Args:
        threads(int) : amount of requested threads
        memory_gb(int) : amount of requested memory
        scheduler(str): sge/slurm/local
        hold(list): list of job depedencies
        submit(bool) : perform the actual submission, when set to False only the submission script is written
    Returns:
        job_id(str) : id of sumbitted job
    """

    qsub_available = (distutils.spawn.find_executable("qsub") is not None)
    sbatch_available = (distutils.spawn.find_executable("sbatch") is not None)

    if scheduler == 'auto':
        if qsub_available:
            scheduler = 'sge'
        elif sbatch_available:
            scheduler = 'slurm'
        else:
            scheduler = 'local'


    if job_alias is None and job_name is None:
        job_name = 'J%s' % str(uuid.uuid4())

    # If no file was specified, we generate a file using the supplied job name
    if job_file_name is None:
        job_alias = job_name


    if working_directory is None:
        working_directory = os.getcwd()

    if submit:
        if scheduler=='sge' and not qsub_available:
            raise ValueError('qsub is not available on the system')
        if scheduler=='slurm' and not sbatch_available:
            if qsub_available:
                print('SBATCH is not available, but QSUB is, reverting to use QSUB')
                scheduler='sge'
            else:
                raise ValueError('sbatch is not available on the system')

    jobfile,stderr, stdout, _job_file_name = create_job_file_paths(target_directory,job_alias=job_alias,prefix=prefix,job_file_name=job_file_name)
    if job_file_name is None:
        job_file_name=_job_file_name
    else:
        if job_file_name!=_job_file_name and not silent:
            print(f'Job file name changed from {job_file_name} to {_job_file_name}')

    job_data = generate_job_script(scheduler=scheduler, jobfile=jobfile,
    stderr=stderr, stdout=stdout,
    job_name=job_name,
    memory_gb=memory_gb, working_directory=working_directory,
    time_h=time_h, threads_n=threads_n,
    email=email, mail_when_finished=mail_when_finished,
     copy_env= copy_env,slurm_scratch_space_size=slurm_scratch_space_size)

    qs = generate_submission_command( jobfile, hold, scheduler)
    write_cmd_to_submission_file(command, job_data, jobfile, scheduler)

    if submit:
        if scheduler=='slurm':
            job_id = os.popen(qs).read().replace('Submitted batch job ','').strip()
            return job_id
        elif scheduler=='sge':
            rd  = os.popen(qs).read()
            job_id = rd.split(' ')[2]
            return job_id.strip()
        elif scheduler=='local':
            # Run the job now:
            os.system(f'bash {jobfile} 2>{stderr} >{stdout}')

    else:
        print('# use the command below to submit your job:')
        print(qs)


## ##### Dependency handling  ##### ##
if __name__ == '__main__':
    import argparse
    username = os.getenv('USER')
    defaultEmail = os.getenv('EMAIL')


    qsub_available = (distutils.spawn.find_executable("qsub") is not None)

    PY36ENV = os.getenv('PY36ENV')
    if PY36ENV is None:
        PY36ENV = 'source /hpc/hub_oudenaarden/bdebarbanson/virtualEnvironments/py36/bin/activate'

    basepath = '/hpc/hub_oudenaarden/%s/cluster' % username if os.path.exists(
        '/hpc/hub_oudenaarden/%s' % username) else os.path.dirname(
        os.path.abspath(__file__)) + '/cluster/'
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Send job to cluster')
    argparser.add_argument(
        '-w', help="Working directory, current dir when not supplied")
    argparser.add_argument(
        '-N',
        type=str,
        help="Job alias, slashes will be removed",
        default="%s-dobby" %
        (username if username is not None else 'anon'))
    argparser.add_argument(
        '-jp',
        type=str,
        help="Job file prefix, this will be the name of the file where the standard error, standard out and job script will be saved. By default it is the current time to prevent collisions. When supplying this be careful that simultanious jobs cannot share the same prefix!",
        default=None)

    argparser.add_argument(
        '-t',
        type=int,
        default=1,
        help="Threads, amount of CPUs requested (PE). Cluster Worker count")
    argparser.add_argument(
        '-time',
        type=int,
        default=24,
        help="Runtime in hours")
    argparser.add_argument('-m', type=int, default=4, help="Memory in gigs")
    argparser.add_argument('-y', action="store_true", help="Submit jobs")

    argparser.add_argument(
        '-sched',
        default='slurm',
        help="scheduler: sge, slurm, local")


    argparser.add_argument(
        '-e',
        type=str,
        help="How to execute the job; submit, local",
        default="submit")
    argparser.add_argument(
        '-hold',
        type=str,
        help="Wait for job(s) with this name to be finished",
        default=None)
    argparser.add_argument(
        '-email',
        type=str,
        help="Send emails to this adress (by default Only kill messages)",
        default=os.getenv('EMAIL'))
    argparser.add_argument(
        '--mf',
        help="Mail when finished",
        action='store_true')
    argparser.add_argument(
        '--nenv',
        help="Do not copy current environment",
        action='store_true')
    argparser.add_argument(
        '--py36',
        help="Source python 3.6 (set PY36ENV variable to change the path)",
        action='store_true')

    argparser.add_argument(
        '-sz',
        help="SLURM Scratch space request in GB",
        type=int)

    argparser.add_argument('c', metavar='command', type=str, nargs='*')
    argparser.add_argument(
        '--silent',
        help="Try to print less",
        action='store_true')
    argparser.add_argument(
        '-s',
        type=str,
        help="Submission data storage path (stdout/stderr)",
        default=os.path.abspath('./cluster'))
    args = argparser.parse_args()

    if args.email == 'none':
        args.email = None

    working_directory = args.w if args.w is not None else os.getcwd()

    jid = submit_job(' '.join(args.c), job_name=args.N, target_directory=args.s,
                    job_file_name = args.jp,
                   working_directory=working_directory,
                   threads_n=args.t, memory_gb=args.m, time_h=args.time, scheduler=args.sched, copy_env=not args.nenv,
                   email=args.email, mail_when_finished=args.mf, hold=(args.hold.split(',') if args.hold is not None else None),
                   submit=args.y,
                   prefix=None,
                   slurm_scratch_space_size=args.sz)
    if jid is not None:
        print(jid)
