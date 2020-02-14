#!/usr/bin/env python3
import sys
import os
import re
import itertools as it
import glob
import time
import datetime
import subprocess
import distutils.spawn

def create_job_file_paths(target_directory,job_alias='dobby', jobName=None):

    if not os.path.exists(target_directory):
        os.makedirs(target_directory)

    if jobName is None:
        jobName = '%s-%s' % (time.strftime("%d_%m_%Y_%H_%M_%S"), job_alias)

    jobfile = args.s + '/%s.sh' % jobName
    stderr = args.s + '/%s.stderr' % jobName
    stdout = args.s + '/%s.stdout' % jobName

    if args.jp is None:
        while os.path.exists(jobfile):
            jobName = '%s-%s' % (time.strftime("%d_%m_%Y_%H_%M_%S"), args.N)
            jobfile = args.s + '/%s.sh' % jobName
            stderr = args.s + '/%s.stderr' % jobName
            stdout = args.s + '/%s.stdout' % jobName
            time.sleep(1)
    else:
        if os.path.exists(jobfile):
            print(
                "Job %s already exists. Files might become corrupted if previous job is still running" %
                jobfile)

    return jobfile,stderr, stdout


def generate_job_script(manager, job_name, memory_gb, working_directory, time_h, threads_n, email, mail_when_finished=False )
    
    if manager=='slurm':
        jobData = [
            '#!/bin/sh',
            '#SBATCH -J %s' % args.N, # Sets job name
            '#SBATCH -n %s' % threads_n,
            '#SBATCH --time %s:00:00' % str(time_h).zfill(2),
            '#SBATCH --mem %sG' % memory_gb,
            '#SBATCH --chdir %s' % (working_directory),
            '#SBATCH -o %s' % stdout,
            '#SBATCH -e %s' % stderr
        ]

        if email is not None:
            jobData.append('#SBATCH --mail-type=FAIL')
            jobData.append('#SBATCH --mail-user=%s' % email)
    elif manager=='sge':
        jobData = [
            '#!/bin/sh',
            '#$ -S /bin/bash',
            '#$ -N %s' % job_name,
            '#$ -l h_rt=%s:00:00' % time_h,
            '#$ -l h_vmem=%sG' % memory_gb,
            # '#$ -l hostname=\'!n00[18]*\'',
            '#$ -wd %s' % (workingDirectory),
            '#$ -o %s' % stdout,
            '#$ -e %s' % stderr,
            '#$ -q all.q'
        ]

        if args.email is not None:
            jobData.append('#$ -M %s' % email)
            jobData.append('#$ -m %sas' % ('e' if mail_when_finished else ''))

        if not args.nenv:
            jobData.append('#$ -V')

        if args.t > 1:
            jobData.append('#$ -pe threaded %s' % threads_n)

    return jobData


#jobData.append( 'export PYTHONPATH="";\nsource /hpc/hub_oudenaarden/bdebarbanson/virtualEnvironments/py36/bin/activate;')
#jobData.append('module load python%s' % args.v)

if args.slurm:
    jobData.append('%s' % cmd)
else:
    jobData.append('%s' % cmd)
if not args.silent:
    print('\n'.join(jobData))

with open(jobfile, 'w') as f:
    f.write('\n'.join(jobData) + '\n')

if args.slurm:
    qs = 'sbatch %s %s' % ((('-d %s' % args.hold)
                          if (args.hold is not None and args.hold != 'none') else ''), jobfile)
else:
    qs = 'qsub %s %s' % ((('-hold_jid %s' % args.hold)
                          if (args.hold is not None and args.hold != 'none') else ''), jobfile)



def submit_job(command, job_alias, target_directory, threads=1, memory_gb=8, time_h=8, manager='sge', mail_when_finished=False, hold=None,submit=True):
    """
    Submit a job

    Args:
        threads(int) : amount of requested threads
        memory_gb(int) : amount of requested memory
        manager(str): sge/slurm/local
        hold(list): list of job depedencies
        submit(bool) : perform the actual submission, when set to False only the submission script is written
    Returns:
        job_id(str) : id of sumbitted job
    """

    qsub_available = (distutils.spawn.find_executable("qsub") is not None)

    sbatch_available = (distutils.spawn.find_executable("sbatch") is not None)

    jobfile,stderr, stdout = create_job_file_paths(target_directory,job_alias)

    create_job_script()



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
        '--slurm',
        action="store_true",
        help="Send to SLURM scheduler")

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


    if args.i:
        for job in glob.glob("%s/cluster/*.sh" % args.s):
            jobName = os.path.basename(job).replace('.sh', '')
            print(jobName)
            t = datetime.datetime.strptime(
                jobName.split('-')[0], "%d_%m_%Y_%H_%M_%S")

            print(
                '%s executed %s ago' %
                (jobName.split(
                    '-',
                    1)[1],
                    datetime.datetime.now() -
                    t))
            if os.path.exists('%s.stderr' % job.replace('.sh', '')):
                print('Running or finished')
            else:
                print('Crashed or queued')
        exit()
    cmd = ' '.join(args.c)
    print('Command to execute: %s' % cmd)

    if not os.path.exists(args.s):
        os.makedirs(args.s)
    if args.jp is None:
        jobName = '%s-%s' % (time.strftime("%d_%m_%Y_%H_%M_%S"), args.N)
    else:
        jobName = args.jp
    jobfile = args.s + '/%s.sh' % jobName
    stderr = args.s + '/%s.stderr' % jobName
    stdout = args.s + '/%s.stdout' % jobName

    if args.jp is None:
        while os.path.exists(jobfile):
            jobName = '%s-%s' % (time.strftime("%d_%m_%Y_%H_%M_%S"), args.N)
            jobfile = args.s + '/%s.sh' % jobName
            stderr = args.s + '/%s.stderr' % jobName
            stdout = args.s + '/%s.stdout' % jobName
            time.sleep(1)
    else:
        if os.path.exists(jobfile):
            print(
                "Job %s already exists. Files might become corrupted if previous job is still running" %
                jobfile)
    workingDirectory = args.w if args.w is not None else os.getcwd()

    # '#$ -r yes', # This allows jobs to be rescheduled
    if args.slurm:

        jobData = [
            '#!/bin/sh',
            '#SBATCH -J %s' % args.N, # Sets job name
            '#SBATCH -n %s' % args.t,
            '#SBATCH --time %s:00:00' % str(args.time).zfill(2),
            '#SBATCH --mem %sG' % args.m,
            '#SBATCH --chdir %s' % (workingDirectory),
            '#SBATCH -o %s' % stdout,
            '#SBATCH -e %s' % stderr
        ]

        if args.email is not None:
            jobData.append('#SBATCH --mail-type=FAIL')
            jobData.append('#SBATCH --mail-user=%s' % args.email)
    else:
        jobData = [
            '#!/bin/sh',
            '#$ -S /bin/bash',
            '#$ -N %s' % args.N,
            '#$ -l h_rt=%s:00:00' % args.time,
            '#$ -l h_vmem=%sG' % args.m,
            # '#$ -l hostname=\'!n00[18]*\'',
            '#$ -wd %s' % (workingDirectory),
            '#$ -o %s' % stdout,
            '#$ -e %s' % stderr,
            '#$ -q all.q'
        ]

        if args.email is not None:
            jobData.append('#$ -M %s' % args.email)
            jobData.append('#$ -m %sas' % ('e' if args.mf else ''))

        if not args.nenv:
            jobData.append('#$ -V')

        if args.t > 1:
            jobData.append('#$ -pe threaded %s' % args.t)

    if not args.silent:
        print('using %s as working directory' % workingDirectory)
    if args.py36:
        jobData.append(PY36ENV)
    jobData.append('cd %s' % workingDirectory)

    #jobData.append( 'export PYTHONPATH="";\nsource /hpc/hub_oudenaarden/bdebarbanson/virtualEnvironments/py36/bin/activate;')
    #jobData.append('module load python%s' % args.v)

    if args.slurm:
        jobData.append('%s' % cmd)
    else:
        jobData.append('%s' % cmd)
    if not args.silent:
        print('\n'.join(jobData))

    with open(jobfile, 'w') as f:
        f.write('\n'.join(jobData) + '\n')

    if args.slurm:
        qs = 'sbatch %s %s' % ((('-d %s' % args.hold)
                              if (args.hold is not None and args.hold != 'none') else ''), jobfile)
    else:
        qs = 'qsub %s %s' % ((('-hold_jid %s' % args.hold)
                              if (args.hold is not None and args.hold != 'none') else ''), jobfile)
    if not args.silent:
        print(qs)
    if args.y:

        if args.e == 'qsub' and qsub_available:

            if args.slurm:
                job_id = os.popen(qs).read().replace('Submitted batch job ','')
            else:
                rd  = os.popen(qs).read()
                job_id = rd.split(' ')[2]

            return job_id
        else:
            cmd = 'sh %s > %s 2> %s' % (jobfile, stdout, stderr)
            if not args.silent:
                print("Local execution: %s" % cmd)

            subprocess.call(cmd, shell=True)
