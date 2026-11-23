#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Cluster mode tests for bamtagmultiome (--cluster).

The cluster sub job commands are rebuilt from sys.argv, therefore these tests
invoke the CLI as a subprocess; running run_multiome_tagging in-process under
pytest would make the worker jobs re-invoke pytest itself.

The SLURM tests are optional and skipped when sbatch is not available:
    pytest tests/test_cluster.py -m slurm -v
"""
import glob
import os
import re
import shutil
import subprocess
import sys
import time
import unittest

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bam_test_helpers import write_misc_alt_bam, count_reads_by_name

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
READ_NAMES = ['main1', 'main2', 'alt1', 'alt2', 'alt3', 'unmapped1']
# The amount of seconds to wait for the cluster jobs to finish:
SLURM_TIMEOUT_S = int(os.environ.get('SCMO_TEST_SLURM_TIMEOUT', 600))


def find_cli():
    """Locate the bamtagmultiome CLI script"""
    for name in ('bamtagmultiome.py', 'bamtagmultiome'):
        path = shutil.which(name)
        if path is not None:
            return path
    return os.path.join(REPO_ROOT, 'singlecellmultiomics', 'universalBamTagger', 'bamtagmultiome.py')


CLI_PATH = find_cli()


def cli_env():
    """Environment for the CLI subprocess; cluster job scripts re-invoke the CLI via its shebang,
    make sure it resolves to this interpreter and this package"""
    env = dict(os.environ)
    env['PATH'] = os.path.dirname(sys.executable) + os.pathsep + env.get('PATH', '')
    env['PYTHONPATH'] = REPO_ROOT + os.pathsep + env.get('PYTHONPATH', '')
    return env


def run_cli(args, timeout=900):
    return subprocess.run([sys.executable, CLI_PATH] + args, capture_output=True, text=True,
                          timeout=timeout, env=cli_env())


class TestClusterLocal(unittest.TestCase):
    """ Round trip through the cluster job submission code using the local scheduler,
    runs on any machine without SLURM """

    write_path = './data/write_test_cluster_local.bam'
    cluster_dir = './data/cluster_test_local'

    @classmethod
    def setUpClass(cls):
        cls.bam_path = write_misc_alt_bam('./data/cluster_local_test.bam')

    @classmethod
    def tearDownClass(cls):
        for path in [cls.bam_path, f'{cls.bam_path}.bai', cls.write_path, f'{cls.write_path}.bai',
                     cls.write_path.replace('.bam', '.status.txt')]:
            if os.path.exists(path):
                os.remove(path)
        shutil.rmtree(cls.cluster_dir, ignore_errors=True)

    def test_cluster_local_roundtrip(self):
        result = run_cli([self.bam_path, '-method', 'qflag', '--cluster', '-sched', 'local',
                          '-clusterdir', self.cluster_dir, '-o', self.write_path])
        self.assertEqual(result.returncode, 0, msg=result.stderr[-2000:])
        # The local scheduler runs all jobs synchronously, the merge job writes:
        with open(self.write_path.replace('.bam', '.status.txt')) as status_file:
            self.assertIn('All done', status_file.read())
        # Every read must occur exactly once in the merged file, unmapped reads included:
        self.assertEqual(count_reads_by_name(self.write_path, READ_NAMES),
                         {name: 1 for name in READ_NAMES})
        # The temporary contig bams and their status markers must be cleaned up:
        temp_dir = os.path.dirname(os.path.abspath(self.write_path))
        self.assertEqual(glob.glob(os.path.join(temp_dir, 'SCMO_*')), [])


@pytest.mark.slurm
@unittest.skipIf(shutil.which('sbatch') is None, 'SLURM not available')
class TestClusterSlurm(unittest.TestCase):
    """ Round trip through real SLURM submission, opt-in with: -m slurm """

    @classmethod
    def setUpClass(cls):
        cls.bam_path = write_misc_alt_bam('./data/cluster_slurm_test.bam')
        cls.bam_no_unmapped = write_misc_alt_bam('./data/cluster_slurm_test_no_unmapped.bam',
                                                 include_unmapped=False)

    @classmethod
    def tearDownClass(cls):
        for path in [cls.bam_path, f'{cls.bam_path}.bai',
                     cls.bam_no_unmapped, f'{cls.bam_no_unmapped}.bai']:
            if os.path.exists(path):
                os.remove(path)

    def submit_and_wait(self, bam_path, write_path, cluster_dir):
        result = run_cli([bam_path, '-method', 'qflag', '--cluster', '-sched', 'slurm',
                          '-clusterdir', cluster_dir, '-o', write_path])
        self.assertEqual(result.returncode, 0, msg=result.stderr[-2000:])
        # Best effort cleanup of submitted jobs if the test fails while they are pending:
        for job_id in re.findall(r'submitted with job id:\s*(\d+)', result.stdout):
            self.addCleanup(lambda jid=job_id: subprocess.run(
                ['scancel', jid], capture_output=True))
        status_path = write_path.replace('.bam', '.status.txt')
        deadline = time.time() + SLURM_TIMEOUT_S
        while time.time() < deadline:
            if os.path.exists(status_path):
                with open(status_path) as status_file:
                    if 'All done' in status_file.read():
                        return result
            time.sleep(10)
        cluster_errors = '\n'.join(
            f'{path}: {open(path).read()[-500:]}' for path in glob.glob(f'{cluster_dir}/*.stderr'))
        self.fail(f'Cluster jobs did not finish within {SLURM_TIMEOUT_S}s. '
                  f'Status file content: {open(status_path).read() if os.path.exists(status_path) else None}. '
                  f'Cluster stderr: {cluster_errors}')

    def cleanup_output(self, write_path, cluster_dir):
        for path in [write_path, f'{write_path}.bai', write_path.replace('.bam', '.status.txt')]:
            if os.path.exists(path):
                os.remove(path)
        shutil.rmtree(cluster_dir, ignore_errors=True)
        temp_dir = os.path.dirname(os.path.abspath(write_path))
        for path in glob.glob(os.path.join(temp_dir, 'SCMO_*')):
            if os.path.isfile(path):
                os.remove(path)
            else:
                shutil.rmtree(path, ignore_errors=True)

    def test_cluster_slurm_roundtrip(self):
        write_path = './data/write_test_cluster_slurm.bam'
        cluster_dir = './data/cluster_test_slurm'
        self.addCleanup(self.cleanup_output, write_path, cluster_dir)
        self.submit_and_wait(self.bam_path, write_path, cluster_dir)
        # Every read must occur exactly once in the merged file, unmapped reads included:
        self.assertEqual(count_reads_by_name(write_path, READ_NAMES),
                         {name: 1 for name in READ_NAMES})
        # The per-contig temporary files and status markers are removed by the merge job:
        temp_dir = os.path.dirname(os.path.abspath(write_path))
        self.assertEqual(glob.glob(os.path.join(temp_dir, 'SCMO_*')), [])

    def test_cluster_slurm_no_unmapped_job_without_unmapped_reads(self):
        write_path = './data/write_test_cluster_slurm_no_unmapped.bam'
        cluster_dir = './data/cluster_test_slurm_no_unmapped'
        self.addCleanup(self.cleanup_output, write_path, cluster_dir)
        self.submit_and_wait(self.bam_no_unmapped, write_path, cluster_dir)
        # No dedicated '*' job should have been created for an input without unmapped reads:
        for script_path in glob.glob(os.path.join(cluster_dir, '*.sh')):
            with open(script_path) as script_file:
                self.assertNotIn("-contig '*'", script_file.read(),
                                 msg=f'unmapped job submitted for input without unmapped reads')
        counts = count_reads_by_name(write_path, READ_NAMES)
        self.assertEqual(counts['unmapped1'], 0)
        self.assertEqual({name: counts[name] for name in ['main1', 'main2', 'alt1', 'alt2', 'alt3']},
                         {name: 1 for name in ['main1', 'main2', 'alt1', 'alt2', 'alt3']})


if __name__ == '__main__':
    unittest.main()
