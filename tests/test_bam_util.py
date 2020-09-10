#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import unittest
import singlecellmultiomics.molecule
import singlecellmultiomics.fragment
import pysam
import pysamiterators.iterators
from singlecellmultiomics.bamProcessing import sorted_bam_file,write_program_tag,verify_and_fix_bam
from singlecellmultiomics.bamProcessing.bamExtractSamples import extract_samples
import os
import sys
from shutil import copyfile,rmtree
from singlecellmultiomics.bamProcessing import get_contigs_with_reads
from singlecellmultiomics.utils.sequtils import pick_best_base_call

class TestFunctions(unittest.TestCase):

    def test_get_contigs_with_reads_1(self):
        cwr = list(get_contigs_with_reads('./data/mini_nla_test.bam'))
        self.assertEqual(len(cwr), 1)
        self.assertIn('chr1', cwr)

    def test_get_contigs_with_reads_2(self):
        cwr = list(get_contigs_with_reads('./data/chic_test_region.bam'))
        self.assertEqual(len(cwr), 1)
        self.assertIn('8', cwr)

class TestSorted(unittest.TestCase):

    def test_verify_and_fix_bam_autoindex(self):
        file_path_without_index = './data/temp_without_index.bam'
        copyfile('./data/mini_nla_test.bam',file_path_without_index)

        try:
            os.remove(file_path_without_index+'.bai')
        except Exception as e:
            pass

        verify_and_fix_bam(file_path_without_index)

        self.assertTrue( os.path.exists(file_path_without_index+'.bai') )

        with pysam.AlignmentFile(file_path_without_index) as f:
            i =0
            # Test if the file has reads.
            for read in f.fetch(contig='chr1'):
                if read.is_read1:
                    i+=1
            self.assertEqual(i, 293)

        try:
            os.remove(file_path_without_index)
        except Exception as e:
            pass

        try:
            os.remove(file_path_without_index+'.bai')
        except Exception as e:
            pass


    def test_write_to_sorted(self):
        write_path = './data/write_test.bam'
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            with sorted_bam_file(write_path, origin_bam=f) as out:
                for molecule in singlecellmultiomics.molecule.MoleculeIterator(
                    alignments=f,
                    molecule_class=singlecellmultiomics.molecule.NlaIIIMolecule,
                    fragment_class=singlecellmultiomics.fragment.NlaIIIFragment,
                    fragment_class_args={'umi_hamming_distance':0},
                    pooling_method=0,
                    yield_invalid=True
                ):
                    molecule.write_pysam(out)

        self.assertTrue(os.path.exists(write_path))
        try:
            os.remove(write_path)
            os.remove(write_path+'.bai')
        except Exception as e:
            pass

    def test_write_to_sorted_custom_compression(self):
        write_path = './data/write_test.bam'
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            with sorted_bam_file(write_path, origin_bam=f,fast_compression=True) as out:
                for molecule in singlecellmultiomics.molecule.MoleculeIterator(
                    alignments=f,
                    molecule_class=singlecellmultiomics.molecule.NlaIIIMolecule,
                    fragment_class=singlecellmultiomics.fragment.NlaIIIFragment,
                    fragment_class_args={'umi_hamming_distance':0},
                    pooling_method=0,
                    yield_invalid=True
                ):
                    molecule.write_pysam(out)

        self.assertTrue(os.path.exists(write_path))
        try:
            os.remove(write_path)
            os.remove(write_path+'.bai')
        except Exception as e:
            pass

    def test_write_to_sorted_non_existing_folder(self):
        write_folder =  './data/non_yet_existing_folder/'
        write_path = write_folder + 'write_test.bam'
        if os.path.exists(write_path):
            os.remove(write_path)

        rmtree(write_folder, ignore_errors=True)

        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:
            with sorted_bam_file(write_path, origin_bam=f) as out:
                for molecule in singlecellmultiomics.molecule.MoleculeIterator(
                    alignments=f,
                    molecule_class=singlecellmultiomics.molecule.NlaIIIMolecule,
                    fragment_class=singlecellmultiomics.fragment.NlaIIIFragment,
                    fragment_class_args={'umi_hamming_distance':0},
                    pooling_method=0,
                    yield_invalid=True
                ):
                    molecule.write_pysam(out)

        self.assertTrue(os.path.exists(write_path))

        with pysam.AlignmentFile(write_path) as f:
            i =0
            # Test if the file has reads.
            for read in f:
                if read.is_read1:
                    i+=1
            self.assertEqual(i, 293)

        try:
            os.remove(write_path)
            os.remove(write_path+'.bai')
        except Exception as e:
            pass

        rmtree(write_folder, ignore_errors=True)

    def test_write_to_read_grouped_sorted(self):
        write_path = './data/write_test_rg.bam'
        read_groups = set()
        with pysam.AlignmentFile('./data/mini_nla_test.bam') as f:

            input_header = f.header.as_dict()
            write_program_tag(input_header,
                program_name='test_bam_util_test1',
                command_line = " ".join(sys.argv),
                version = singlecellmultiomics.__version__,
                description = f'a description'
                )

            write_program_tag(input_header,
                program_name='test_bam_util_test2',
                command_line = " ".join(sys.argv),
                version = singlecellmultiomics.__version__,
                description = f'a description'
                )
            #print([x for x in input_header['PG'] if not 'bwa mem' in x.get('CL','')])
            with sorted_bam_file(write_path, header=input_header,read_groups=read_groups) as out:
                for molecule in singlecellmultiomics.molecule.MoleculeIterator(
                    alignments=f,
                    molecule_class=singlecellmultiomics.molecule.NlaIIIMolecule,
                    fragment_class=singlecellmultiomics.fragment.NlaIIIFragment,
                    fragment_class_args={'umi_hamming_distance':0},
                    pooling_method=0,
                    yield_invalid=True
                ):
                    molecule.write_pysam(out)
                    for frag in molecule:
                        read_groups.add( frag.get_read_group() )

        self.assertTrue(os.path.exists(write_path))

        # Now test if the program tag is there...
        with pysam.AlignmentFile(write_path) as f:
            self.assertTrue( 1==len([x for x in f.header['PG'] if 'test_bam_util_test1' in x.get('PN','')]) )
            self.assertTrue( 1==len([x for x in f.header['PG'] if 'test_bam_util_test2' in x.get('PN','')]) )
            i =0

            # Test if the file has reads.
            for read in f:
                if read.is_read1:
                    i+=1
            self.assertEqual(i, 293)

        try:
            os.remove(write_path)
        except Exception as e:
            pass

        try:
            os.remove(write_path+'.bai')
        except Exception as e:
            pass


    def  test_sample_extraction(self):

        output_path= './data/write_test_extract.bam'
        final_output_paths= ['./data/write_test_extract0.bam', './data/write_test_extract1.bam']

        capture_samples = {'0':['APKS1P25-NLAP2L1_30','APKS1P25-NLAP2L1_31'],
        # 4 reads expected
         '1':['APKS1P25-NLAP2L2_86']}
         #samtools view ./data/mini_nla_test.bam | grep 'APKS1P25-NLAP2L2_86' | wc -l -> 8
        extract_samples( './data/mini_nla_test.bam', output_path, capture_samples )

        with pysam.AlignmentFile(final_output_paths[0]) as f:
            for i,_ in enumerate(f):
                pass
            self.assertEqual(i, 4-1)

        with pysam.AlignmentFile(final_output_paths[1]) as f:
            for i,_ in enumerate(f):
                pass
            self.assertEqual(i, 8-1)


        for p in final_output_paths:
            try:
                os.remove(p)
                os.remove(f'{p}.bai')
            except Exception as e:
                pass

class TestBaseCalling(unittest.TestCase):

    def test_pick_best(self):

        self.assertTrue( pick_best_base_call( ('A',32), ('C',22) ) == ('A', 32) )
        self.assertTrue( pick_best_base_call( ('C',22), ('A',32) ) == ('A', 32))
        self.assertTrue( pick_best_base_call( ('C',32), ('A',32) ) == ('N',0) )





if __name__ == '__main__':
    unittest.main()
