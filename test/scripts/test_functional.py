#!/usr/bin/python3

"""
Wochenende - Functional Tests
=============================

These tests are for checking the overall functionality of the Wochenende 
pipeline.

For this we generate tests as follows
0. Set up the parameters for the pipeline run
1. Create a temprorary folder for the new pipeline run
2. Start the actual run using the predefined parameters
3. Collect information about the run
4. Sum up the collected information

Not done yet
------------
* Full PE support
* Tests for reporting
* Tests for plots
* Test for only partially successful runs
"""

import pytest
import subprocess
import os
from itertools import combinations, product

class TestFunctional():
    p_arguments = [ # fixed arguments that need to be set like this
        '--metagenome test', # TODO(B1T0): change run_Wochenende.py accordingly
        '--readType SE', # TODO(B1T0): add PE support
        '--debug',
    ]
    p_options = { # arguments that take an input
        '--aligner': ['bwamem'  ],#  , 'minimap2', 'ngmlr'],
        '--threads': [1, 2  ],#  , 4, 8, 16, 32, 56, 80],
        '--remove_mismatching': [1  ]#  ,2,3,4,5],
    }
    p_flags = [ # optional arguments that do not require an additional input
        '--fastp',
        '--nextera',
        # '--trim_galore',
        # '--longread',
        # '--no_duplicate_removal',
        # '--no_prinseq',
        # '--no_fastqc',
        # '--no_abra',
        # '--mq30',
        # '--force_restart' # TODO(B1T0): implement tests that make use of this
    ]


    def gen_arguments(self):
        """Returns the joined p_arguments"""
        return list(self.p_arguments)


    def gen_options(self):
        """Generates all combinations of p_options"""
        for length in range(1, len(self.p_options)+1):
            for combination in combinations(self.p_options, length):
                for option, arguments in combination:
                    for argument in arguments:
                        yield f'{option} {argument}'


    def gen_flags(self):
        """Generates all combinations of p_flags"""
        for length in range(1, len(self.p_flags)+1):
            for combination in combinations(self.p_flags, length):
                yield list(combination)


    def gen_pipeline_arguments(self):
        """
        Generates pipeline argument strings

        Yield:
            Pipeline argument combinations for a full feature test
        """
        print("#### Arguments ####")
        print(self.gen_arguments())
        print("#### Options ####")
        print(self.gen_options())
        for o in self.gen_options():
            print(o)
        print("#### Flags ####")
        print(self.gen_flags())
        for f in self.gen_flags():
            print(f)
        return list(product(self.gen_arguments(), self.gen_options(), self.gen_flags()))
    

    def get_test_data(self):
        """
        Returns absolute paths of references (.fa files) and reads (.fastq 
        files) from /path/to/wochenende/test/data.

        Note:
            For testing, we provide explicit reference files, no reference from
            one of the refseq_dict in the script's configuration section.
        """
        # TODO(B1T0): add PE support

        cwd = os.getcwd()
        assert cwd.endswith('test') or cwd.endswith('test/')

        datadir = os.path.join(cwd, 'data')

        datalist = os.listdir(datadir)
        refs = [os.path.join(datadir, f) for f in datalist if f.endswith('.fa')]
        reads = [os.path.join(datadir, f) for f in datalist 
                 if f.endswith('_R1.fastq')]

        return refs, reads


    def test_testfiles(self):
        """
        Checks if /path/to/wochenende/test/data contains at least the necessary
        two files, as it should contain at least one reference (.fa) and one 
        reads file (.fastq).

        Note:
            For testing, we provide explicit reference files, no reference from
            one of the refseq_dict in the script's configuration section.
        """

        refs, reads = self.get_test_data()

        assert refs
        assert reads

        print('# Found references:' + ', '.join(list(refs)))
        print('# Found reads:' + ', '.join(list(reads)))


    def pytest_generate_tests(self, metafunc):
        # TODO(B1T0): add PE support

        refs, reads = self.get_test_data()

        if 'reference' in metafunc.fixturenames:
            metafunc.parametrize('reference', refs)

        if 'read' in metafunc.fixturenames:
            metafunc.parametrize('read', reads)

        if 'pipeline_arguments' in metafunc.fixturenames:
            metafunc.parametrize(argnames = 'pipeline_arguments', argvalues = self.gen_pipeline_arguments())


    @pytest.fixture()
    def setup_tmpdir(self, tmpdir, reference, read):
        """
        Sets up the temporary working directory by creating it and symlinking 
        the necessary files.

        Args:
            tmpdir (py.path.local): fixture that creates a temporary directory
            ref (str): absolute path to reference file
            read (str): absolute path to read file
        """
        # TODO(B1T0): add PE support
        # TODO(B1T0): modify for reporting
        # TODO(B1T0): modify for plots

        os.symlink(reference, os.path.join(tmpdir, 'ref.fa'))
        os.symlink(read, os.path.join(tmpdir, 'reads_R1.fastq'))
        # if os.path.exists(read.replace('_R1', '_R2')):
        #     os.symlink(read.replace('_R1', '_R2'), 
        #                 os.path.join(tmpdir, 'reads_R1.fastq'))

        main_we_path = os.path.join(tmpdir, "..", "..")
        # assumes we have main_we_path/test/tmpdir

        files_needed = [
            'run_Wochenende.py'
        ]
        for f in files_needed:
            os.symlink(os.path.join(main_we_path, f), os.path.join(tmpdir, f))

        return tmpdir


    def test_pipeline(self, setup_tmpdir, pipeline_arguments):
        """Main pipeline test"""
        os.chdir(setup_tmpdir)

        cmd = ['python3', 'run_Wochenende.py'] + list(pipeline_arguments) + ['reads_R1.fastq']

        print(f'# Current working directory: {os.getcwd()}')
        print(f'The following files exist here: {os.listdir(setup_tmpdir)}')
        print(f'# The following line will be executed {cmd}')
        # assert 0

        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        print('Return Code: ' + str(proc.returncode))
        print('stdout:\n' + str(proc.stdout))
        print('stderr:\n' + str(proc.stderr))

        assert proc.returncode
        assert 0




""" old test script (bash)
conda activate wochenende

python3 run_Wochenende.py --metagenome testdb --testWochenende ../testdb/reads_R1.fastq --force_restart
"""



"""get_wochenende.sh
path_we=/mnt/ngsnfs/tools/dev/Wochenende

cp $path_we/README* .
cp $path_we/*.sh .
cp $path_we/*.py .
cp -R $path_we/plots/ .
cp -R $path_we/reporting/ .
"""