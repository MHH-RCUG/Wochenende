#!/usr/bin/python3

"""
Wochenende - Functional Tests
=============================

These tests are for checking the overall functionality of the Wochenende pipeline.

For this we generate tests as follows
0. Set up the parameters for the pipeline run
1. Create a temporary folder for the new pipeline run
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
from itertools import combinations, product, chain


class TestFunctional:
    p_arguments = [  # fixed arguments that need to be set like this for testing
        '--metagenome', 'test',  # TODO(B1T0): change run_Wochenende.py accordingly
        '--debug',
    ]
    p_options = {  # optional arguments that take an input
        '--aligner': ['bwamem'],  # , 'minimap2', 'ngmlr'],
        '--readType': ['SE', 'PE'], # TODO: still need to distinguish for PE
        '--threads': ['1', '2'],  # , 4, 8, 16, 32, 56, 80],
        # '--remove_mismatching': ['1']  # ,2,3,4,5],
    }
    p_flags = [  # optional arguments that do not require an additional input
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
        """Generates all combinations of p_options

        Please don't ask how. It works. Trust me.

        """
        ds = list(chain(*[list(map(dict, combinations(self.p_options.items(), i)))
                          for i in range(len(self.p_options) + 1)]))
        return chain(*[[list(chain(*x)) for x in product(*[[[k, v] for v in d[k]] for k in d.keys()])] for d in ds])

        # for length in range(2, len(self.p_options) + 1):
        #     args = self.p_options.keys()
        #     iss = [self.p_options[arg] for arg in args]
        #     combos = product(*iss)
        #     for combination in product(*iss):
        #         # combinations(self.p_options.keys(), length):
        #         yield list(chain())
        #         precombine = product(*[[[o, a] for a in self.p_options[o]] for o in combination])
        #         yield list(chain(*precombine))

    def gen_flags(self):
        """Generates all combinations of p_flags"""
        for length in range(len(self.p_flags) + 1):
            for combination in combinations(self.p_flags, length):
                yield list(combination)

    def gen_pipeline_arguments(self):
        """ Generates pipeline argument strings

        Yield
        -----
        Pipeline argument combinations for a full feature test

        """
        print("#### Arguments ####")
        print(self.gen_arguments())
        print("#### Options ####")
        print(self.gen_options())
        # for o in self.gen_options():
        #     print(o)
        print("#### Flags ####")
        print(self.gen_flags())
        # for f in self.gen_flags():
        #     print(f)
        return [self.gen_arguments() + a + b for a,b in product(self.gen_options(), self.gen_flags())]

    def get_test_data(self):
        """Returns absolute paths of reference (.fa) and read (.fastq) test files.

        Returns
        -------
        refs : list of os.path
            absolute paths to reference file(s)
        reads : list of os.path
            absolute paths to read file(s)

        Note
        ----
        For testing, we provide explicit reference files. No reference from one of the
        refseq_dict in the script's configuration section is used. All test files need to
        be stored in /path/to/wochenende/test/data.

        """
        # TODO(B1T0): add PE support

        # print(f'### {os.getcwd()} ###')
        cwd = os.getcwd()
        # print(f'### Working in {cwd} ###')
        assert cwd.endswith('test') or cwd.endswith('test/')

        datadir = os.path.join(cwd, 'data')

        datalist = os.listdir(datadir)
        refs = [os.path.join(datadir, f) for f in datalist if f.endswith('.fa')]
        reads = [os.path.join(datadir, f) for f in datalist
                 if f.endswith('_R1.fastq')]

        return refs, reads

    def get_files(self):
        """Returns absolute path of the run_Wochenende.py script.

        Returns
        -------
        script : os.path
            absolute paths to script file

        Note
        ----
        All test files need to be stored in /path/to/wochenende/test/data.

        """

        # print(f'### {os.getcwd()} ###')
        cwd = os.getcwd()
        # print(f'### Working in {cwd} ###')
        assert cwd.endswith('test') or cwd.endswith('test/')
        needed_files = [
            'data/TruSeq3-PE.fa',
            '../run_Wochenende.py',
            '../dependencies'
        ]
        return [os.path.join(cwd, f) for f in needed_files]

    def test_testfiles(self):
        """Checks existence of test files.

        Checks if /path/to/wochenende/test/data contains at least the necessary two
        files, as it should contain at least one reference (.fa) and one reads file
        (.fastq).

        Notes
        -----
        For testing, we provide explicit reference files, no reference from one of the
        refseq_dict in the script's configuration section.

        """

        refs, reads = self.get_test_data()

        assert refs
        assert reads

        print('# Found references:' + ', '.join(list(refs)))
        print('# Found reads:' + ', '.join(list(reads)))

    def pytest_generate_tests(self, metafunc):
        """Main test generator"""
        # TODO: documentation
        # TODO: add PE support

        refs, reads = self.get_test_data()
        files = self.get_files()

        if 'files' in metafunc.fixturenames:
            metafunc.parametrize('files', [files])

        if 'reference' in metafunc.fixturenames:
            metafunc.parametrize('reference', refs)

        if 'read' in metafunc.fixturenames:
            metafunc.parametrize('read', reads)

        if 'pipeline_arguments' in metafunc.fixturenames:
            metafunc.parametrize(argnames='pipeline_arguments',
                                 argvalues=self.gen_pipeline_arguments())

    @pytest.fixture()
    def setup_tmpdir(self, files, read, reference, tmpdir):
        """Sets up the temporary working directory.

        Creates the directory and symlinks the necessary files into it.

        Parameters
        ----------
        files : list of str
            absolute paths to necessary files
        read : str
            absolute path to read file
        reference : str
            absolute path to reference file
        tmpdir : py.path.local
            fixture that creates a temporary directory

        Returns
        -------
        tmpdir : py.path.local
            path to the temporary directory

        """
        # TODO: add PE support
        # TODO: modify for reporting
        # TODO: modify for plots

        os.symlink(reference, os.path.join(tmpdir, 'ref.fa'))
        os.symlink(read, os.path.join(tmpdir, 'reads_R1.fastq'))
        if os.path.exists(read.replace('_R1', '_R2')):
            os.symlink(read.replace('_R1', '_R2'), os.path.join(tmpdir, 'reads_R2.fastq'))

        for f in files:
            # print(f.split('/'))
            os.symlink(f, os.path.join(tmpdir, f.split('/')[-1]))

        # main_we_path = os.path.join(os.getcwd(), "..", "..")
        # print(f'##### {main_we_path} #####')

        return tmpdir

    def test_pipeline(self, setup_tmpdir, pipeline_arguments):
        """Main pipeline test

        Parameters
        ----------
        setup_tmpdir
        pipeline_arguments

        Returns
        -------

        """
        # TODO: documentation

        os.chdir(setup_tmpdir)
        print(f'# Current working directory: {os.getcwd()}')
        # print(f'The following files exist here: {os.listdir(setup_tmpdir)}')

        cmd = ['python3', 'run_Wochenende.py'] + list(pipeline_arguments) + \
              ['reads_R1.fastq']
        # cmd = ['conda', 'info']
        print(f'# The following line will be executed\n {" ".join(cmd)}')
        print(f'# The following line will be executed\n {cmd}')

        # proc = subprocess.run(['ls','-lrtah'], stdout=subprocess.PIPE)
        #
        # print('Return Code: ' + str(proc.returncode))
        # print('stdout:\n' + str(proc.stdout))
        # print('stderr:\n' + str(proc.stderr))

        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            # shell=True # Seems to be necessary (?)
        )

        stdout = str(proc.stdout).replace("\\n", "\n")
        stderr = str(proc.stderr).replace("\\n", "\n")

        print(f'Return Code: {str(proc.returncode)}')
        print(f'stdout: {stdout}')
        print(f'stderr: {stderr}')

        assert proc.returncode == 0
