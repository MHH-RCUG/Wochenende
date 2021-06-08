import pytest
import subprocess
import os
# import re
from itertools import combinations, product, chain


class TestFunctional:
    __test__ = False
    """Tests for the overall functionality of the Wochenende pipeline.

    Attributes
    ----------
    _p_fixed_arguments : dict [str, str]
        Fixed arguments that need to be set like this for testing. Value is empty if the
        argument just needs to exist
    _p_options : dict [str, list of str]
        Optional arguments that take an input. Each list gives the inputs for testing.
    _p_flags : list of str
        Optional arguments that do not require an additional input.

    Notes
    -----
    To run only a few combinations, just comment the irrelevant ones.
    This test class relies on an existing slurm installation. It might be too much to run
    this due to combinatorial explosion. For a simple test of installation, one can rely
    on `test_installation.py`. 
    """

    # TODO: PE support
    # TODO: Test partially successful runs, continuation and --force-restart
    # TODO: Collect test run information
    # TODO: Sum up test run information

    _p_fixed_arguments = [
        '--metagenome', 'test',
        '--debug'
    ]
    _p_options = {
        '--aligner': ['bwamem', 'minimap2', 'ngmlr'],
        #'--readType': ['SE', 'PE'],
        #'--threads': ['8', '4', '8', '16', '32', '56', '80'],
        #'--remove_mismatching': ['1', '2', '3', '4', '5'],
    }
    _p_flags = [
        #'--fastp',
        #'--nextera',
        #'--trim_galore',
        #'--longread',
        #'--no_duplicate_removal',
        #'--no_prinseq',
        #'--no_fastqc',
        #'--no_abra',
        #'--mq30',
        #'--force_restart'
    ]

    def _gen_fixed_arguments(self):
        """Returns the fixed arguments. Mainly for consistency."""
        return self._p_fixed_arguments

    def _gen_options(self):
        """Generates all combinations of p_options

        Please don't ask how. It works. Trust me. Sorry.

        """
        ds = list(chain(*[list(map(dict, combinations(self._p_options.items(), i)))
                          for i in range(len(self._p_options) + 1)]))
        return chain(*[[list(chain(*x)) for x in product(*[[[k, v] for v in d[k]] for k in d.keys()])] for d in ds])

    def _gen_flags(self):
        """Generates all combinations of p_flags"""
        for length in range(len(self._p_flags) + 1):
            for combination in combinations(self._p_flags, length):
                yield list(combination)

    def _gen_pipeline_arguments(self):
        """ Generates pipeline argument strings

        Yield
        -----
        Pipeline argument combinations for a full feature test

        """
        print("\n# Arguments to test")
        fixed = list(self._gen_fixed_arguments())
        print(fixed)

        print("# Options to test")
        options = list(self._gen_options())
        print(options)

        print("# Flags to test")
        flags = list(self._gen_flags())
        # print(flags)

        return [fixed + opts + flgs for opts, flgs in product(options, flags)]

    def _get_refs_and_reads_paths(self):
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

        cwd = os.getcwd()
        assert cwd.endswith('test') or cwd.endswith('test/')

        datadir = os.path.join(cwd, 'data')
        datalist = os.listdir(datadir)

        refs = [os.path.join(datadir, f) for f in datalist if f.endswith('.fa')]
        reads = [os.path.join(datadir, f) for f in datalist if f.endswith('_R1.fastq')]

        return refs, reads

    def _get_fixed_files_paths(self):
        """Returns singleton list of list of absolute paths of necessary files."""

        cwd = os.getcwd()
        assert cwd.endswith('test') or cwd.endswith('test/')

        main_we_dir = os.path.join(cwd, '..')

        files = [os.path.join(main_we_dir, f) for f in os.listdir(main_we_dir)]

        return [files]

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

    def pytest_generate_tests(self, metafunc):
        """Main test generator"""

        refs, reads = self._get_refs_and_reads_paths()
        files = self._get_fixed_files_paths()

        if 'files' in metafunc.fixturenames:
            metafunc.parametrize('files', files)

        if 'reference' in metafunc.fixturenames:
            metafunc.parametrize('reference', refs)

        if 'read' in metafunc.fixturenames:
            metafunc.parametrize('read', reads)

        if 'pipeline_arguments' in metafunc.fixturenames:
            metafunc.parametrize(argnames='pipeline_arguments',
                                 argvalues=self._gen_pipeline_arguments())

    @pytest.mark.slow
    def test_pipeline_call(self, setup_tmpdir, pipeline_arguments):
        """Main pipeline test using subprocesses"""

        # change to temporary directory
        os.chdir(setup_tmpdir)
        print(f'\n# Running Tests in: {os.getcwd()}')
        print(f'# Seeing files: {os.listdir(setup_tmpdir)}')

        # Create slurm job script
        slurm_template = os.path.join('test', 'slurm-template.sh')
        with open(slurm_template, 'r') as f:
            slurm_file = f.read()

        slurm_cmd = ['python3', 'run_Wochenende.py'] + list(pipeline_arguments) + \
                    ['$fastq']
        slurm_file = slurm_file + ' '.join(slurm_cmd) + '\n wait'

        sbatch_test_filename = 'run_Wochenende_SLURM_test.sh'
        with open(sbatch_test_filename, 'w') as f:
            f.write(slurm_file)
            f.flush()

        # print(slurm_file)

        # Start slurm subprocess
        cmd = ['srun', sbatch_test_filename, 'reads_R1.fastq']
        # print(' '.join(cmd))

        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            # shell=True # Seems to be necessary (?)
        )

        stdout = str(proc.stdout).replace("\\n", "\n\t")
        stderr = str(proc.stderr).replace("\\n", "\n\t")

        print(f'# Return Code: {str(proc.returncode)}')
        print(f'# stdout:\n\t{stdout}')
        print(f'# stderr:\n\t{stderr}')

        assert proc.returncode == 0
