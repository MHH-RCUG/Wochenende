import pytest
import run_Wochenende as we
import os
from os.path import isfile
from shutil import which, copytree
from subprocess import run, PIPE

@pytest.fixture()
def setup_tmpdir(request, tmpdir):
    """Fixture to set up a temporary testing directory"""
    os.chdir(request.config.invocation_dir)
    # ON UPDATE to python >= v3.8, use dirs_exist_ok=True for copytree instead of rmdir
    os.rmdir(tmpdir)
    copytree("..", tmpdir, symlinks=True)
    return tmpdir

@pytest.mark.parametrize("readType,longread", [('SE',False), ('PE',False), ('SE',True)])
def test_run(request, setup_tmpdir, readType, longread):
    """Runs the pipeline for a functional test."""
    os.chdir(setup_tmpdir)
    print(os.listdir())
    cmd = ['python3', 'run_Wochenende.py',
           '--aligner', 'minimap2' if longread else 'bwamem',
           '--readType', readType,
           '--metagenome', 'test',
           '--threads', '8',
           '--debug']
    if longread:
        cmd += ['--longread']
    cmd += ['test/data/reads_R1.fastq']

    proc = run(cmd, stdout=PIPE, stderr=PIPE)

    stdout = str(proc.stdout).replace("\\n", "\n\t")
    stderr = str(proc.stderr).replace("\\n", "\n\t")

    print(f'# Return Code: {str(proc.returncode)}')
    print(f'# stdout:\n\t{stdout}')
    print(f'# stderr:\n\t{stderr}')

    assert proc.returncode == 0
    for f in os.listdir('test/data'):
        assert os.path.getsize(f"test/data/{f}") > 0

    os.chdir(request.config.invocation_dir)

@pytest.mark.parametrize("path", [p for _, p in we.path_refseq_dict.items()])
def test_reference_paths(path):
    """Tests all reference paths for existence"""
    assert isfile(f"../{path}") or isfile(path)

@pytest.mark.parametrize("path", [p for varname, p in vars(we).items() if ('path' in varname) and not varname in "path_tmpdir,path_refseq_dict"])
def test_program_paths(path):
    """Tests all program / binary paths for existence"""
    assert which(path) != None or isfile(f"../{path}") or isfile(path)

@pytest.mark.parametrize("path", [p for varname, p in vars(we).items() if 'adapter' in varname])
def test_adapter_paths(path):
    """Tests all adapter paths for existence"""
    assert isfile(f"../{path}") or isfile(path)
