import subprocess
import os
from importlib.util import spec_from_loader, module_from_spec
from importlib.machinery import SourceFileLoader
import pytest
from argparse import Namespace
from tempfile import TemporaryDirectory

ROOT = os.path.dirname(os.path.dirname(__file__))
peakhood_path = os.path.join(ROOT, "bin", "peakhood")

spec = spec_from_loader("peakhood",
                        SourceFileLoader("peakhood", peakhood_path))
peakhood = module_from_spec(spec)
spec.loader.exec_module(peakhood)


@pytest.fixture
def temp_batch():
    tmpdir = TemporaryDirectory()
    yield tmpdir.name


@pytest.fixture
def default_batch(temp_batch):
    args = Namespace(bam_pp_mode=1,
                     eib_ratio_mode=1,
                     eib_width=20,
                     f2_mode=1,
                     in_2bit="/home/rabsch/Documents/denbi_mount/rabsch/DATA/hg38.2bit",
                     in_folder='tests/testdata/batch',
                     in_gtf='/home/rabsch/Documents/denbi_mount/rabsch/DATA/Homo_sapiens.GRCh38.104.gtf.gz',
                     isr_ext_mode=1,
                     isr_max_reg_len=10,
                     isrn_prefilter=False,
                     list_f1_filter=False,
                     list_f2_filter=False,
                     list_tbt_filter_ids=False,
                     max_len=250,
                     merge_ext=0,
                     merge_gtf=False,
                     merge_mode=1,
                     min_ei_ratio=4,
                     min_eib_ratio=4,
                     min_eol=0.9,
                     min_exbs_isr_c=2,
                     min_n_tis_sites=3,
                     new_ids=True,
                     no_eibr_filter=False,
                     no_eibr_wt_filter=False,
                     no_eir_wt_filter=False,
                     no_f1_filter=False,
                     no_tbt_filter=False,
                     no_tis_filter=False,
                     out_folder=temp_batch,
                     pre_merge=True,
                     read_pos_mode=1,
                     report=True,
                     rev_filter=False,
                     rnaseq_bam=False,
                     rnaseq_bam_rev=False,
                     sc_thr=None,
                     seq_ext=0,
                     seq_ext_mode=1,
                     which="batch")
    return args


def _compare_folders(expected_folder, test_folder):
    expected = os.listdir(expected_folder)
    test = os.listdir(test_folder)
    for file in expected:
        expected_file = os.path.join(expected_folder, file)
        test_file = os.path.join(test_folder, file)
        if ".log" in file or "settings" in file or not ".bed" in file or not ".fa" in file:
            continue
        assert os.path.exists(test_file)
        if os.path.isdir(expected_file):
            _compare_folders(expected_file, test_file)
        elif os.path.isfile(expected_file):
            with open(expected_file, "r") as exp, open(test_file, "r") as tes:
                test_set = set()
                for exp_line in exp:
                    test_set.add(exp_line)
                for test_line in tes:
                    if test_line not in test_set:
                        assert test_line in test_set


def test_foo(default_batch):
    peakhood.main_batch(default_batch)
    expeted = os.path.join(os.path.dirname(__file__), "expected", "batch_expected")
    _compare_folders(os.path.join(expeted), default_batch.out_folder)

    x = 0
    p = 0