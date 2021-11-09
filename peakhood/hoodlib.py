#!/usr/bin/env python3

from distutils.spawn import find_executable
import matplotlib.pyplot as plt
# import plotly.express as px
import seaborn as sns
import pandas as pd
import numpy as np
import subprocess
import statistics
import random
import math
import gzip
import uuid
import sys
import re
import os


"""

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~ OPEN FOR BUSINESS ~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


~~~~~~~~~~~~~
Run doctests
~~~~~~~~~~~~~

python3 -m doctest hoodlib.py


"""


################################################################################

"""
Define additional chromosome names to be supported.

From Drosophila melanogaster:
chr2L
chr2R
chr3L
chr3R
2L
2R
3L
3R

"""

add_chr_names_dic = {
    "chr2L" : 1,
    "chr2R" : 1,
    "chr3L" : 1,
    "chr3R" : 1,
    "2L" : 1,
    "2R" : 1,
    "3L" : 1,
    "3R" : 1
}


################################################################################

def is_tool(name):
    """Check whether tool "name" is in PATH."""
    return find_executable(name) is not None


################################################################################

def dir_get_files(file_dir,
                  file_ending=False,
                  check=True):
    """
    Return list of files from given file_dir.
    E.g. file_ending="bed" to filter for .bed files.

    >>> test_dir = "test_data"
    >>> dir_get_files(test_dir, file_ending="bam")
    ['empty.bam', 'test_reads_chrM.bam']

    """

    from os import listdir
    from os.path import isfile, join
    dir_files = [f for f in listdir(file_dir) if isfile(join(file_dir, f))]
    if check:
        assert dir_files, "given directory \"%s\" contains no files" %(file_dir)
    # If filter for file ending true.
    if file_ending:
        new_files = []
        for df in dir_files:
            if re.search(".+\.%s" %(file_ending), df):
                new_files.append(df)
        if check:
            assert new_files, "no files left after filtering by file ending \"%s\"" %(file_ending)
        return sorted(new_files)
    else:
        return sorted(dir_files)


################################################################################

def shutil_copy_file(from_file, to_file):
    """
    Copy a file to another destination.
    This will overwrite to_file (!).

    """
    assert os.path.exists(from_file), "given file %s does not exist" %(from_file)
    from shutil import copyfile
    copyfile(from_file, to_file)


################################################################################

def get_filter_lists(list_f1_filter, list_f2_filter,
                     valid_filters_dic=False):
    """
    Check and get peakhood extract filter lists.

    """

    """
    Filter setting checks.

    valid_filters_dic:
        1 : transcript level filter
        2 : exon level filter

    """
    if not valid_filters_dic:
        valid_filters_dic = {"TSC": 1,
                             "EIR": 2,
                             "ISRN": 2,
                             "ISR": 1,
                             "ISRFC" : 1,
                             "SEO" : 1,
                             "FUCO": 1,
                             "TCOV": 1,
                             "TSL": 1
                             }

    f1_filters = []
    f2_filters = []

    if list_f1_filter:
        assert not list_found_duplicates(list_f1_filter), "--f1-filter list contains duplicates. Please provide each filter ID only once"
        for fid in list_f1_filter:
            assert fid in valid_filters_dic, "invalid --f1-filter ID given (%s)" %(fid)
            f1_filters.append(fid)
    else:
        f1_filters = ['TSC']

    if list_f2_filter:
        assert not list_found_duplicates(list_f2_filter), "--f2-filter list contains duplicates. Please provide each filter ID only once"
        for fid in list_f2_filter:
            assert fid in valid_filters_dic, "invalid --f2-filter ID given (%s)" %(fid)
            f2_filters.append(fid)
    else:
        f2_filters = ['EIR', 'ISRN', 'ISR', 'ISRFC', 'SEO', 'FUCO', 'TCOV']

    return f1_filters, f2_filters


################################################################################

def get_tsl_score(tsl, gc_basic, ccds):
    """
    Get score from TSL flag.

    Quality tags in (Ensembl) GTF:
    CCDS:
        Member of the consensus CDS gene set, confirming coding regions
        between ENSEMBL, UCSC, NCBI and HAVANA.
    basic:
        Identifies a subset of representative transcripts for each gene;
        prioritises full-length protein coding transcripts over partial
        or non-protein coding transcripts within the same gene, and
        intends to highlight those transcripts that will be useful
        to the majority of users.

    >>> tsl = "3 (assigned to previous version 5)"
    >>> get_tsl_score(tsl, False, False)
    14
    >>> tsl = "1"
    >>> get_tsl_score(tsl, True, True)
    32
    >>> tsl = "NA"
    >>> get_tsl_score(tsl, False, False)
    4

    """
    assert tsl, "given tsl empty"

    tag2sc_dic = {
        "1" : 24,
        "2" : 20,
        "3" : 16,
        "4" : 12,
        "5" : 8,
        "NA" : 4 }

    tsl_sc = 0
    if re.search('assigned', tsl):
        tsl_sc -= 2
        m = re.search('(.+) \(assigned', tsl)
        tsl_sc += tag2sc_dic[m.group(1)]
    else:
        tsl_sc += tag2sc_dic[tsl]
    if gc_basic:
        tsl_sc += 3
    if ccds:
        tsl_sc += 5
    return tsl_sc


################################################################################

def bed_extract_sequences_from_2bit(in_bed, out_fa, in_2bit,
                                    lc_repeats=False,
                                    convert_to_rna=False):
    """
    Extract sequences from genome (provide genome .2bit file).
    twoBitToFa executable needs to be in PATH. Store extracted
    sequences in out_fa.

    convert_to_rna:
        If true, read in extracted sequences and convert to RNA.
    lc_repeats:
        If True, do not convert repeat regions to uppercase and output.

    >>> in_bed = "test_data/test_seq_extr.sites.bed"
    >>> tmp_2bit_fa = "test_data/test_seq_extr.sites.2bit.tmp.fa"
    >>> tmp_seq_fa = "test_data/test_seq_extr.sites.seq.tmp.fa"
    >>> exp_fa = "test_data/test_seq_extr.sites.exp.fa"
    >>> in_fa = "test_data/test_seq_extr.sequences.fa"
    >>> in_2bit = "test_data/test_seq_extr.sequences.2bit"
    >>> id2row_dic = bed_read_rows_into_dic(in_bed)
    >>> seqs_dic = read_fasta_into_dic(in_fa, dna=True)
    >>> id2seq_dic = extract_transcript_sequences(id2row_dic, seqs_dic, revcom=True)
    >>> fasta_output_dic(id2seq_dic, tmp_seq_fa)
    >>> bed_extract_sequences_from_2bit(in_bed, tmp_2bit_fa, in_2bit)
    >>> diff_two_files_identical(tmp_seq_fa, exp_fa)
    True
    >>> diff_two_files_identical(tmp_2bit_fa, exp_fa)
    True

    """
    # Check for twoBitToFa.
    assert is_tool("twoBitToFa"), "twoBitToFa not in PATH"

    # Run twoBitToFa and check.
    check_cmd = "twoBitToFa"
    if not lc_repeats:
        check_cmd += " -noMask"
    check_cmd += " -bed=" + in_bed + " " + in_2bit + " " + out_fa
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "twoBitToFa is complaining:\n%s\n%s" %(check_cmd, output)
    if convert_to_rna:
        # Read in tmp_fa into dictionary (this also converts sequences to RNA).
        seqs_dic = read_fasta_into_dic(out_fa)
        # Output RNA sequences.
        fasta_output_dic(seqs_dic, out_fa,
                         split=True)


################################################################################

def get_chr_ids_from_bam_file(in_bam, k_top_reads=1000):
    """
    Return dictionary with chromosome IDs from BAM file in_bam.
    Look only at the # k_top_reads read entries for
    extracting chromosome IDs.

    >>> in_bam = "test_data/test_reads_chrM.bam"
    >>> get_chr_ids_from_bam_file(in_bam)
    {'chrM': 1}

    """

    assert os.path.exists(in_bam), "in_bam %s does not exist" %(in_bam)

    check_cmd = "samtools view " + in_bam + " | head -" + str(k_top_reads) + " | cut -f 3 | sort | uniq "
    output = subprocess.getoutput(check_cmd)

    chr_ids_dic = {}

    for line in output.split('\n'):
        cols = line.strip().split("\t")
        chr_id = cols[0]
        chr_ids_dic[chr_id] = 1

    return chr_ids_dic


################################################################################

def get_chromosome_lengths_from_2bit(in_2bit, out_lengths,
                                     std_chr_filter=False):
    """
    Get chromosome lengths from in_2bit .2bit file. Write lengths
    to out_lengths, with format:
    chr1	248956422
    chr10	133797422
    chr11	135086622
    ...
    Also return a dictionary with key=chr_id and value=chr_length.

    std_chr_filter:
        Filter / convert chromosome IDs with function check_convert_chr_id(),
        removing non-standard chromosomes, and convert IDs like 1,2,X,MT ..
        to chr1, chr2, chrX, chrM.

    """

    # Check for twoBitInfo.
    assert is_tool("twoBitInfo"), "twoBitInfo not in PATH"

    # Run twoBitInfo and check.
    check_cmd = "twoBitInfo " + in_2bit + " " + out_lengths
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "twoBitInfo is complaining:\n%s\n%s" %(check_cmd, output)

    # Read in lengths into dictionary.
    chr_len_dic = {}
    with open(out_lengths) as f:
        for line in f:
            row = line.strip()
            cols = line.strip().split("\t")
            chr_id = cols[0]
            chr_l = int(cols[1])
            # Check ID.
            if std_chr_filter:
                new_chr_id = check_convert_chr_id(chr_id)
                # If not standard chromosome ID or conversion failed, skip.
                if not new_chr_id:
                    continue
                else:
                    chr_id = new_chr_id
            assert chr_id not in chr_len_dic, "non-unique chromosome ID \"%s\" encountered in \"%s\"" %(chr_id, out_lengths)
            chr_len_dic[chr_id] = chr_l
    f.closed
    assert chr_len_dic, "chr_len_dic empty (\"%s\" empty? Chromosome IDs filter activated?)" %(out_lengths)

    return chr_len_dic


################################################################################

def gtf_get_transcript_lengths(in_gtf,
                               chr_ids_dic=None,
                               tr2exc_dic=None):
    """
    Get transcript lengths (= length of their exons, not unspliced length!)
    from GTF file.

    tr2exc_dic:
    Optionally provide a transcript ID to exon count dictionary for counting
    transcript exons.

    >>> in_gtf = "test_data/map_test_in.gtf"
    >>> gtf_get_transcript_lengths(in_gtf)
    {'ENST001': 2000, 'ENST002': 2000}

    """
    # Transcript ID to exonic length dictionary.
    tr2len_dic = {}
    # Open GTF either as .gz or as text file.
    if re.search(".+\.gz$", in_gtf):
        f = gzip.open(in_gtf, 'rt')
    else:
        f = open(in_gtf, "r")
    for line in f:
        # Skip header.
        if re.search("^#", line):
            continue
        cols = line.strip().split("\t")
        chr_id = cols[0]
        feature = cols[2]
        feat_s = int(cols[3])
        feat_e = int(cols[4])
        infos = cols[8]
        if not feature == "exon":
            continue
        # Extract transcript ID.
        m = re.search('transcript_id "(.+?)"', infos)
        assert m, "transcript_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        tr_id = m.group(1)
        # Sum up length.
        ex_len = feat_e - feat_s + 1
        if not tr_id in tr2len_dic:
            tr2len_dic[tr_id] = ex_len
        else:
            tr2len_dic[tr_id] += ex_len
        if tr2exc_dic is not None:
            if not tr_id in tr2exc_dic:
                tr2exc_dic[tr_id] = 1
            else:
                tr2exc_dic[tr_id] += 1
        if chr_ids_dic is not None:
            chr_ids_dic[chr_id] = 1

    f.close()
    assert tr2len_dic, "No IDs read into dictionary (input file \"%s\" empty or malformatted?)" % (in_gtf)
    return tr2len_dic


################################################################################

def rem_exb_pairs_exb_dist(id2ids_dic, id2exids_dic,
                           id2gen_se_dic, id2gen_cp_dic, exid2gen_se_dic,
                           max_exb_dist=10,
                           exid2exnr_dic=False,
                           exid2trid_dic=False):
    """
    Remove exon border pairs from id2ids_dic if the sites are too far
    away from the matching exon borders (supporting the pair).

    >>> id2ids_dic = {'id1': ['id2'], 'id2': ['id1']}
    >>> id2exids_dic = {'id1': ['t1_e1'], 'id2': ['t1_e2']}
    >>> id2gen_se_dic = {'id1': [1980, 1990], 'id2': [3005, 3025]}
    >>> id2gen_cp_dic = {'id1': 1985, 'id2': 3015}
    >>> exid2gen_se_dic = {'t1_e1': [1000, 2000], 't1_e2': [3000, 4000], 't1_e3': [5000, 6000]}
    >>> exid2exnr_dic = {'t1_e1': 1, 't1_e2': 2, 't1_e3' : 3}
    >>> exid2trid_dic = {'t1_e1': 't1', 't1_e2': 't1', 't1_e3': 't1'}
    >>> rem_exb_pairs_exb_dist(id2ids_dic, id2exids_dic, id2gen_se_dic, id2gen_cp_dic, exid2gen_se_dic, max_exb_dist=10, exid2exnr_dic=exid2exnr_dic, exid2trid_dic=exid2trid_dic)
    {'id1': ['id2'], 'id2': ['id1']}
    >>> rem_exb_pairs_exb_dist(id2ids_dic, id2exids_dic, id2gen_se_dic, id2gen_cp_dic, exid2gen_se_dic, max_exb_dist=9, exid2exnr_dic=exid2exnr_dic, exid2trid_dic=exid2trid_dic)
    {}
    >>> id2exids_dic = {'id1': ['t1_e1'], 'id2': ['t1_e3']}
    >>> rem_exb_pairs_exb_dist(id2ids_dic, id2exids_dic, id2gen_se_dic, id2gen_cp_dic, exid2gen_se_dic, max_exb_dist=10, exid2exnr_dic=exid2exnr_dic, exid2trid_dic=exid2trid_dic)
    {}
    >>> id2ids_dic = {'id1': ['id2', 'id3'], 'id2': ['id1'], 'id3': ['id1']}
    >>> id2exids_dic = {'id1': ['t1_e1'], 'id2': ['t1_e2'], 'id3': ['t1_e1']}
    >>> id2gen_se_dic = {'id1': [1980, 1990], 'id2': [3005, 3025], 'id2': [1970, 1980]}
    >>> id2gen_cp_dic = {'id1': 1985, 'id2': 3015, 'id1': 1975}
    >>> rem_exb_pairs_exb_dist(id2ids_dic, id2exids_dic, id2gen_se_dic, id2gen_cp_dic, exid2gen_se_dic, max_exb_dist=10, exid2exnr_dic=exid2exnr_dic, exid2trid_dic=exid2trid_dic)
    {'id1': ['id2'], 'id2': ['id1']}

    """

    assert exid2exnr_dic, "exid2exnr_dic empty"
    assert exid2trid_dic, "exid2trid_dic empty"

    rem_sids_list = []

    for sid1 in id2ids_dic:
        new_con_list = []
        exl1 = id2exids_dic[sid1]
        for sid2 in id2ids_dic[sid1]:
            exl2 = id2exids_dic[sid2]
            # Compare the two exon ID lists.
            for exid1 in exl1:
                exnr1 = exid2exnr_dic[exid1]
                trid1 = exid2trid_dic[exid1]
                for exid2 in exl2:
                    exnr2 = exid2exnr_dic[exid2]
                    trid2 = exid2trid_dic[exid2]
                    # Check distance of site to exon borders.
                    if trid1 == trid2 and abs(exnr1-exnr2) == 1:
                        # Orientation of the two sites.
                        cp_diff = id2gen_cp_dic[sid1] - id2gen_cp_dic[sid2]
                        sid1_us = False  # sid1 upstream of sid2 ?
                        if cp_diff == 0:
                            assert False, "two exon border pair site IDs %s,%s with same genomic center positions encountered (%i)" %(sid1, sid2, id2gen_cp_dic[sid1])
                        elif cp_diff < 1:
                            sid1_us = True

                        # SID1 distance to exon okay?
                        sid1_check = False
                        if sid1_us:
                            # Distance of site and exon end.
                            end_dist = exid2gen_se_dic[exid1][1] - id2gen_se_dic[sid1][1]
                            if end_dist <= max_exb_dist:
                                sid1_check = True
                        else:
                            # Distance of site and exon start.
                            end_dist = id2gen_se_dic[sid1][0] - exid2gen_se_dic[exid1][0]
                            if end_dist <= max_exb_dist:
                                sid1_check = True
                        if not sid1_check:
                            continue

                        # SID2 distance to exon okay?
                        sid2_check = False
                        if sid1_us:
                            # Distance of site and exon start.
                            end_dist = id2gen_se_dic[sid2][0] - exid2gen_se_dic[exid2][0]
                            if end_dist <= max_exb_dist:
                                sid2_check = True
                        else:
                            # Distance of site and exon end.
                            end_dist = exid2gen_se_dic[exid2][1] - id2gen_se_dic[sid2][1]
                            if end_dist <= max_exb_dist:
                                sid2_check = True
                        if sid1_check and sid2_check:
                            new_con_list.append(sid2)

        if new_con_list:
            id2ids_dic[sid1] = new_con_list
        else:
            rem_sids_list.append(sid1)

    for rem_sid in rem_sids_list:
        del id2ids_dic[rem_sid]

    return id2ids_dic


################################################################################

def check_neighbor_ex_lists(exl1, exl2,
                            exid2exnr_dic=False,
                            exid2trid_dic=False):
    """
    Check if two exon ID lists contain neighboring exons.

    >>> exl1 = ["t1_e5", "t2_e4", "t3_e5"]
    >>> exl2 = ["t1_e5", "t2_e5"]
    >>> check_neighbor_ex_lists(exl1, exl2)
    True
    >>> exl1 = ["t1_e5", "t2_e4", "t4_e5"]
    >>> exl2 = ["t1_e5", "t2_e2"]
    >>> check_neighbor_ex_lists(exl1, exl2)
    False

    """

    if exid2exnr_dic and exid2trid_dic:
        for exid1 in exl1:
            exnr1 = exid2exnr_dic[exid1]
            trid1 = exid2trid_dic[exid1]
            for exid2 in exl2:
                exnr2 = exid2exnr_dic[exid2]
                trid2 = exid2trid_dic[exid2]
                if trid1 == trid2:
                    if abs(exnr1-exnr2) == 1:
                        return True
    else:
        for exid1 in exl1:
            m = re.search("(.+)_e(\d+)$", exid1)
            trid1 = m.group(1)
            exnr1 = int(m.group(2))
            for exid2 in exl2:
                m = re.search("(.+)_e(\d+)$", exid2)
                trid2 = m.group(1)
                exnr2 = int(m.group(2))
                if trid1 == trid2:
                    if abs(exnr1-exnr2) == 1:
                        return True
    return False


################################################################################

def rem_exb_pairs_no_neighbor_ex(id2ids_dic, id2exids_dic,
                                 exid2exnr_dic=False,
                                 exid2trid_dic=False):
    """
    Remove exon border pairs from id2ids_dic which are not
    on neighboring exons.

    >>> id2ids_dic = {'id1': ['id2'], 'id2': ['id1', 'id3'], 'id3': ['id2']}
    >>> id2exids_dic = {'id1' : ['t1_e5', 't2_e4'], 'id2': ['t1_e5'], 'id3': ['t1_e6']}
    >>> rem_exb_pairs_no_neighbor_ex(id2ids_dic, id2exids_dic)
    {'id2': ['id3'], 'id3': ['id2']}

    """

    rem_sids_list = []

    for sid1 in id2ids_dic:
        new_con_list = []
        exl1 = id2exids_dic[sid1]
        for sid2 in id2ids_dic[sid1]:
            exl2 = id2exids_dic[sid2]
            check = check_neighbor_ex_lists(exl1, exl2,
                                    exid2exnr_dic=exid2exnr_dic,
                                    exid2trid_dic=exid2trid_dic)
            if check:
                new_con_list.append(sid2)
        if new_con_list:
            id2ids_dic[sid1] = new_con_list
        else:
            rem_sids_list.append(sid1)

    for rem_sid in rem_sids_list:
        del id2ids_dic[rem_sid]

    return id2ids_dic


################################################################################

def extract_transcript_sequences(bed_dic, seq_dic,
                                 ext_mode=1,
                                 ext_lr=False,
                                 revcom=False,
                                 ids_dic=False,
                                 out_bed=False,
                                 out_bed_add_sc_dic=False,
                                 full_hits_only=False):
    """
    Given a dictionary with bed regions (region ID -> BED row) and a
    sequence dictionary (Sequence ID -> sequence), extract the BED region
    sequences and return in new dictionary (region ID -> region sequence).

    ext_mode:
        1: whole site
        2: center position site
        3: upstream end position site
    ext_lr:
        Optionally, extend regions by ext_lr nt (up- and downstream).
        In case full extension is not possible, use maximum extension possible.
    revcom:
        if revcom=True and strand of bed_dic region is "-", return the reverse
        complement of the region sequence.
    ids_dic:
        IDs to extract sequences for from bed_dic.
    full_hits_only:
        Set full_hits_only=True to only recover full hits.
    out_bed:
        Output BED file path to store regions for which sequences were
        extracted in.
    out_bed_add_sc_dic:
        Region ID to score mapping, with score to be added to out_bed ID column
        4, so column 4 ID changes from "id1" to "id1,5"

    >>> seq_dic = {"T1" : "AAAACCCCGGGGTTTT", "T2" : "ATATACACAGAGCGCGCTCTGTGT"}
    >>> bed_dic = {"S1" : "T1\\t4\\t8\\tS1\\t0\\t+", "S2" : "T2\\t6\\t8\\tS2\\t0\\t+"}
    >>> extract_transcript_sequences(bed_dic, seq_dic, ext_lr=2)
    {'S1': 'AACCCCGG', 'S2': 'ACACAG'}
    >>> extract_transcript_sequences(bed_dic, seq_dic, ext_lr=5, full_hits_only=True)
    {'S2': 'TATACACAGAGC'}
    >>> bed_dic = {"S1" : "T1\\t4\\t8\\tS1\\t0\\t+", "S2" : "T2\\t6\\t8\\tS2\\t0\\t+"}
    >>> extract_transcript_sequences(bed_dic, seq_dic, ext_lr=2, ext_mode=2)
    {'S1': 'CCCCG', 'S2': 'CACAG'}

    """
    id2seq_dic = {}

    if out_bed:
        OUT2BED = open(out_bed,"w")

    # Process .bed regions.
    for reg_id in bed_dic:
        cols = bed_dic[reg_id].split("\t")
        seq_id = cols[0]
        reg_s = int(cols[1])
        reg_e = int(cols[2])
        reg_sc = cols[4]
        reg_pol = cols[5]
        if ids_dic:
            if reg_id not in ids_dic:
                continue
        assert seq_id in seq_dic, "sequence ID \"%s\" not found in given sequence dictionary" %(seq_id)
        seq = seq_dic[seq_id]

        # Update region lengths.
        if ext_mode == 1:
            new_s = reg_s
            new_e = reg_e
        elif ext_mode == 2:
            new_e = get_center_position(reg_s, reg_e)
            new_s = new_e - 1
        elif ext_mode == 3:
            new_s = reg_s
            new_e = reg_s + 1
            if reg_pol == "-":
                new_s = reg_e - 1
                new_e = reg_e
        else:
            assert False, "invalid ext_mode set"

        # Expected length.
        exp_l = new_e - new_s

        # Adjust if given start or end is out of bounds.
        if new_s < 0:
            new_s = 0
        if new_e > len(seq):
            new_e = len(seq)

        # If region should be extended up- and downstream by ext_lr.
        if ext_lr:
            new_s = new_s - ext_lr
            new_e = new_e + ext_lr
            exp_l = new_e - new_s
            # If start or end is out of bounds after extension.
            if new_s < 0:
                new_s = 0
            if new_e > len(seq):
                new_e = len(seq)
        reg_seq = seq[new_s:new_e]
        reg_l = len(reg_seq)
        if full_hits_only:
            if not reg_l == exp_l:
                continue
        if revcom:
            if reg_pol == "-":
                id2seq_dic[reg_id] = revcom_seq(reg_seq)
            else:
                id2seq_dic[reg_id] = reg_seq
        else:
            id2seq_dic[reg_id] = reg_seq
        if out_bed:
            if out_bed_add_sc_dic:
                reg_id = reg_id + "," + str(out_bed_add_sc_dic[reg_id])
            OUT2BED.write("%s\t%i\t%i\t%s\t%s\t%s\n" %(seq_id, new_s, new_e, reg_id, reg_sc, reg_pol))
    if out_bed:
        OUT2BED.close()
    assert id2seq_dic, "no sequences extracted"
    return id2seq_dic


################################################################################

def check_string_in_file(in_file, string):
    """
    Check whether a string is found inside a text file.
    Return True if found, else False.

    >>> in_file = "test_data/test_run.log"
    >>> check_string_in_file(in_file, "AssertionError")
    True
    >>> in_file = "test_data/empty_file"
    >>> check_string_in_file(in_file, "AssertionError")
    False

    """

    found = False

    f = open(in_file, "r")
    for line in f:
        if string in line:
            found = True
            break
    f.close()

    return found


################################################################################

def revcom_seq(seq,
               upper=False,
               convert_to_rna=False):
    """
    Return reverse complement to seq. By default, convert seq to uppercase
    and translate to DNA.

    # Convert to RNA.
    if convert_to_rna:
        new_seq_rna = new_seq.replace("T","U").replace("t","u")
        new_seq = new_seq_rna

    >>> seq = "AAACAGatt"
    >>> revcom_seq(seq)
    'aatCTGTTT'
    >>> revcom_seq(seq, upper=True)
    'AATCTGTTT'
    >>> revcom_seq(seq, convert_to_rna=True)
    'aauCUGUUU'

    """
    assert seq, "given sequence empty"
    # Make uppercase and convert to DNA.
    if upper:
        seq = seq[::-1].upper().replace("U","T")
    else:
        seq = seq[::-1].replace("U","T").replace("u","t")
    intab = "ACGTacgt"
    outtab = "TGCAtgca"
    # If RNA revcom should be output.
    if convert_to_rna:
        seq = seq.replace("T","U").replace("t","u")
        intab = "ACGUacgu"
        outtab = "UGCAugca"
    # Make revcom.
    transtab = str.maketrans(intab, outtab)
    rc_seq = seq.translate(transtab)
    return rc_seq


################################################################################

def fasta_output_dic(fasta_dic, fasta_out,
                     split=False,
                     out_ids_dic=False,
                     header_add_sc_dic=False,
                     to_upper=False,
                     split_size=60):
    """
    Output FASTA sequences dictionary (sequence_id -> sequence) to fasta_out.

    split:
        Split FASTA sequence for output to file
    split_size:
        Split size (row width)
    to_upper:
        Convert sequences to uppercase.
    out_ids_dic:
        IDs to output dictionary.
    header_add_sc_dic:
        ID to scoring mapping.
        Add a score to the header, so header format changes from "id1"
        to "id1,10"

    >>> fasta_dic = {'seq1': 'ACGTACGTACGTAC', 'seq2': 'ACGT'}
    >>> split_size = 4
    >>> fasta_exp = "test_data/test5.exp.fa"
    >>> fasta_out = "test_data/test5.tmp.fa"
    >>> fasta_output_dic(fasta_dic, fasta_out, split=True, split_size=split_size)
    >>> diff_two_files_identical(fasta_exp, fasta_out)
    True

    """
    # Check.
    assert fasta_dic, "given dictionary fasta_dic empty"
    # Write sequences to FASTA file.
    OUTFA = open(fasta_out,"w")
    for seq_id in fasta_dic:
        seq = fasta_dic[seq_id]
        if out_ids_dic:
            if seq_id not in out_ids_dic:
                continue
        if to_upper:
            seq = seq.upper()
        out_id = seq_id
        if header_add_sc_dic:
            out_id = out_id + "," + str(header_add_sc_dic[seq_id])
        if split:
            OUTFA.write(">%s\n" %(out_id))
            for i in range(0, len(seq), split_size):
                OUTFA.write("%s\n" %((seq[i:i+split_size])))
        else:
            OUTFA.write(">%s\n%s\n" %(out_id, seq))
    OUTFA.close()


################################################################################

def read_fasta_into_dic(fasta_file,
                        seqs_dic=False,
                        ids_dic=False,
                        dna=False,
                        report=1,
                        all_uc=False,
                        empty_check=True,
                        id_check=True,
                        skip_data_id="set",
                        skip_n_seqs=True):
    """
    Read in FASTA sequences, store in dictionary and return dictionary.
    FASTA file can be plain text or gzipped (watch out for .gz ending).

    >>> test_fasta = "test_data/test.fa"
    >>> read_fasta_into_dic(test_fasta)
    {'seq1': 'acguACGUacgu', 'seq2': 'ugcaUGCAugcaACGUacgu'}

    """
    if not seqs_dic:
        seqs_dic = {}
    seq_id = ""

    # Open FASTA either as .gz or as text file.
    if re.search(".+\.gz$", fasta_file):
        f = gzip.open(fasta_file, 'rt')
    else:
        f = open(fasta_file, "r")
    for line in f:
        if re.search(">.+", line):
            m = re.search(">(.+)", line)
            seq_id = m.group(1)
            if id_check:
                assert seq_id not in seqs_dic, "non-unique FASTA header \"%s\" in \"%s\"" % (seq_id, fasta_file)
            if ids_dic:
                if seq_id in ids_dic:
                    seqs_dic[seq_id] = ""
            else:
                seqs_dic[seq_id] = ""
        elif re.search("[ACGTUN]+", line, re.I):
            m = re.search("([ACGTUN]+)", line, re.I)
            seq = m.group(1)
            if seq_id in seqs_dic:
                if dna:
                    # Convert to DNA, concatenate sequence.
                    seq = seq.replace("U","T").replace("u","t")
                else:
                    # Convert to RNA, concatenate sequence.
                    seq = seq.replace("T","U").replace("t","u")
                if all_uc:
                    seq = seq.upper()
                seqs_dic[seq_id] += seq
    f.close()

    # Check if sequences read in.
    if empty_check:
        assert seqs_dic, "no sequences read in (input FASTA file \"%s\" empty or mal-formatted?)" %(fasta_file)
    # If sequences with N nucleotides should be skipped.
    c_skipped_n_ids = 0
    if skip_n_seqs:
        del_ids = []
        for seq_id in seqs_dic:
            seq = seqs_dic[seq_id]
            if re.search("N", seq, re.I):
                if report == 1:
                    print ("WARNING: sequence with seq_id \"%s\" in file \"%s\" contains N nucleotides. Discarding sequence ... " % (seq_id, fasta_file))
                c_skipped_n_ids += 1
                del_ids.append(seq_id)
        for seq_id in del_ids:
            del seqs_dic[seq_id]
        assert seqs_dic, "no sequences remaining after deleting N containing sequences (input FASTA file \"%s\")" %(fasta_file)
        if c_skipped_n_ids:
            if report == 2:
                print("# of N-containing %s regions discarded:  %i" %(skip_data_id, c_skipped_n_ids))
    return seqs_dic


################################################################################

def get_seqs_dic_repeat_region_ratios(seqs_dic):
    """
    Given a dictionary of sequences, calculate the repeat region content /
    ratio for each site. Return dictionary of repeat region ratios.

    >>> seqs_dic = {'s1' : "ACGUacgu", 's2' : "acnacnACNacn", 's3' : "A", 's4' : "u"}
    >>> get_seqs_dic_repeat_region_ratios(seqs_dic)
    {'s1': 0.5, 's2': 0.75, 's3': 0.0, 's4': 1.0}

    """
    assert seqs_dic, "seqs_dic empty"

    ratios_dic = {}
    for seq_id in seqs_dic:
        seq = seqs_dic[seq_id]
        seq_l = len(seq)
        c_r = 0
        for nt in seq:
            if nt.islower():
                c_r += 1
        r_ratio = c_r / seq_l
        ratios_dic[seq_id] = r_ratio
    assert ratios_dic, "ratios_dic empty"
    return ratios_dic


################################################################################

def bed_get_chromosome_ids(bed_file):
    """
    Read in .bed file, return chromosome IDs (column 1 IDs).
    Return dic with chromosome ID -> count mapping.

    >>> test_file = "test_data/test6.bed"
    >>> bed_get_chromosome_ids(test_file)
    {'chr1': 2, 'chr2': 2, 'chr3': 1}

    """
    ids_dic = {}
    with open(bed_file) as f:
        for line in f:
            row = line.strip()
            cols = line.strip().split("\t")
            chr_id = cols[0]
            if chr_id in ids_dic:
                ids_dic[chr_id] += 1
            else:
                ids_dic[chr_id] = 1
    f.closed
    assert ids_dic, "No chromosome IDs read into dictionary (input file \"%s\" empty or malformatted?)" % (bed_file)
    return ids_dic


################################################################################

def bed_get_region_lengths(bed_file):
    """
    Read in .bed file, store and return region lengths in dictionary.
    key   :  region ID (.bed col4)
    value :  region length (.bed col3-col2)

    >>> test_file = "test_data/test4.bed"
    >>> bed_get_region_lengths(test_file)
    {'tr1_e1': 14, 'tr2_e1': 30}

    """
    id2len_dic = {}
    with open(bed_file) as f:
        for line in f:
            row = line.strip()
            cols = line.strip().split("\t")
            site_s = int(cols[1])
            site_e = int(cols[2])
            site_id = cols[3]
            site_l = site_e - site_s
            assert site_id not in id2len_dic, "column 4 IDs not unique in given .bed file \"%s\"" %(bed_file)
            id2len_dic[site_id] = site_l
    f.closed
    assert id2len_dic, "No IDs read into dictionary (input file \"%s\" empty or malformatted?)" % (in_bed)
    return id2len_dic


################################################################################

def bed_convert_transcript_to_genomic_sites(in_bed, in_gtf, out_bed,
                                            site2hitc_dic=None,
                                            chr_id_style=0,
                                            out_folder=False):
    """
    Dependencies:
    bedtools (tested with 2.29.0)
    gzip

    Convert in_bed .bed file with transcript sites into genomic coordinates
    sites file. in_bed column 1 transcript IDs have to be present in
    in_gtf GTF file, from which genomic exon coordinates of the transcript
    will get extracted.

    site2hitc_dic:
    A site2hitc_dic can be given, where site ID to hit count will be stored
    for usage outside the function.

    Output:
    By default output to out_bed file, using id_p1, id_p2 IDs.
    If out_folder=True, use out_bed name as folder name.
    In this case, output these files to folder:
    exon_regions_genome.bed
    exon_regions_transcript.bed
    complete_hits.bed
    split_hits.bed
    all_hits.bed

    >>> test_gtf = "test_data/test_tr2gen.gtf"
    >>> test_in_bed = "test_data/test_tr2gen.bed"
    >>> test_out_exp_bed = "test_data/test_tr2gen.exp.bed"
    >>> test_out_tmp_bed = "test_data/test_tr2gen.tmp.bed"
    >>> bed_convert_transcript_to_genomic_sites(test_in_bed, test_gtf, test_out_tmp_bed, chr_id_style=1)
    >>> diff_two_files_identical(test_out_exp_bed, test_out_tmp_bed)
    True
    >>> test_out = "test_data/tr2gen_tmp_out"
    >>> test_out_tmp_bed = "test_data/tr2gen_tmp_out/all_hits.bed"
    >>> bed_convert_transcript_to_genomic_sites(test_in_bed, test_gtf, test_out, out_folder=True, chr_id_style=1)
    >>> diff_two_files_identical(test_out_exp_bed, test_out_tmp_bed)
    True

    """

    # Generate .tmp files.
    random_id = uuid.uuid1()
    tmp_bed = str(random_id) + ".tmp.bed"
    random_id = uuid.uuid1()
    tmp_out = str(random_id) + ".tmp.out"

    # Output files if output_folder=True.
    if out_folder:
        if not os.path.exists(out_bed):
            os.makedirs(out_bed)
    out_exon_regions_genome_bed = out_bed + "/" + "exon_regions_genome.bed"
    out_exon_regions_transcript_bed = out_bed + "/" + "exon_regions_transcript.bed"
    out_unique_hits_bed = out_bed + "/" + "unique_hits.bed"
    out_split_hits_bed = out_bed + "/" + "split_hits.bed"
    out_all_hits_bed = out_bed + "/" + "all_hits.bed"

    # Transcript IDs dic.
    tr_ids_dic = bed_get_chromosome_ids(in_bed)

    # Extract transcript exon regions from GTF and store as BED.
    gtf_extract_exon_bed(in_gtf, tmp_bed, tr_ids_dic=tr_ids_dic,
                         chr_id_style=chr_id_style)
    if out_folder:
        make_file_copy(tmp_bed, out_exon_regions_transcript_bed)

    # Get exon region lengths.
    exid2len_dic = bed_get_region_lengths(tmp_bed)

    # Get exon numbers for each transcript.
    tr_exc_dic = bed_get_transcript_exon_numbers(tmp_bed)

    # Read in exon region stats.
    id2chr_dic = {}
    id2s_dic = {}
    id2e_dic = {}
    id2pol_dic = {}
    exid2trid_dic = {}
    with open(tmp_bed) as f:
        for line in f:
            row = line.strip()
            cols = line.strip().split("\t")
            chr_id = cols[0]
            site_s = int(cols[1])
            site_e = int(cols[2])
            site_id = cols[3]
            site_pol = cols[5]
            id2chr_dic[site_id] = chr_id
            id2s_dic[site_id] = site_s
            id2e_dic[site_id] = site_e
            id2pol_dic[site_id] = site_pol
            if re.search(".+_e\d", site_id):
                m = re.search("(.+)_e\d", site_id)
                tr_id = m.group(1)
                exid2trid_dic[site_id] = tr_id
            else:
                assert False, "site ID \"%s\" missing added _e exon number" %(site_id)
    f.close()

    # Output exon regions with transcript coordinates.
    OUTBED = open(tmp_bed, "w")
    for tr_id in tr_exc_dic:
        ex_c = tr_exc_dic[tr_id]
        new_s = 0
        for i in range(ex_c):
            i += 1
            ex_id = tr_id + "_e" + str(i)
            gen_s = id2s_dic[ex_id]
            gen_e = id2e_dic[ex_id]
            ex_len = gen_e - gen_s
            tr_s = new_s
            tr_e = new_s + ex_len
            OUTBED.write("%s\t%i\t%i\t%s\t0\t+\n" % (tr_id,tr_s,tr_e,ex_id))
            new_s = tr_e
    OUTBED.close()

    if out_folder:
        make_file_copy(tmp_bed, out_exon_regions_genome_bed)

    # Overlap in_bed with tmp_bed.
    params = "-wb"
    intersect_bed_files(in_bed, tmp_bed, params, tmp_out,
                        sorted_out=True)

    # Read in transcript site overlaps with transcript exon regions.
    site2c_dic = {}
    # Dictionaries for later outputting unique + split hits separately.
    siteid2pol_dic = {}
    siteid2sc_dic = {}
    partid2chrse_dic = {}
    with open(tmp_out) as f:
        for line in f:
            row = line.strip()
            cols = line.strip().split("\t")
            tr_id = cols[0]
            part_s = int(cols[1])
            part_e = int(cols[2])
            site_id = cols[3]
            site_sc = cols[4]
            ex_s = int(cols[7])
            ex_e = int(cols[8])
            ex_id = cols[9]
            ex_pol = id2pol_dic[ex_id]
            siteid2pol_dic[site_id] = ex_pol
            siteid2sc_dic[site_id] = site_sc
            if site_id in site2c_dic:
                site2c_dic[site_id] += 1
            else:
                site2c_dic[site_id] = 1
            # Hit part number.
            hit_c = site2c_dic[site_id]
            # Calculate genomic hit coordinates.
            # Plus strand case.
            gen_s = id2s_dic[ex_id] + part_s - ex_s
            gen_e = id2s_dic[ex_id] + part_e - ex_s
            # Minus strand case.
            if ex_pol == "-":
                gen_s = id2e_dic[ex_id] - part_e + ex_s
                gen_e = id2e_dic[ex_id] - part_s + ex_s
            # part ID.
            part_id = site_id + "_p" + str(hit_c)
            # Store chrse for each part ID.
            chrse = "%s\t%i\t%i" %(id2chr_dic[ex_id],gen_s,gen_e)
            partid2chrse_dic[part_id] = "%s\t%i\t%i" %(id2chr_dic[ex_id],gen_s,gen_e)

    # Produce seperate output files for unique + split hits.
    all_hits_bed = out_bed
    if out_folder:
        all_hits_bed = out_all_hits_bed
    ALLBED = open(all_hits_bed, "w")
    if out_folder:
        UNIBED = open(out_unique_hits_bed, "w")
        SPLBED = open(out_split_hits_bed, "w")
        for site_id in site2c_dic:
            hit_c = site2c_dic[site_id]
            if site2hitc_dic is not None:
                site2hitc_dic[site_id] = hit_c
            site_pol = siteid2pol_dic[site_id]
            site_sc = siteid2sc_dic[site_id]
            # For unique hit use site ID, for split hits use part IDs.
            if hit_c == 1:
                # Unique hits.
                part_id = site_id + "_p1"
                UNIBED.write("%s\t%s\t%s\t%s\n" %(partid2chrse_dic[part_id],site_id,site_sc,site_pol))
            else:
                # Split hits.
                for i in range(hit_c):
                    i += 1
                    part_id = site_id + "_p" + str(i)
                    SPLBED.write("%s\t%s\t%s\t%s\n" %(partid2chrse_dic[part_id],part_id,site_sc,site_pol))
    # Output all hits.
    for site_id in site2c_dic:
        hit_c = site2c_dic[site_id]
        if site2hitc_dic is not None:
            site2hitc_dic[site_id] = hit_c
        site_pol = siteid2pol_dic[site_id]
        site_sc = siteid2sc_dic[site_id]
        # For unique hit use site ID, for split hits use part IDs.
        if hit_c == 1:
            # Unique hits.
            part_id = site_id + "_p1"
            ALLBED.write("%s\t%s\t%s\t%s\n" %(partid2chrse_dic[part_id],site_id,site_sc,site_pol))
        else:
            # Split hits.
            for i in range(hit_c):
                i += 1
                part_id = site_id + "_p" + str(i)
                ALLBED.write("%s\t%s\t%s\t%s\n" %(partid2chrse_dic[part_id],part_id,site_sc,site_pol))

    # Delete tmp files.
    if os.path.exists(tmp_bed):
        os.remove(tmp_bed)
    if os.path.exists(tmp_out):
        os.remove(tmp_out)


################################################################################

def move_rename_file(in_file, out_file):
    """
    Move / rename in_file to out_file.

    """
    check_cmd = "mv " + in_file + " " + out_file
    assert in_file != out_file, "mv does not like to mv file into same file (%s)" %(check_cmd)
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "mv did not like your input (in_file: %s, out_file: %s):\n%s" %(in_file, out_file, output)


################################################################################

def touch_file(in_file):
    """
    Create an empty file.

    """
    assert not os.path.exists(in_file), "file %s to create already exists" %(in_file)
    check_cmd = "touch " + file
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "touch did not like your input (in_file: %s):\n%s" %(in_file, output)


################################################################################

def create_empty_file(in_file, check=True):
    """
    Create an empty file.

    check:
        If False do not complain if file exists, but overwrite the file.

    """
    if check:
        assert not os.path.exists(in_file), "file %s to create already exists" %(in_file)
    EMPTYOUT = open(in_file, "w")
    EMPTYOUT.close()
    if check:
        assert os.path.exists(in_file), "created file %s not found" %(in_file)


################################################################################

def read_settings_into_dic(settings_file,
                           check=True,
                           val_col2=False):
    """
    Read settings file content into dictionary.
    Each row expected to have following format:
    setting_id<tab>setting_value
    Skip rows with > 2 entries.
    Dictionary format: str(col1) -> str(col2)

    >>> test_in = "test_data/test_settings.out"
    >>> read_settings_into_dic(test_in)
    {'peyote': '20.5', 'china_white': '43.1', 'bolivian_marching_powder': '1000.0'}

    """
    assert os.path.isfile(settings_file), "file %s does not exist" %(settings_file)
    set_dic = {}
    with open(settings_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            settings_id = cols[0]
            if val_col2:
                settings_val = cols[2]
            else:
                settings_val = cols[1]
            if settings_id not in set_dic:
                set_dic[settings_id] = settings_val
            else:
                assert False, "settings ID %s appears > 1 in given settings file" %(settings_id)
    f.closed
    if check:
        assert set_dic, "set_dic empty (nothing read in?)"
    return set_dic


################################################################################

def get_transcript_sequences_from_gtf(in_gtf, in_2bit,
                                      lc_repeats=False,
                                      correct_min_ex_order=False,
                                      tr2exc_dic=False,
                                      tr_ids_dic=False,
                                      chr_id_style=1,
                                      tmp_out_folder=False):
    """
    Get spliced transcript sequences based on in_gtf annotations. For
    transcripts with > 1 exon, concatenate the exon sequences to build
    the transcript sequence. If one exon is missing / not extracted or
    if extracted lengths don't fit, the transcript sequence will be
    skipped / not output.
    Return dictionary with transcript_id -> sequence mapping.

    correct_min_ex_order:
        Set True if minus strand exon order should be corrected.
    tr2exc_dic:
        Transcript ID to exon count mapping.
    tr_ids_dic:
        Defines transcript IDs for which sequence should be extracted.

    """
    # Generate .tmp files.
    random_id = uuid.uuid1()
    tmp_bed = str(random_id) + ".tmp.bed"
    random_id = uuid.uuid1()
    tmp_fa = str(random_id) + ".tmp.fa"
    if tmp_out_folder:
        tmp_bed = tmp_out_folder + "/" + tmp_bed
        tmp_fa = tmp_out_folder + "/" + tmp_fa

    # Transcript sequences dic.
    tr_seqs_dic = {}

    # Extract transcript exon regions from GTF and store as BED.
    gtf_extract_exon_bed(in_gtf, tmp_bed,
                         chr_id_style=chr_id_style,
                         correct_min_ex_order=correct_min_ex_order,
                         tr2exc_dic=tr2exc_dic,
                         tr_ids_dic=tr_ids_dic)

    # Extract exon region sequences from .2bit.
    bed_extract_sequences_from_2bit(tmp_bed, tmp_fa, in_2bit,
                                    lc_repeats=lc_repeats)

    # Get transcript lengths from tmp_bed for comparison.
    tr_len_dic = bed_get_transcript_lengths_from_exon_regions(tmp_bed)
    # Get exon numbers for each transcript.
    tr_exc_dic = bed_get_transcript_exon_numbers(tmp_bed)

    # Read in sequences.
    exon_seqs_dic = read_fasta_into_dic(tmp_fa,
                                        skip_n_seqs=False)

    # Concatenate exon region sequences.
    for tr_id in tr_exc_dic:
        ex_c = tr_exc_dic[tr_id]
        for i in range(ex_c):
            i += 1
            ex_id = tr_id + "_e" + str(i)
            if ex_id in exon_seqs_dic:
                ex_seq = exon_seqs_dic[ex_id]
                if tr_id not in tr_seqs_dic:
                    tr_seqs_dic[tr_id] = ex_seq
                else:
                    tr_seqs_dic[tr_id] += ex_seq
            else:
                print("WARNING: no sequence extracted for exon ID \"%s\". Skipping \"%s\" .. " %(ex_id, tr_id))
                if tr_id in tr_seqs_dic:
                    del tr_seqs_dic[tr_id]
                break
    # Checks.
    assert tr_seqs_dic, "tr_seqs_dic empty (no FASTA sequences extracted?)"
    for tr_id in tr_seqs_dic:
        tr_len = len(tr_seqs_dic[tr_id])
        exp_len = tr_len_dic[tr_id]
        assert tr_len == exp_len, "BED transcript length != FASTA transcript length for \"%s\"" %(tr_id)

    # Delete tmp files.
    if os.path.exists(tmp_bed):
        os.remove(tmp_bed)
    if os.path.exists(tmp_fa):
        os.remove(tmp_fa)

    # Return transcript sequences dic constructed from exon sequences.
    return tr_seqs_dic


################################################################################

def concatenate_files(file1, file2):
    """
    Add file 1 content to file 2 (using cat).

    """
    assert os.path.isfile(file1), "file %s does not exist" %(file1)
    assert os.path.isfile(file2), "file %s does not exist" %(file2)
    check_cmd = "cat " + file1 + " >> " + file2
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "cat did not like your input (file 1 %s cat to file 2 %s):\n%s" %(file1, file2, output)


################################################################################

def koma_sepp(n):
    """
    Take input integer n and return comma-separated string, separating
    1000s.

    >>> koma_sepp(131032047)
    '131,032,047'
    >>> koma_sepp(18781)
    '18,781'
    >>> koma_sepp(666)
    '666'

    """
    return '{:,}'.format(n)


################################################################################

def bed_get_transcript_exon_numbers(in_bed):
    """
    Get number of exons for each transcript from in_bed BED file with
    transcript exon regions, with ID format:
    transcriptid_e1 (exon 1), transcriptid_e1 (exon 2)
    This is the output format from gtf_extract_exon_bed(), so both can
    be used in combination.

    >>> in_bed = "test_data/test6.bed"
    >>> bed_get_transcript_exon_numbers(in_bed)
    {'ENST1': 2, 'ENST2': 2, 'ENST3': 1}

    """
    tr_exc_dic = {}
    # Open input .bed file.
    with open(in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            site_id = cols[3]
            if re.search(".+_e\d", site_id):
                m = re.search("(.+)_e\d", site_id)
                tr_id = m.group(1)
                if tr_id not in tr_exc_dic:
                    tr_exc_dic[tr_id] = 1
                else:
                    tr_exc_dic[tr_id] += 1
            else:
                assert False, "site ID \"%s\" missing added _e exon number" %(site_id)
    f.close()
    assert tr_exc_dic, "nothing was read in (\"%s\" empty or malformatted?)" %(in_bed)
    return tr_exc_dic


################################################################################

def bed_get_transcript_lengths_from_exon_regions(in_bed):
    """
    Get spliced transcript lengths from in_bed BED file with transcript
    exon regions, with ID format:
    transcriptid_e1 (exon 1), transcriptid_e1 (exon 2)
    This is the output format from gtf_extract_exon_bed(), so both can
    be used in combination.

    >>> in_bed = "test_data/test6.bed"
    >>> bed_get_transcript_lengths_from_exon_regions(in_bed)
    {'ENST1': 4000, 'ENST2': 1500, 'ENST3': 2500}

    """
    tr_len_dic = {}
    # Open input .bed file.
    with open(in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            site_s = int(cols[1])
            site_e = int(cols[2])
            site_id = cols[3]
            site_len = site_e - site_s
            if re.search(".+_e\d", site_id):
                m = re.search("(.+)_e\d", site_id)
                tr_id = m.group(1)
                if tr_id not in tr_len_dic:
                    tr_len_dic[tr_id] = site_len
                else:
                    tr_len_dic[tr_id] += site_len
            else:
                assert False, "site ID \"%s\" missing added _e exon number" %(site_id)
    f.close()
    assert tr_len_dic, "nothing was read in (\"%s\" empty or malformatted?)" %(in_bed)
    return tr_len_dic


################################################################################

def make_file_copy(in_file, out_file,
                   delete_in=False):
    """
    Make a file copy by copying in_file to out_file.

    """
    check_cmd = "cat " + in_file + " > " + out_file
    assert in_file != out_file, "cat does not like to cat file into same file (%s)" %(check_cmd)
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "cat did not like your input (in_file: %s, out_file: %s):\n%s" %(in_file, out_file, output)
    # Delete in_file.
    if delete_in:
        if os.path.exists(in_file):
            os.remove(in_file)


################################################################################

def bed_read_reg_str_into_dic(in_bed):
    """
    Read BED rows into dictionary of region strings.
    mapping: 'chr1,10,20,+' -> [reg_id].

    >>> test_bed = "test_data/test3.bed"
    >>> bed_read_reg_str_into_dic(test_bed)
    {'chr1,10,20,+': ['CLIP1'], 'chr1,30,45,-': ['CLIP2']}

    """
    reg_str_dic = {}
    with open(in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            chr_id = cols[0]
            site_s = cols[1]
            site_e = cols[2]
            site_id = cols[3]
            site_pol = cols[5]
            reg_str = "%s,%s,%s,%s" %(chr_id, site_s, site_e, site_pol)
            if reg_str in reg_str_dic:
                reg_str_dic[reg_str].append(site_id)
            else:
                reg_str_dic[reg_str] = [site_id]
    assert reg_str_dic, "reg_str_dic empty"
    return reg_str_dic


################################################################################

def bed_read_rows_into_dic(in_bed,
                           new_stem_id="CLIP",
                           max_len=None,
                           sc_thr=None,
                           rev_filter=False,
                           new_ids=False,
                           id2sc_dic=None,
                           id2len_dic=None,
                           to_list=False,
                           id2gen_se_dic=None,
                           check_chr_id_format=True,
                           chr_id_style=1,
                           int_whole_nr=True,
                           remove_id_count=False,
                           chr_ids_dic=None,
                           filt_stats_dic=None):
    """
    Read in .bed file rows into dictionary.
    Mapping is region ID -> bed row.

    id2sc_dic:
        Store column 5 scores for each site ID.
    id2len_dic:
        Store site ID -> site length.
    check_chr_id_format:
        If set check and filter chromosome IDs.
    new_ids:
        If True, generate new IDs.
    to_list:
        Store BED region column values as list.
    remove_id_count:
        If site IDs have format like site_id,count, remove count from ID before
        storing.
    chr_ids_dic:
        Store chromosome IDs in dictionary.


    >>> test_bed = "test_data/test3.bed"
    >>> bed_read_rows_into_dic(test_bed)
    {'CLIP1': 'chr1\\t10\\t20\\tCLIP1\\t0\\t+', 'CLIP2': 'chr1\\t30\\t45\\tCLIP2\\t0\\t-'}

    """

    id2row_dic = {}

    c_read = 0
    c_chr_filt = 0
    c_max_len = 0
    c_sc_thr = 0
    c_out = 0

    with open(in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            chr_id = cols[0]
            site_s = cols[1]
            site_e = cols[2]
            site_id = cols[3]
            site_sc = float(cols[4])
            site_pol = cols[5]
            site_l = int(site_e) - int(site_s)

            if id2gen_se_dic is not None:
                id2gen_se_dic[site_id] = [int(site_s), int(site_e)]

            if remove_id_count:
                m = re.search("(.+),\d", site_id)
                assert m, "remove_id_count True but side ID %s does not contain count" %(site_id)
                site_id = m.group(1)
            # ID count.
            c_read += 1
            # Apply various filters.
            cont = False
            if max_len is not None:
                if site_l > max_len:
                    cont = True
                    c_max_len += 1
            if sc_thr is not None:
                if rev_filter:
                    if site_sc > sc_thr:
                        c_sc_thr += 1
                        cont = True
                else:
                    if site_sc < sc_thr:
                        c_sc_thr += 1
                        cont = True
            if check_chr_id_format:
                new_chr_id = check_convert_chr_id(chr_id,
                                            id_style=chr_id_style)
                if not new_chr_id:
                    c_chr_filt += 1
                    cont = True
                else:
                    chr_id = new_chr_id
            if cont:
                continue

            if not new_ids:
                assert site_id not in id2row_dic, "non-unique site ID (\"%s\") found in \"%s\". Please set --new-site-id or provide unique column 4 --in site IDs" %(site_id, in_bed)
            else:
                if new_stem_id:
                    site_id = new_stem_id + "_" + str(c_read)
                else:
                    site_id = "CLIP_" + str(c_read)

            if id2sc_dic is not None:
                id2sc_dic[site_id] = site_sc
            if id2len_dic is not None:
                id2len_dic[site_id] = site_l

            if int_whole_nr:
                if not site_sc % 1:
                    site_sc = int(site_sc)

            c_out += 1

            if chr_ids_dic is not None:
                chr_ids_dic[chr_id] = 1

            row = "%s\t%s\t%s\t%s\t%s\t%s" %(chr_id, site_s, site_e, site_id, str(site_sc), site_pol)
            if to_list:
                id2row_dic[site_id] = [chr_id, int(site_s), int(site_e), site_id, str(site_sc), site_pol]
            else:
                id2row_dic[site_id] = row

    f.closed

    if filt_stats_dic is not None:
        filt_stats_dic["max_len"] = c_max_len
        filt_stats_dic["sc_thr"] = c_sc_thr
        filt_stats_dic["chr_filt"] = c_chr_filt
        filt_stats_dic["c_in"] = c_read
        filt_stats_dic["c_out"] = c_out

    return id2row_dic


################################################################################

def get_col_count(in_file):
    """
    Count number of columns (separated by TAB) in first row and return
    number.

    >>> in_file = "test_data/test1.bed"
    >>> get_col_count(in_file)
    6
    >>> in_file = "test_data/tr_list.txt"
    >>> get_col_count(in_file)
    1
    >>> in_file = "test_data/empty_file"
    >>> get_col_count(in_file)
    0

    """

    col_c = 0
    with open(in_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            col_c = len(cols)
            break
    f.closed
    return col_c


################################################################################

def bed_write_row_dic_into_file(id2row_dic, out_bed,
                                ext_mode=1,
                                ext_lr=0,
                                id2sc_dic=None,
                                zero_scores=False,
                                id2filt_dic=False,
                                chr_len_dic=False,
                                out_bed_add_sc_dic=False,
                                id2out_dic=False):
    """
    Write .bed row dictionary (column 5 ID as key, .bed row as string)
    into .bed file.
    Example dictionary:
    {'reg1_e1': 'chr1\t1000\t1100\treg1_e1\t0\t+', ... }

    ext_mode:
        1: whole site
        2: center position site
        3: upstream end position site
    ext_lr:
        Extend site (defined by ext_mode) on both sides.
    zero_scores:
        Set to output zero scores (BED column 5), instead of stored scores.
    id2out_dic:
        IDs dictionary for which to output regions.
    id2filt_dic:
        IDs dictionary for which to NOT output regions.
    chr_len_dic:
        Chromosome ID to chromosome length dic for length checking.
    out_bed_add_sc_dic:
        Site ID to score mapping, with score to be added to out_bed ID column
        4, so column 4 ID changes from "id1" to "id1,5"

    """
    assert id2row_dic, "given id2row_dic empty"
    OUTBED = open(out_bed, "w")
    c_out = 0
    for site_id in id2row_dic:
        if id2out_dic:
            if site_id not in id2out_dic:
                continue
        if id2filt_dic:
            if site_id in id2filt_dic:
                continue
        c_out += 1

        cols = id2row_dic[site_id].split("\t")
        chr_id = cols[0]
        reg_s = int(cols[1])
        reg_e = int(cols[2])
        reg_id = cols[3]
        reg_sc = cols[4]
        reg_pol = cols[5]

        # scores to zero needed for twoBitToFa extraction.
        if zero_scores:
            reg_sc = "0"

        # Update region lengths.
        if ext_mode == 1:
            new_s = reg_s
            new_e = reg_e
        elif ext_mode == 2:
            new_e = get_center_position(reg_s, reg_e)
            new_s = new_e - 1
        elif ext_mode == 3:
            new_s = reg_s
            new_e = reg_s + 1
            if reg_pol == "-":
                new_s = reg_e - 1
                new_e = reg_e
        else:
            assert False, "invalid ext_mode set"
        if ext_lr:
            new_s = new_s - ext_lr
            new_e = new_e + ext_lr
            if new_s < 0:
                new_s = 0
            if chr_len_dic:
                assert chr_id in chr_len_dic, "chromosome ID not in chr_len_dic (--gen file)" %(chr_id)
                if new_e > chr_len_dic[chr_id]:
                    new_e = chr_len_dic[chr_id]

        if out_bed_add_sc_dic:
            reg_id += "," + str(out_bed_add_sc_dic[reg_id])

        OUTBED.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(chr_id, str(new_s), str(new_e), reg_id, reg_sc, reg_pol))

    OUTBED.close()
    assert c_out, "nothing was output"


################################################################################

def calc_site_distances(sites_list, id2cp_dic):
    """

    sites_list:
        List with site IDs on reference for which to calculate distances.
    id2cp_dic:
        site_id -> center position (1-based) on reference.

    >>> sites_list = ["sid1", "sid2", "sid3"]
    >>> id2cp_dic = {"sid1": 100, "sid2": 150, "sid3": 300}
    >>> calc_site_distances(sites_list, id2cp_dic)
    {'sid1,sid2': 50, 'sid1,sid3': 200, 'sid2,sid3': 150}

    """
    ids_seen_dic = {}
    sids2dist_dic = {}
    for sid1 in sites_list:
        for sid2 in sites_list:
            if sid1 == sid2:
                continue
            sid1sid2 = sid1 + "," + sid2
            sid2sid1 = sid2 + "," + sid1
            if sid1sid2 in ids_seen_dic:
                continue
            if sid2sid1 in ids_seen_dic:
                continue
            dist = abs(id2cp_dic[sid1] - id2cp_dic[sid2])
            sids2dist_dic[sid1sid2] = dist
            ids_seen_dic[sid1sid2] = 1
            ids_seen_dic[sid2sid1] = 1
    assert sids2dist_dic, "sids2dist_dic empty"
    return sids2dist_dic


################################################################################

def calc_full_site_distances(sites_list, id2coords_dic):
    """
    Calculate distances of start/ends of sites, instead of center positions
    like done in calc_site_distances().

    sites_list:
        List of sitetrids to calculate distances for.

    id2coords_dic:
        site_id,tr_id -> [transcript_s, transcript_e]


    >>> sites_list = ["sid1", "sid2", "sid3", "sid4"]
    >>> id2coords_dic = {"sid1": [80,120], "sid2": [100,110], "sid3": [110,130], "sid4": [150,170]}
    >>> calc_full_site_distances(sites_list, id2coords_dic)
    {'sid1;sid2': 0, 'sid1;sid3': 0, 'sid1;sid4': 30, 'sid2;sid3': 0, 'sid2;sid4': 40, 'sid3;sid4': 20}

    """
    ids_seen_dic = {}
    sids2dist_dic = {}
    for sid1 in sites_list:
        for sid2 in sites_list:
            if sid1 == sid2:
                continue
            sid1sid2 = sid1 + ";" + sid2
            sid2sid1 = sid2 + ";" + sid1
            if sid1sid2 in ids_seen_dic:
                continue
            if sid2sid1 in ids_seen_dic:
                continue

            sid1_s = id2coords_dic[sid1][0]
            sid1_e = id2coords_dic[sid1][1]
            sid2_s = id2coords_dic[sid2][0]
            sid2_e = id2coords_dic[sid2][1]

            dist = get_site_ends_distance(sid1_s, sid1_e, sid2_s, sid2_e)

            sids2dist_dic[sid1sid2] = dist
            ids_seen_dic[sid1sid2] = 1
            ids_seen_dic[sid2sid1] = 1

    assert sids2dist_dic, "sids2dist_dic empty"
    return sids2dist_dic


################################################################################

def get_site_ends_distance(s1, e1, s2, e2):
    """
    Get nearest distance between two site ends / starts.

    >>> get_site_ends_distance(100, 120, 130, 150)
    10
    >>> get_site_ends_distance(130, 150, 100, 120)
    10
    >>> get_site_ends_distance(100, 120, 120, 140)
    0
    >>> get_site_ends_distance(120, 140, 110, 160)
    0

    """
    diff1 = s2 - e1
    diff2 = s1 - e2
    dist = diff2
    if diff1 >= diff2:
        dist = diff1
    if diff1 <= 0 and diff2 <= 0:
        dist = 0
    return dist


################################################################################

def get_center_position(start, end):
    """
    Get center position (1-based), given a (genomic) start (0-based) and
    end coordinate (1-based).

    >>> get_center_position(10, 11)
    11
    >>> get_center_position(1000,2000)
    1501
    >>> get_center_position(11, 20)
    17

    """
    # If region has length of 1, return end position.
    center_pos = end
    # Otherwise calculate new center position.
    if not end - start == 1:
        center_pos = round( ( (end - start) / 2 ) + start ) + 1
    return center_pos


################################################################################

def get_ei_border_ratio_from_exon_id(exon_id, regid2nc_dic,
                                     exid2eibrs_dic=None,
                                     ratio_mode=1,
                                     last_exon_dic=None,
                                     last_exon_ratio=2.5,
                                     min_reg_cov=5,
                                     min_reg_mode=1):
    """

    Ratio is average of ratios at both exon ends (if embedded in introns),
    or if first / last exon, only one ratio.
    Assign -1, if only exon, or if both exon and intron border region read
    count below min_reg_cov.

    min_reg_cov:
        Minimum region read coverage. If both exon and intron border region
        have < min_reg_cov, return ratio of -1.
    regid2nc_dic:
        Contains exon/intron/border region ID -> [norm_cov, coverage, reg_len]
    exid2eibrs_dic:
        Exon ID to all EIB ratios list mapping.
    ratio_mode:
        How to calculate the returned EIBR ratio.
        1: Return the exon-intro border ratio with the higher coverage.
        2: Average the two exon-intron border ratios of the exon,
           if both have more than > min_reg_cov
    last_exon_dic:
        Last transcript exon ID -> polarity
        Used for prioritizing the inner exon intron border for multi-exon
        transcript last exons. Only effective for ratio_mode 1.
    last_exon_ratio:
        If the outer last exon read count is higher last_exon_ratio, prioritize
        the outter border again, i.e. select the outter ratio
        for EIB ratio calculation.

    >>> regid2nc_dic = {"t1_e1_ebe2" : [0.5, 10, 20], "t1_e1_ebi2" : [0.2, 4, 20]}
    >>> get_ei_border_ratio_from_exon_id("t1_e1", regid2nc_dic)
    (2.5, 'first_exon')
    >>> get_ei_border_ratio_from_exon_id("t2_e1", regid2nc_dic)
    (-1, 'single_exon')
    >>> regid2nc_dic = {"t1_e2_ebe1" : [0.5, 10, 20], "t1_e2_ebi1" : [0.25, 5, 20], "t1_e2_ebe2" : [1.0, 20, 20], "t1_e2_ebi2" : [0.25, 5, 20]}
    >>> get_ei_border_ratio_from_exon_id("t1_e2", regid2nc_dic, ratio_mode=2)
    (3.0, 'inner_exon')
    >>> regid2nc_dic = {"t1_e2_ebe1" : [0.5, 10, 20], "t1_e2_ebi1" : [0.25, 5, 20], "t1_e2_ebe2" : [0.1, 2, 20], "t1_e2_ebi2" : [0.1, 2, 20]}
    >>> get_ei_border_ratio_from_exon_id("t1_e2", regid2nc_dic, ratio_mode=2)
    (2.0, 'inner_exon_ds_lc')
    >>> regid2nc_dic = {"t1_e2_ebe1" : [0.1, 2, 20], "t1_e2_ebi1" : [0.1, 2, 20], "t1_e2_ebe2" : [0.5, 10, 20], "t1_e2_ebi2" : [0.25, 5, 20]}
    >>> get_ei_border_ratio_from_exon_id("t1_e2", regid2nc_dic, ratio_mode=2)
    (2.0, 'inner_exon_us_lc')
    >>> regid2nc_dic = {"t1_e1_ebe2" : [0.5, 10, 20], "t1_e1_ebi2" : [0.0, 0, 20]}
    >>> get_ei_border_ratio_from_exon_id("t1_e1", regid2nc_dic)
    (10, 'first_exon')
    >>> regid2nc_dic = {"t1_e1_ebe2" : [0.0, 0, 20], "t1_e1_ebi2" : [0.5, 10, 20]}
    >>> get_ei_border_ratio_from_exon_id("t1_e1", regid2nc_dic)
    (0.0, 'first_exon')
    >>> regid2nc_dic = {"t1_e2_ebe1" : [0.5, 10, 20], "t1_e2_ebi1" : [0.25, 5, 20], "t1_e2_ebe2" : [1.0, 20, 20], "t1_e2_ebi2" : [0.25, 5, 20]}
    >>> get_ei_border_ratio_from_exon_id("t1_e2", regid2nc_dic, ratio_mode=1)
    (4.0, 'inner_exon')


    """

    exb_id_e1 = exon_id + "_ebe1"
    exb_id_i1 = exon_id + "_ebi1"
    exb_id_e2 = exon_id + "_ebe2"
    exb_id_i2 = exon_id + "_ebi2"

    # For single-exon transcripts.
    if exb_id_e1 not in regid2nc_dic and exb_id_e2 not in regid2nc_dic:
        if exid2eibrs_dic is not None:
            assert exon_id not in exid2eibrs_dic, "exon ID %s already stored in exid2eibrs_dic"
            exid2eibrs_dic[exon_id] = [-1]
        return -1, "single_exon"

    # Last exon.
    if exb_id_e1 in regid2nc_dic and exb_id_e2 not in regid2nc_dic:
        assert exb_id_i1 in regid2nc_dic, "exb_id_e1 %s in regid2nc_dic, but not exb_id_i1 %s" %(exb_id_e1, exb_id_i1)
        ratio1 = -1
        sel_crit = "last_exon"
        if regid2nc_dic[exb_id_e1][1] >= min_reg_cov or regid2nc_dic[exb_id_i1][1] >= min_reg_cov:
            if regid2nc_dic[exb_id_i1][0]:
                ratio1 = regid2nc_dic[exb_id_e1][0] / regid2nc_dic[exb_id_i1][0]
            else:
                ratio1 = regid2nc_dic[exb_id_e1][1]
        if exid2eibrs_dic is not None:
            assert exon_id not in exid2eibrs_dic, "exon ID %s already stored in exid2eibrs_dic"
            exid2eibrs_dic[exon_id] = [ratio1]
        return ratio1, sel_crit

    # First exon.
    if exb_id_e1 not in regid2nc_dic and exb_id_e2 in regid2nc_dic:
        assert exb_id_i2 in regid2nc_dic, "exb_id_e2 %s in regid2nc_dic, but not exb_id_i2 %s" %(exb_id_e2, exb_id_i2)
        ratio2 = -1
        sel_crit = "first_exon"
        if regid2nc_dic[exb_id_e2][1] >= min_reg_cov or regid2nc_dic[exb_id_i2][1] >= min_reg_cov:
            if regid2nc_dic[exb_id_i2][0]:
                ratio2 = regid2nc_dic[exb_id_e2][0] / regid2nc_dic[exb_id_i2][0]
            else:
                ratio2 = regid2nc_dic[exb_id_e2][1]
        if exid2eibrs_dic is not None:
            assert exon_id not in exid2eibrs_dic, "exon ID %s already stored in exid2eibrs_dic"
            exid2eibrs_dic[exon_id] = [ratio2]
        return ratio2, sel_crit

    # In-between exons.
    if exb_id_e1 in regid2nc_dic and exb_id_e2 in regid2nc_dic:
        assert exb_id_i1 in regid2nc_dic, "exb_id_e1 %s in regid2nc_dic, but not exb_id_i1 %s" %(exb_id_e1, exb_id_i1)
        assert exb_id_i2 in regid2nc_dic, "exb_id_e2 %s in regid2nc_dic, but not exb_id_i2 %s" %(exb_id_e2, exb_id_i2)
        ratio1 = -1
        ratio2 = -1

        # if exon_id == "ENST00000366553.3_e2":
        #     print(exon_id)
        #     print("regid2nc_dic[exb_id_i1][1]:", regid2nc_dic[exb_id_i1][1])
        #     print("regid2nc_dic[exb_id_e1][1]:", regid2nc_dic[exb_id_e1][1])
        #     print("regid2nc_dic[exb_id_e2][1]:", regid2nc_dic[exb_id_e2][1])
        #     print("regid2nc_dic[exb_id_i2][1]:", regid2nc_dic[exb_id_i2][1])
        #     print("regid2nc_dic[exb_id_i1][0]:", regid2nc_dic[exb_id_i1][0])
        #     print("regid2nc_dic[exb_id_e1][0]:", regid2nc_dic[exb_id_e1][0])
        #     print("regid2nc_dic[exb_id_e2][0]:", regid2nc_dic[exb_id_e2][0])
        #     print("regid2nc_dic[exb_id_i2][0]:", regid2nc_dic[exb_id_i2][0])

        sel_crit = "inner_exon"
        if regid2nc_dic[exb_id_e1][1] >= min_reg_cov or regid2nc_dic[exb_id_i1][1] >= min_reg_cov:
            if regid2nc_dic[exb_id_i1][0]:
                ratio1 = regid2nc_dic[exb_id_e1][0] / regid2nc_dic[exb_id_i1][0]
            else:
                ratio1 = regid2nc_dic[exb_id_e1][1]
        else:
            sel_crit += "_us_lc"
        if regid2nc_dic[exb_id_e2][1] >= min_reg_cov or regid2nc_dic[exb_id_i2][1] >= min_reg_cov:
            if regid2nc_dic[exb_id_i2][0]:
                ratio2 = regid2nc_dic[exb_id_e2][0] / regid2nc_dic[exb_id_i2][0]
            else:
                ratio2 = regid2nc_dic[exb_id_e2][1]
        else:
            sel_crit += "_ds_lc"

        if exid2eibrs_dic is not None:
            assert exon_id not in exid2eibrs_dic, "exon ID %s already stored in exid2eibrs_dic" %(exon_id)
            exid2eibrs_dic[exon_id] = [ratio1, ratio2]

        if ratio1 == -1 and ratio2 != -1:
            avg_ratio = ratio2
        elif ratio1 != -1 and ratio2 == -1:
            avg_ratio = ratio1
        elif ratio1 == -1 and ratio2 == -1:
            avg_ratio = -1
        else:
            if ratio_mode == 1:
                cov_b1 = regid2nc_dic[exb_id_i1][0] + regid2nc_dic[exb_id_e1][0]
                cov_b2 = regid2nc_dic[exb_id_i2][0] + regid2nc_dic[exb_id_e2][0]
                if cov_b1 > cov_b2:
                    avg_ratio = ratio1
                else:
                    avg_ratio = ratio2

                if last_exon_dic is not None:
                    if exon_id in last_exon_dic:
                        sel_crit = "last_exon"
                        exon_pol = last_exon_dic[exon_id]
                        # Define inner borders.
                        cov_inner = cov_b1
                        ratio_inner = ratio1
                        cov_outer = cov_b2
                        ratio_outer = ratio2
                        if exon_pol == "-":
                            cov_inner = cov_b2
                            ratio_inner = ratio2
                            cov_outer = cov_b1
                            ratio_outer = ratio1
                        if cov_inner*last_exon_ratio >= cov_outer:
                            avg_ratio = ratio_inner
                            sel_crit += "_inner"
                        else:
                            avg_ratio = ratio_outer
                            sel_crit += "_outer"

            elif ratio_mode == 2:
                avg_ratio = statistics.mean([ratio1, ratio2])
            else:
                assert False, "invalid ratio_mode (%i)" %(ratio_mode)

        return avg_ratio, sel_crit

    assert False, "invalid get_ei_border_ratio_from_exon_id()"


################################################################################

def get_ei_ratio_from_exon_id(exon_id, regid2nc_dic,
                              isr_intron_reg_dic=False,
                              dummy_int_cov=0.001):
    """
    Given exon ID and coverage information of all exon / intron IDs,
    calculate the exon/intron coverage ratio for the exon.
    Also return the calculated intron coverage.

    regid2nc_dic:
        Exon/Intron ID -> [coverage, overlap read count, region length]
        With coverage = read count / region length.
    dummy_int_cov:
        For intronless transcripts, add a pseudo intron coverage here.
        e.g. 0.001 == 1 read over 1,000 nt

    >>> regid2nc_dic = {"id_e2" : [100.0, 10000, 100], "id_i1" : [2.0, 200, 100], "id_i2" : [6.0, 600, 100], "id2_e1" : [100.0, 10000, 100]}
    >>> get_ei_ratio_from_exon_id("id_e2", regid2nc_dic)
    (25.0, 4.0, 800)
    >>> get_ei_ratio_from_exon_id("id2_e1", regid2nc_dic)
    (100000.0, 0.001, 0)
    >>> regid2nc_dic = {"id_e2" : [100.0, 10000, 100], "id_i1" : [4.0, 400, 100], "id_i2" : [0.0, 0, 100]}
    >>> get_ei_ratio_from_exon_id("id_e2", regid2nc_dic)
    (50.0, 2.0, 400)
    >>> regid2nc_dic = {"id_e2" : [100.0, 10000, 100], "id_i1" : [0.0, 0, 100], "id_i2" : [0.0, 0, 100]}
    >>> get_ei_ratio_from_exon_id("id_e2", regid2nc_dic)
    (100.0, 0.0, 0)

    """
    assert re.search(".+_e\d", exon_id), "exon ID %s has invalid format"
    m = re.search("(.+)_e(\d+)", exon_id)
    tr_id = m.group(1)
    ex_nr = int(m.group(2))
    ex_cov = regid2nc_dic[exon_id][0]
    us_intron_id = tr_id + "_i" + str(ex_nr-1)
    ds_intron_id = tr_id + "_i" + str(ex_nr)
    int_count = 0
    int_cov = 0.0
    div = 0
    if us_intron_id in regid2nc_dic:
        int_cov += regid2nc_dic[us_intron_id][0]
        int_count += regid2nc_dic[us_intron_id][1]
        div += 1
    if ds_intron_id in regid2nc_dic:
        int_cov += regid2nc_dic[ds_intron_id][0]
        int_count += regid2nc_dic[ds_intron_id][1]
        div += 1
    if div:
        int_cov = int_cov / div
    else:
        # For intronless transcripts (exon has no neighboring introns).
        int_cov = dummy_int_cov
    if int_cov == 0.0:
        ei_ratio = ex_cov
    else:
        ei_ratio = ex_cov / int_cov
    return ei_ratio, int_cov, int_count


################################################################################

def check_transcript_eir_tc_state(tr_id, exid2cov_dic, trid2exc_dic,
                                  filter_ei_ratio=1.0,
                                  min_ei_cov=20,
                                  min_ei_cov_sum=False,
                                  exid2isrn_dic=False,
                                  max_isrn_c=0,
                                  min_occ_exons_ratio=0.2):
    """
    Check whether a given transcript (tr_id) should be assigned transcript
    context (return True), or genomic context (return False). Base this
    decision on EIR ratios present in the transcript. Further controlled
    by parameters:

    filter_ei_ratio:
        Exons with EIR <= filter_ei_ratio are counted for tr_id.
    min_occ_exons_ratio:
        Set ratio of the total number of transcript exons required for
        transcript to be assigned to genomic context (return False).
        So if # of exons <= filter_ei_ratio is >= then the ratio of
        exons (round() to next integer, or math.ceil), return False.
    min_ei_cov:
        Minimum coverage (== read count) the exon or surrounding intron
        region(s) need to have to be included into exon counting.
    min_ei_cov_sum:
        Instead of OR, use the sum of exon and intron region for min_ei_cov.
    exid2isrn_dic:
        Exon ID -> ISR count to neighborhood.
        If this is given, use only exons with ISRN counts <= max_isrn_c.
    max_isrn_c:
        Define maximum ISRN count (exid2isrn_dic) (default: 0).
    exid2cov_dic:
        exon ID -> [ei_ratio, ex_cov, int_cov, ex_read_count, int_read_count]
    trid2exc_dic:
        Transcript ID -> exon count

    >>> exid2cov_dic = {'t1_e1': [1.0, 1, 1, 15, 25], 't1_e2': [2.0, 2, 1, 35, 25], 't1_e3': [1.0, 2, 2, 35, 45]}
    >>> trid2exc_dic = {'t1': 3}
    >>> check_transcript_eir_tc_state("t1", exid2cov_dic, trid2exc_dic)
    False
    >>> check_transcript_eir_tc_state("t1", exid2cov_dic, trid2exc_dic, min_ei_cov=50)
    True

    """

    assert min_occ_exons_ratio > 0 and min_occ_exons_ratio < 1, "min_occ_exons_ratio needs to be > 0 and < 1"

    tr_exc = trid2exc_dic[tr_id]
    if tr_exc == 1:
        return True

    c_gc_eir = 0  # count exons with EIR <= filter_ei_ratio.
    for i in range(tr_exc):
        ex_nr = i + 1
        ex_id = tr_id + "_e" + str(ex_nr)
        if exid2isrn_dic:
            isrn_c = exid2isrn_dic[ex_id]
            if isrn_c > max_isrn_c:
                continue
        ei_ratio = exid2cov_dic[ex_id][0]
        ex_rc = exid2cov_dic[ex_id][3]
        int_rc = exid2cov_dic[ex_id][4]
        # print("ex_id:", ex_id, "ei_ratio:", ei_ratio)
        # print("ex_rc:", ex_rc, "int_rc:", int_rc)
        check_count = max(ex_rc, int_rc)
        if min_ei_cov_sum:
            check_count = ex_rc + int_rc
        if check_count >= min_ei_cov:
            if ei_ratio <= filter_ei_ratio:
                c_gc_eir += 1

    min_occ_exons_set = math.ceil(tr_exc*min_occ_exons_ratio)
    if min_occ_exons_set < 1:
        min_occ_exons_set = 1

    if c_gc_eir >= min_occ_exons_set:
        return False
    else:
        return True


################################################################################

def check_transcript_eibr_tc_state(tr_id, exid2eibr_dic, trid2exc_dic,
                                  filter_eib_ratio=2.0,
                                  min_occ_exons_ratio=0.2):
    """
    Check whether a given transcript (tr_id) should be assigned transcript
    context (return True), or genomic context (return False). Base this
    decision on EIR ratios present in the transcript. Further controlled
    by parameters:

    filter_ei_ratio:
        Exons with EIBR <= filter_eib_ratio are counted as genomic context
        evidence for tr_id. So the higher filter_ei_ratio is set, the more
        likely tr_id gets assigned to genomic context (returns False).

    min_occ_exons_ratio:
        Set ratio of the total number of transcript exons required for
        transcript to be assigned to genomic context (return False).
        So if # of exons <= filter_eib_ratio is >= then the ratio of
        exons (round() to next integer, or math.ceil), return False.
    min_occ_exons:
        If there are n exons EIBR <= filter_eib_ratio, and n >= min_occ_exons,
        report this transcript as genomic context (return False).
        If min_occ_exons_ratio is set (e.g. 0.25), then use the max of
        min_occ_exons and round(min_occ_exons_ratio*number_of_exons)
    min_occ_exons_ratio:
        Set ratio of total number transcript exons for min_occ_exons.
    exid2eibr_dic:
        exon ID -> exon-intron border ratio
    trid2exc_dic:
        Transcript ID -> exon count

    >>> exid2eibr_dic = {'t1_e1': -1, 't1_e2': 4, 't1_e3': 2, 't1_e4': 1.75}
    >>> trid2exc_dic = {'t1': 4}
    >>> check_transcript_eibr_tc_state("t1", exid2eibr_dic, trid2exc_dic)
    False
    >>> check_transcript_eibr_tc_state("t1", exid2eibr_dic, trid2exc_dic, filter_eib_ratio=1.5)
    True
    >>> check_transcript_eibr_tc_state("t1", exid2eibr_dic, trid2exc_dic, min_occ_exons_ratio=0.75)
    True

    """

    assert min_occ_exons_ratio > 0 and min_occ_exons_ratio < 1, "min_occ_exons_ratio needs to be > 0 and < 1"

    tr_exc = trid2exc_dic[tr_id]
    if tr_exc == 1:
        return True

    c_gc_eibr = 0  # count exons with EIBR <= filter_eib_ratio.
    for i in range(tr_exc):
        ex_nr = i + 1
        ex_id = tr_id + "_e" + str(ex_nr)
        ex_eibr = exid2eibr_dic[ex_id]
        #print("ex_id:", ex_id, "ex_eibr:", ex_eibr)
        if ex_eibr == -1:
            continue
        if ex_eibr <= filter_eib_ratio:
            c_gc_eibr += 1

    #print("c_gc_eibr:", c_gc_eibr)

    min_occ_exons_set = round(tr_exc*min_occ_exons_ratio)
    if min_occ_exons_set < 1:
        min_occ_exons_set = 1

    if c_gc_eibr >= min_occ_exons_set:
        return False
    else:
        return True


################################################################################

def select_highest_conf_tr_id(tr_ids_list,
                              trid2isrc_dic, trid2tslsc_dic,
                              trid2nc_dic, trid2len_dic,
                              sel_mode=1):
    """
    Select transcript ID with highest confidence for an exonic --in site.

    sel_mode:
        Selection mode
        1: # IS reads > transcript coverage > TSL > transcript length > transcript ID
        2: # IS reads > TSL > transcript coverage > transcript length > transcript ID

    Priority (sel_mode 1):
    1) trid2isrc_dic
    2) trid2nc_dic
    3) trid2tslsc_dic
    4) trid2len_dic
    5) transcript ID

    Priority (sel_mode 2):
    1) trid2isrc_dic
    2) trid2tslsc_dic
    3) trid2nc_dic
    4) trid2len_dic
    5) transcript ID

    tr_ids_list:
        List of transcript IDs, from which to determine transcript ID
        with highest confidence.
    trid2isrc_dic:
        transcript ID to IS reads count supporting the transcript.
    trid2tslsc_dic:
        transcript ID to TSL score.
    trid2nc_dic:
        transcript ID to normalized coverage (sum of exon coverages).
        Normalized coverage: # read starts / region length.
    trid2len_dic:
        transcript ID to transcript length.

    If there is a draw after 4 rounds, select transcript ID with lexicographic
    filtering ID and choosing first one.

    """
    # The chosen one.
    tco_id = ""
    tco_isrc = -666
    tco_tslsc = 0
    tco_nc = 0.0
    tco_len = 0
    # Select based on IS read counts.
    sel_ids = get_highest_scoring_ids(tr_ids_list, trid2isrc_dic)
    if len(sel_ids) == 1:
        return sel_ids[0], "is_read_c"
    if sel_mode == 1:
        # Select based on normalized coverage.
        sel_ids = get_highest_scoring_ids(sel_ids, trid2nc_dic)
        if len(sel_ids) == 1:
            return sel_ids[0], "norm_cov"
        # Select based on TSL scores.
        sel_ids = get_highest_scoring_ids(sel_ids, trid2tslsc_dic)
        if len(sel_ids) == 1:
            return sel_ids[0], "tsl_sc"
    elif sel_mode == 2:
        # Select based on TSL scores.
        sel_ids = get_highest_scoring_ids(sel_ids, trid2tslsc_dic)
        if len(sel_ids) == 1:
            return sel_ids[0], "tsl_sc"
        # Select based on normalized coverage.
        sel_ids = get_highest_scoring_ids(sel_ids, trid2nc_dic)
        if len(sel_ids) == 1:
            return sel_ids[0], "norm_cov"
    else:
        assert False, "invalid sel_mode set"
    # Select based on transcript length.
    sel_ids = get_highest_scoring_ids(sel_ids, trid2len_dic)
    if len(sel_ids) == 1:
        return sel_ids[0], "tr_len"
    # Select based on sorting IDs.
    sel_ids.sort()
    return sel_ids[0], "id_sort"


################################################################################

def get_highest_scoring_ids(ids_set, id2sc_dic,
                            scfp_list=False):
    """
    Given a set of IDs, and a ID to score mapping, return highest
    scoring ID(s) in list.

    scfp_list:
        If score is stored in first position of list, with ID
        mapping to list in id2sc_dic.

    >>> ids_set = ['s1', 's2', 's3',  's4']
    >>> id2sc_dic = {'s1' : 10, 's2' : 5, 's3' : 10, 's4' : 7}
    >>> get_highest_scoring_ids(ids_set, id2sc_dic)
    ['s1', 's3']
    >>> id2sc_dic = {'s1' : 10, 's2' : 5, 's3' : 10, 's4' : 17}
    >>> get_highest_scoring_ids(ids_set, id2sc_dic)
    ['s4']
    >>> ids_set = ['s1', 's2', 's3']
    >>> id2sc_dic = {'s1' : 0, 's2' : 0, 's3' : 0}
    >>> get_highest_scoring_ids(ids_set, id2sc_dic)
    ['s1', 's2', 's3']

    """
    max_ids = []
    max_sc = -6666
    for set_id in ids_set:
        if scfp_list:
            if max_sc < id2sc_dic[set_id][0]:
                max_sc = id2sc_dic[set_id][0]
        else:
            if max_sc < id2sc_dic[set_id]:
                max_sc = id2sc_dic[set_id]
    for set_id in ids_set:
        if scfp_list:
            if id2sc_dic[set_id][0] == max_sc:
                max_ids.append(set_id)
        else:
            if id2sc_dic[set_id] == max_sc:
                max_ids.append(set_id)
    assert max_ids, "max_ids empty"
    max_ids = list(set(max_ids))
    max_ids.sort()
    return max_ids


################################################################################

def check_a_in_b_list(a_list, b_list):
    """
    Check if any entry in list A is present in list B. If true, return True,
    unless False.

    """
    check = False
    for a in a_list:
        if a in b_list:
            check = True
            break
    return check


################################################################################

def two_lists_get_intersect(l1,l2):
    """
    Given two lists, return list with common elements.

    >>> l1 = [1,5,10,20,30]
    >>> l2 = [5,20,40]
    >>> two_lists_get_intersect(l1, l2)
    [5, 20]
    >>> l3 = [50]
    >>> two_lists_get_intersect(l1, l3)
    []

    """

    assert l1, "given l1 empty"
    assert l2, "given l2 empty"
    l3 = [element for element in l1 if element in l2]
    return l3


################################################################################

def bed_check_unique_ids(bed_file):
    """
    Check whether .bed file (6 column format with IDs in column 4)
    has unique column 4 IDs.

    >>> test_bed = "test_data/test1.bed"
    >>> bed_check_unique_ids(test_bed)
    True
    >>> test_bed = "test_data/test2.bed"
    >>> bed_check_unique_ids(test_bed)
    False

    """

    check_cmd = "cut -f 4 " + bed_file + " | sort | uniq -d"
    output = subprocess.getoutput(check_cmd)
    if output:
        return False
    else:
        return True


################################################################################

def diff_two_files_identical(file1, file2):
    """
    Check whether two files are identical. Return true if diff reports no
    differences.

    >>> file1 = "test_data/file1"
    >>> file2 = "test_data/file2"
    >>> diff_two_files_identical(file1, file2)
    True
    >>> file1 = "test_data/test1.bed"
    >>> diff_two_files_identical(file1, file2)
    False

    """
    same = True
    check_cmd = "diff " + file1 + " " + file2
    output = subprocess.getoutput(check_cmd)
    if output:
        same = False
    return same


################################################################################

def gtf_get_intron_exon_cov(in_gtf, in_bam, out_bed,
                            correct_min_ex_order=False,
                            tr2exc_dic=False,
                            read_pos_mode=1,
                            eib_width=10,
                            border_mode=1,
                            count_isr_double=True,
                            chr_id_style=1,
                            add_isr_bed=False,
                            reg2cov_dic=None,
                            isr_sub_count=True,
                            isr_intron_reg_dic=False,
                            tmp_out_folder=False,
                            tr_ids_dic=False):
    """

    Get intron-exon region coverage of exonic site transcripts.

    read_pos_mode:
        Defines what part of read to use for coverage calculation.
        1: full-length read (thus can be counted > 1)
        2: 5' end of read
        3: center position of read
    border_mode:
        Mode for exon-intron border extration.
        1: on each side of each exon
        2: on first and last exon only the inner sites, where
           exon is connected to intron.
    add_isr_bed:
        Provide ISR end position BED from earlier computations, to
        count intron-spanning reads twice for coverage calculations.
        Does not include --rnaseq-bam reads!
    reg2cov_dic:
        ["chr,s,e,pol"] -> read overlap count / coverage.
    regid2reg_dic:
        Region ID (exon intron) to genomic region "chr_id,s,e,pol"
    tr_ids_dic:
        Transcript IDs to keep dictionary.

    Output intron/exon regions of specified transcripts from GTF to BED.
    E.g.
    chr1	1000	2000	ENST001_e1	0	+
    chr1	2000	3000	ENST001_i1	0	+
    chr1	3000	4000	ENST001_e2	0	+
    chr1	6000	7000	ENST002_e2	0	-
    chr1	7000	8000	ENST002_i1	0	-
    chr1	8000	9000	ENST002_e1	0	-

    Then convert in_bam to BED regions, overlap with intron/exon regions
    and get overlap counts normalized by intron/exon lengths.
    Then overlap in_bam BAM reads with intron/exon regions,
    counting overlaps for intron / exon regions.

    bamToBed -i TAF13_gene_reads.bam -split

    """

    # For reverse ording we need to have exon numbers.
    if correct_min_ex_order:
        assert tr2exc_dic, "tr2exc_dic needed if correct_min_ex_order True"

    random_id = uuid.uuid1()
    tmp_bed = str(random_id) + ".intron_exon.tmp.bed"
    if tmp_out_folder:
        tmp_bed = tmp_out_folder + "/" + tmp_bed
    OUTBED = open(tmp_bed, "w")

    # Read in exon features from GTF file.
    c_gtf_ex_feat = 0
    # Start end coordinates of exons.
    exon_e_dic = {}
    exon_s_dic = {}
    # Transcript stats.
    tr2pol_dic = {}
    tr2chr_dic = {}
    # dic for sanity checking exon number order.
    tr2exon_nr_dic = {}
    # Remember me.
    proc_tr_dic = {}

    # Open GTF either as .gz or as text file.
    if re.search(".+\.gz$", in_gtf):
        f = gzip.open(in_gtf, 'rt')
    else:
        f = open(in_gtf, "r")
    for line in f:
        # Skip header.
        if re.search("^#", line):
            continue
        cols = line.strip().split("\t")
        chr_id = cols[0]
        feature = cols[2]
        feat_s = int(cols[3])
        feat_e = int(cols[4])
        feat_pol = cols[6]
        infos = cols[8]

        if feature != "exon":
            continue

        # Restrict to standard chromosomes.
        new_chr_id = check_convert_chr_id(chr_id,
                                id_style=chr_id_style)
        if not new_chr_id:
            continue
        else:
            chr_id = new_chr_id

        # Make start coordinate 0-base (BED standard).
        feat_s = feat_s - 1

        # Extract transcript ID.
        m = re.search('transcript_id "(.+?)"', infos)
        assert m, "transcript_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        transcript_id = m.group(1)

        # Check if transcript ID is in transcript dic.
        if tr_ids_dic:
            if transcript_id not in tr_ids_dic:
                if len(tr_ids_dic) == len(proc_tr_dic):
                    break
                else:
                    continue
        proc_tr_dic[transcript_id] = 1

        # Extract exon number.
        m = re.search('exon_number "(\d+?)"', infos)
        # Try Gencode encoding.
        if not m:
            m = re.search('exon_number (\d+?);', infos)
        assert m, "exon_number entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        exon_nr = int(m.group(1))

        # Store transcript stats.
        tr2pol_dic[transcript_id] = feat_pol
        tr2chr_dic[transcript_id] = chr_id

        # Check whether exon numbers are incrementing for each transcript ID.
        if not transcript_id in tr2exon_nr_dic:
            tr2exon_nr_dic[transcript_id] = exon_nr
        else:
            assert tr2exon_nr_dic[transcript_id] < exon_nr, "transcript ID \"%s\" without increasing exon number order in GTF file \"%s\"" %(transcript_id, in_gtf)
            tr2exon_nr_dic[transcript_id] = exon_nr

        # Reverse ordering for minus strand.
        if correct_min_ex_order and feat_pol == "-":
            assert transcript_id in tr2exc_dic, "transcript ID %s not in tr2exc_dic" %(transcript_id)
            exon_nr = tr2exc_dic[transcript_id] - exon_nr + 1
            assert exon_nr >= 1, "exon number < 1 assigned (%i) for transcript ID %s (# exons: %i)" %(exon_nr, transcript_id, tr2exc_dic[transcript_id])
        # Construct exon ID.
        exon_id = transcript_id + "_e" + str(exon_nr)

        # Store infos.
        exon_s_dic[exon_id] = feat_s
        exon_e_dic[exon_id] = feat_e

        # Count exon entry.
        c_gtf_ex_feat += 1

        # Output genomic exon region.
        OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,feat_s,feat_e,exon_id,feat_pol))

    OUTBED.close()
    f.close()

    # Check for read-in features.
    assert c_gtf_ex_feat, "no exon features read in from \"%s\"" %(in_gtf)

    # Append intron regions.
    tr2intron_nr_dic = {}
    OUTBED = open(tmp_bed, "a")
    for tr_id in tr2pol_dic:
        tr_pol = tr2pol_dic[tr_id]
        chr_id = tr2chr_dic[tr_id]
        tr_c = tr2exon_nr_dic[tr_id]
        tr2intron_nr_dic[tr_id] = 0
        intron_c = 0

        # # 1-exon transcripts, no introns.
        # if tr_c == 1:
        #     continue

        ex_list = []
        for i in range(tr_c):
            ex_nr = i + 1
            ex_id = tr_id + "_e" + str(ex_nr)
            ex_list.append(ex_id)

        # If one exon only.
        if tr_c == 1:
            if border_mode == 1:
                ex_id = ex_list[0]
                ex_s = exon_s_dic[ex_id]
                ex_e = exon_e_dic[ex_id]
                exb_id_i1 = ex_id + "_ebi1"
                exb_id_e1 = ex_id + "_ebe1"
                exb_id_e2 = ex_id + "_ebe2"
                exb_id_i2 = ex_id + "_ebi2"

                i1s, i1e, e1s, e1e, e2s, e2e, i2s, i2e = get_intron_exon_border_coords(
                                                            ex_s, ex_e, eib_width=eib_width)

                OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,i1s,i1e,exb_id_i1,tr_pol))
                OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,e1s,e1e,exb_id_e1,tr_pol))
                OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,e2s,e2e,exb_id_e2,tr_pol))
                OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,i2s,i2e,exb_id_i2,tr_pol))

            continue

        # For multi-exon transcripts, output intronic and border regions.
        intron_len_list = []

        for i in range(len(ex_list)):

            ex1i = i
            ex2i = i + 1

            # As long as not last exon, add introns.
            if ex2i < len(ex_list):

                ex1id = ex_list[ex1i]
                ex2id = ex_list[ex2i]

                ex1s = exon_s_dic[ex1id]
                ex2s = exon_s_dic[ex2id]
                ex1e = exon_e_dic[ex1id]
                ex2e = exon_e_dic[ex2id]

                intron_id = tr_id + "_i" + str(ex2i)
                intron_c += 1

                # Plus case.
                intron_s = ex1e
                intron_e = ex2s
                ex1_len = ex1e - ex1s
                ex2_len = ex2e - ex2s
                # Lengths of exons and embedded intron.
                triple_len_list = [ex1_len, 0, ex2_len]
                triple_ids_list = [ex1id, intron_id, ex2id]
                if tr_pol == "-":
                    intron_s = ex2e
                    intron_e = ex1s
                    triple_len_list[0] = ex2_len
                    triple_len_list[2] = ex1_len

                prev_intron_len = False
                if intron_len_list:
                    ill_len = len(intron_len_list)
                    prev_intron_len = intron_len_list[ill_len-1]

                intron_len = intron_e - intron_s
                intron_len_list.append(intron_len)
                triple_len_list[1] = intron_len

                # Output intron region.
                OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,intron_s,intron_e,intron_id,tr_pol))

                """
                Get intron-exon border regions.
                I.e. regions at intron / exon ends, for which to calculate
                coverage too.
                """

                i1s, i1e, e1s, e1e, e2s, e2e, i2s, i2e = get_intron_exon_border_coords(
                                                            ex1s, ex1e, eib_width=eib_width,
                                                            us_intron_len=prev_intron_len,
                                                            ds_intron_len=intron_len)

                exb_id_i1 = ex1id + "_ebi1"
                exb_id_e1 = ex1id + "_ebe1"
                exb_id_e2 = ex1id + "_ebe2"
                exb_id_i2 = ex1id + "_ebi2"

                if border_mode == 1:
                    OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,i1s,i1e,exb_id_i1,tr_pol))
                    OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,e1s,e1e,exb_id_e1,tr_pol))
                    OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,e2s,e2e,exb_id_e2,tr_pol))
                    OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,i2s,i2e,exb_id_i2,tr_pol))
                elif border_mode == 2:
                    if prev_intron_len:
                        OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,i1s,i1e,exb_id_i1,tr_pol))
                        OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,e1s,e1e,exb_id_e1,tr_pol))
                    OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,e2s,e2e,exb_id_e2,tr_pol))
                    OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,i2s,i2e,exb_id_i2,tr_pol))
                else:
                    assert False, "invalid border_mode set (%i)" %(border_mode)

            elif ex2i == len(ex_list):

                # Last exon.
                ex_id = ex_list[ex1i]
                ex_s = exon_s_dic[ex_id]
                ex_e = exon_e_dic[ex_id]

                # Previous intron length.
                prev_intron_len = intron_len_list[intron_c-1]

                i1s, i1e, e1s, e1e, e2s, e2e, i2s, i2e = get_intron_exon_border_coords(
                                                            ex_s, ex_e, eib_width=eib_width,
                                                            us_intron_len=prev_intron_len,
                                                            ds_intron_len=False)

                exb_id_i1 = ex_id + "_ebi1"
                exb_id_e1 = ex_id + "_ebe1"
                exb_id_e2 = ex_id + "_ebe2"
                exb_id_i2 = ex_id + "_ebi2"

                if border_mode == 1:
                    OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,i1s,i1e,exb_id_i1,tr_pol))
                    OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,e1s,e1e,exb_id_e1,tr_pol))
                    OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,e2s,e2e,exb_id_e2,tr_pol))
                    OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,i2s,i2e,exb_id_i2,tr_pol))
                elif border_mode == 2:
                    OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,i1s,i1e,exb_id_i1,tr_pol))
                    OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,e1s,e1e,exb_id_e1,tr_pol))
                else:
                    assert False, "invalid border_mode set (%i)" %(border_mode)

                break

        tr2intron_nr_dic[tr_id] = intron_c

    OUTBED.close()

    # Sanity check exon + intron numbers.
    for tr_id in tr2exon_nr_dic:
        exon_nr = tr2exon_nr_dic[tr_id]
        intron_nr = tr2intron_nr_dic[tr_id]
        assert (exon_nr-1) == intron_nr, "intron number != exon number - 1 for \"%s\" (%i != %i - 1)" %(tr_id, intron_nr, exon_nr)

    # Sort intron exon BED.
    bed_sort_file(tmp_bed, out_bed)

    # First get start positions for overlap calculation.
    # Minus strand: end, plus strand: start.
    check_cmd = "bamToBed -i " + in_bam + " -split"
    output = subprocess.getoutput(check_cmd)

    GRRSBEDOUT = open(tmp_bed, "w")
    c_read = 0
    for line in output.split('\n'):
        cols = line.strip().split("\t")
        c_read += 1
        chr_id = cols[0]
        reg_s = int(cols[1])
        reg_e = int(cols[2])
        read_pol = cols[5]
        new_id = "r%i" %(c_read)
        if read_pos_mode == 1:
            new_s = reg_s
            new_e = reg_e
        elif read_pos_mode == 2:
            new_s = reg_s
            new_e = new_s + 1
            if read_pol == "-":
                new_s = reg_e - 1
                new_e = reg_e
        elif read_pos_mode == 3:
            new_e = get_center_position(reg_s, reg_e)
            new_s = new_e - 1
        else:
            assert False, "invalid value set for read_pos_mode"
        assert new_e > new_s, "BED region end <= region start coordinate"
        GRRSBEDOUT.write("%s\t%i\t%i\t%s\t0\t%s\n" %(chr_id, new_s, new_e, new_id, read_pol))
    GRRSBEDOUT.close()
    assert c_read, "no output produced from \"%s\"" %(check_cmd)


    """
    Append ISR BED to reads BED, to count ISR reads twice in overlap
    calculations.

    """
    if add_isr_bed and count_isr_double:
        assert os.path.exists(add_isr_bed), "set ISR containing BED file %s does not exist" %(add_isr_bed)
        print("Concatenate ISR reads BED to transcript reads BED ... ")
        concatenate_files(add_isr_bed, tmp_bed)

    # Region (intron exon border) to coverage stats.
    regid2nc_dic = {}

    """
    Example output:

    $ cat a.bed
    chr1	1000	1030	t1	0	+
    chr1	1050	1080	t2	0	+
    chr1	1100	1130	t3	0	+
    $ cat b.bed
    chr1	1000	1001	r1	0	+
    chr1	1010	1011	r2	0	+
    chr1	1020	1021	r3	0	+
    chr1	1050	1051	r4	0	+
    chr1	1060	1061	r5	0	+
    $ intersectBed -a a.bed -b b.bed -s -sorted -c
    chr1	1000	1030	t1	0	+	3
    chr1	1050	1080	t2	0	+	2
    chr1	1100	1130	t3	0	+	0

    NOTE that regions without overlaps also appear here ...

    And with full-length regions:
    a.bed
    chr1	1000	2000	e1	0	+
    chr1	2000	3000	i1	0	+
    chr1	3000	4000	e2	0	+
    b.bed
    chr1	1990	2020	r1	0	+
    chr1	1995	2025	r2	0	+
    chr1	2000	2040	r3	0	+
    intersectBed -a a.bed -b b.bed -s -c
    chr1	1000	2000	e1	0	+	2
    chr1	2000	3000	i1	0	+	3
    chr1	3000	4000	e2	0	+	0

    """

    check_cmd = "sort -k1,1 -k2,2n " + tmp_bed + " | intersectBed -a " + out_bed + " -b stdin -s -sorted -c"
    output = subprocess.getoutput(check_cmd)
    for line in output.split('\n'):
        cols = line.strip().split("\t")
        chr_id = cols[0]
        reg_s = int(cols[1])
        reg_e = int(cols[2])
        reg_id = cols[3]
        reg_pol = cols[5]
        ol_c = int(cols[6])
        reg_l = reg_e - reg_s

        gen_reg = "%s,%i,%i,%s" %(chr_id, reg_s, reg_e, reg_pol)
        if reg2cov_dic is not None:
            reg2cov_dic[gen_reg] = ol_c

        # Add pseudocounts only for exon / intron regions, for coverage calculations.
        if not re.search(".+_eb[ei]\d+$", reg_id):
            ol_c += 1

        # Substract ISR count from intronic read count (if isr_sub_count).
        if isr_intron_reg_dic:
            if gen_reg in isr_intron_reg_dic:
                isr_c = isr_intron_reg_dic[gen_reg]
                if isr_sub_count:
                    ol_c -= isr_c
                    if ol_c < 1:
                        ol_c = 1

        norm_c = ol_c / reg_l
        regid2nc_dic[reg_id] = [norm_c, ol_c, reg_l]

    assert regid2nc_dic, "regid2nc_dic empty (nothing read in)"

    if os.path.exists(tmp_bed):
        os.remove(tmp_bed)

    return regid2nc_dic


################################################################################

def get_intron_exon_border_coords(ex_s, ex_e,
                                  us_intron_len=False,
                                  ds_intron_len=False,
                                  eib_width=20):
    """
    Get intron-exon border region coordinates surrounding one exon.

    ex_s:
        Genomic exon start (0-based)
    ex_e:
        Genomic exon end (1-based)
    us_intron_len:
        Upstream intron length (for single exon transcripts = False)
    ds_intron_len:
        Downstream intron length (for single exon transcripts = False)
    eib_width:
        Width of intron/exon border region.

    >>> get_intron_exon_border_coords(1000, 2000)
    (980, 1000, 1000, 1020, 1980, 2000, 2000, 2020)
    >>> get_intron_exon_border_coords(1000, 1020, eib_width=50)
    (950, 1000, 1000, 1020, 1000, 1020, 1020, 1070)
    >>> get_intron_exon_border_coords(1000, 1020, eib_width=50, us_intron_len=30, ds_intron_len=40)
    (970, 1000, 1000, 1020, 1000, 1020, 1020, 1060)

    """
    # Exon length.
    ex_len = ex_e - ex_s

    # Upstream intron border region.
    i1s = ex_s - eib_width
    if us_intron_len:
        if us_intron_len < eib_width:
            i1s = ex_s - us_intron_len
    i1e = ex_s
    # Upstream exon border region.
    e1s = ex_s
    e1e = ex_s + eib_width
    if ex_len < eib_width:
        e1e = ex_s + ex_len

    # Downstream exon border region.
    e2s = ex_e - eib_width
    if ex_len < eib_width:
        e2s = ex_e - ex_len
    e2e = ex_e

    # Downstream intron border region.
    i2s = ex_e
    i2e = ex_e + eib_width
    if ds_intron_len:
        if ds_intron_len < eib_width:
            i2e = ex_e + ds_intron_len

    return i1s, i1e, e1s, e1e, e2s, e2e, i2s, i2e


################################################################################

def get_gene_infos_from_annot_table(tr_gene_annot_file,
                                    trid2gid_dic, trid2gn_dic, trid2gbt_dic,
                                    trid2tbt_dic):
    """
    Read in gene infos from peakhood extract folder file.
    Store infos with mapping: transcript ID -> gene info

    tr_gene_annot_file content:
    transcript_id	tr_biotype	tr_length	tr_exon_c	tr_gene_id	tr_gene_name	tr_gene_biotype
    ENST00000379370.7	protein_coding	7328	36	MSTRG.88	AGRN	-
    ENST00000620552.4	protein_coding	7394	39	MSTRG.88	AGRN	-
    ...

    """

    with open(tr_gene_annot_file) as f:
        for line in f:
            row = line.strip()
            cols = line.strip().split("\t")
            if cols[0] == "transcript_id":
                continue
            tr_id = cols[0]
            tr_bt = cols[1]
            gene_id = cols[4]
            gene_name = cols[5]
            gene_bt = cols[6]
            trid2gid_dic[tr_id] = gene_id
            trid2gn_dic[tr_id] = gene_name
            trid2gbt_dic[tr_id] = gene_bt
            trid2tbt_dic[tr_id] = tr_bt
    f.close()


################################################################################

def get_gen_motif_hits(motif, gen_fa, gen_bed,
                       hits_out_bed=False,
                       hits_out_fa=False):
    """
    Get motif hits on transcript sequences / sites.

    Return number of unique hits and effective region size.

    """
    c_hits = 0
    c_uniq_hits = 0
    c_sites = 0
    region_size = 0

    gen_seqs_dic = read_fasta_into_dic(gen_fa,
                                       dna=False,
                                       empty_check=False,
                                       skip_n_seqs=False)
    if not gen_seqs_dic:
        return c_hits, c_uniq_hits, c_sites, region_size

    gen_bed_dic = bed_read_rows_into_dic(gen_bed, to_list=True,
                                         check_chr_id_format=False)
    if not gen_bed_dic:
        return c_hits, c_uniq_hits, c_sites, region_size

    assert len(gen_seqs_dic) == len(gen_bed_dic), "len(gen_seqs_dic) != len(gen_bed_dic) for files %s and %s" %(gen_fa, gen_bed)
    c_sites = len(gen_bed_dic)

    region_size = get_uniq_gen_size(gen_bed)

    hit_reg_dic = {}
    hit_seqs_dic = {}
    gen_check_pos_dic = {}

    for site_id in gen_seqs_dic:
        seq = gen_seqs_dic[site_id]
        chr_id = gen_bed_dic[site_id][0]
        gen_s = gen_bed_dic[site_id][1]
        gen_e = gen_bed_dic[site_id][2]
        gen_pol = gen_bed_dic[site_id][5]
        site_id_hit_c = 0

        for match in re.finditer(motif, seq):

            hit = match.group()
            hit_s = match.start()
            hit_e = match.end()

            hit_gen_s = gen_s + hit_s
            hit_gen_e = gen_s + hit_e
            if gen_pol == "-":
                hit_gen_s = gen_e - hit_e
                hit_gen_e = gen_e - hit_s

            c_hits += 1

            check_id = "%s,%i-%i,%s" %(chr_id,hit_gen_s,hit_gen_e,gen_pol)
            if check_id in gen_check_pos_dic:
                continue
            else:
                site_id_hit_c += 1
                c_uniq_hits += 1
                hit_id = site_id + "," + str(site_id_hit_c)
                hit_reg_dic[hit_id] = [chr_id, hit_gen_s, hit_gen_e, gen_pol]
                hit_seqs_dic[hit_id] = hit
                gen_check_pos_dic[check_id] = 1

    if not hit_reg_dic:
        return c_hits, c_uniq_hits, c_sites, region_size

    if hits_out_bed:
        OUTHITBED = open(hits_out_bed, "w")
        for site_id in hit_reg_dic:
            OUTHITBED.write("%s\t%i\t%i\t%s\t0\t%s\n" %(hit_reg_dic[site_id][0], hit_reg_dic[site_id][1], hit_reg_dic[site_id][2], site_id, hit_reg_dic[site_id][3]))
        OUTHITBED.close()
    if hits_out_fa:
        OUTHITFA = open(hits_out_fa, "w")
        for site_id in hit_reg_dic:
            hit_seq = hit_seqs_dic[site_id]
            OUTHITFA.write(">%s\n%s\n" %(site_id, hit_seq))
        OUTHITFA.close()

    return c_hits, c_uniq_hits, c_sites, region_size


################################################################################

def get_tr_motif_hits(motif, tr_fa, tr_bed, gen_exon_bed,
                      hits_tr_con_gen_reg_bed=False,
                      hits_tr_con_gen_reg_split_bed=False,
                      hits_out_bed=False,
                      hits_out_fa=False):
    """
    Get motif hits on transcript sequences / sites.

    Return number of unique hits and effective region size.

    """
    c_hits = 0
    c_uniq_hits = 0
    c_sites = 0
    region_size = 0
    tr_seqs_dic = read_fasta_into_dic(tr_fa,
                                      dna=False,
                                      empty_check=True,
                                      skip_n_seqs=False)
    if not tr_seqs_dic:
        return c_hits, c_uniq_hits, c_sites, region_size

    tr_bed_dic = bed_read_rows_into_dic(tr_bed, to_list=True,
                                        check_chr_id_format=False)
    if not tr_bed_dic:
        return c_hits, c_uniq_hits, c_sites, region_size

    assert len(tr_seqs_dic) == len(tr_bed_dic), "len(tr_seqs_dic) != len(tr_bed_dic) for files %s and %s" %(tr_fa, tr_bed)
    c_sites = len(tr_bed_dic)

    region_size = get_uniq_tr_size(tr_bed, gen_exon_bed)

    hit_reg_dic = {}
    hit_seqs_dic = {}
    tr_check_pos_dic = {}

    for site_id in tr_seqs_dic:
        seq = tr_seqs_dic[site_id]
        tr_id = tr_bed_dic[site_id][0]
        tr_s = tr_bed_dic[site_id][1]
        tr_e = tr_bed_dic[site_id][2]
        site_id_hit_c = 0

        for match in re.finditer(motif, seq):
            hit = match.group()
            hit_s = match.start() # 0-based.
            hit_e = match.end() # 1-based.
            hit_tr_s = tr_s + hit_s
            hit_tr_e = tr_s + hit_e

            check_id = "%s,%i-%i" %(tr_id,hit_tr_s,hit_tr_e)
            if check_id in tr_check_pos_dic:
                continue
            else:
                site_id_hit_c += 1
                hit_id = site_id + "," + str(site_id_hit_c)
                hit_reg_dic[hit_id] = [tr_id, hit_tr_s, hit_tr_e]
                hit_seqs_dic[hit_id] = hit
                tr_check_pos_dic[check_id] = 1

    if not hit_reg_dic:
        return c_hits, c_uniq_hits, c_sites, region_size

    c_hits = len(hit_reg_dic)

    if hits_out_bed:
        OUTHITBED = open(hits_out_bed, "w")
        for site_id in hit_reg_dic:
            OUTHITBED.write("%s\t%i\t%i\t%s\t0\t+\n" %(hit_reg_dic[site_id][0], hit_reg_dic[site_id][1], hit_reg_dic[site_id][2], site_id))
        OUTHITBED.close()
    if hits_out_fa:
        OUTHITFA = open(hits_out_fa, "w")
        for site_id in hit_reg_dic:
            hit_seq = hit_seqs_dic[site_id]
            OUTHITFA.write(">%s\n%s\n" %(site_id, hit_seq))
        OUTHITFA.close()

    # Deduplicate transcript hits.
    dedup_hit_reg_dic = get_uniq_tr_regions(hit_reg_dic, gen_exon_bed,
                                            hits_tr_con_gen_reg_bed=hits_tr_con_gen_reg_bed,
                                            hits_tr_con_gen_reg_split_bed=hits_tr_con_gen_reg_split_bed)

    c_uniq_hits = len(dedup_hit_reg_dic)

    return c_hits, c_uniq_hits, c_sites, region_size


################################################################################

def get_uniq_tr_regions(tr_reg_dic, gen_exon_bed,
                        hits_tr_con_gen_reg_bed=False,
                        hits_tr_con_gen_reg_split_bed=False):
    """
    Get unique transcript regions, by mapping transcript sites on genome
    and remove duplicate regions (i.e., regions with same start+end+chr+pol).
    Also considers/removes identical split regions.

    tr_reg_dic:
        Transcript regions with mapping:
        region_id -> [tr_id, tr_s, tr_e] # tr_s 0-based.

    gen_exon_bed:
        genomic exon regions containing the transcript sites.
        Exon ID with format: transcriptid_e[1...n] with n == number of
        exons for transcript.

    gen_split_hits_bed_out:
        Provide file path to output transcript split hits on genome.

    tr_reg_dic:
    t1	995	1005	s1	0	+
    t2	995	1005	s2	0	+
    t2	995	1010	s3	0	+
    t3	995	1005	s4	0	+

    tr_reg_dic = {
    's1': ['t1', 995, 1005],
    's2': ['t2', 995, 1005],
    's3': ['t2', 995, 1010],
    's4': ['t3', 995, 1005]}

    test_uniq_tr_reg.gen_ex.bed:
    chr1	1000	2000	t1_e1	0	+
    chr1	3000	4000	t1_e2	0	+
    chr1	1000	2000	t2_e1	0	+
    chr1	3000	4000	t2_e2	0	+
    chr1	5000	6000	t2_e3	0	+
    chr1	1000	2000	t3_e2	0	-
    chr1	3000	4000	t3_e1	0	-

    >>> tr_reg_dic = {'s1': ['t1', 995, 1005], 's2': ['t2', 995, 1005], 's3': ['t2', 995, 1010], 's4': ['t3', 995, 1005]}
    >>> gen_exon_bed = "test_data/test_uniq_tr_reg.gen_ex.bed"
    >>> get_uniq_tr_regions(tr_reg_dic, gen_exon_bed)
    {'s1': ['t1', 995, 1005], 's3': ['t2', 995, 1010], 's4': ['t3', 995, 1005]}

    """
    assert tr_reg_dic, "given tr_reg_dic empty"

    # Generate .tmp files.
    random_id = uuid.uuid1()
    tmp1_bed = str(random_id) + ".tmp1.bed"
    random_id = uuid.uuid1()
    tmp2_bed = str(random_id) + ".tmp2.bed"
    random_id = uuid.uuid1()
    tmp3_bed = str(random_id) + ".tmp3.bed"

    # Read in exon region stats.
    id2chr_dic = {}
    id2s_dic = {}
    id2e_dic = {}
    id2pol_dic = {}
    exid2trid_dic = {}
    tr_exc_dic = {}  # count exon numbers.
    with open(gen_exon_bed) as f:
        for line in f:
            row = line.strip()
            cols = line.strip().split("\t")
            chr_id = cols[0]
            site_s = int(cols[1])
            site_e = int(cols[2])
            site_id = cols[3]
            site_pol = cols[5]
            id2chr_dic[site_id] = chr_id
            id2s_dic[site_id] = site_s
            id2e_dic[site_id] = site_e
            id2pol_dic[site_id] = site_pol
            if re.search(".+_e\d", site_id):
                m = re.search("(.+)_e\d", site_id)
                tr_id = m.group(1)
                exid2trid_dic[site_id] = tr_id
                if tr_id not in tr_exc_dic:
                    tr_exc_dic[tr_id] = 1
                else:
                    tr_exc_dic[tr_id] += 1
            else:
                assert False, "site ID \"%s\" missing added _e exon number" %(site_id)
    f.close()
    assert exid2trid_dic, "exid2trid_dic empty (nothing read in from %s)" %(gen_exon_bed)

    # Output transcript sites.
    OUTTBED = open(tmp1_bed, "w")
    for reg_id in tr_reg_dic:
        OUTTBED.write("%s\t%i\t%i\t%s\t0\t+\n" %(tr_reg_dic[reg_id][0], tr_reg_dic[reg_id][1], tr_reg_dic[reg_id][2], reg_id))
    OUTTBED.close()

    # Output exon regions with transcript coordinates.
    OUTTBED = open(tmp2_bed, "w")
    for tr_id in tr_exc_dic:
        ex_c = tr_exc_dic[tr_id]
        new_s = 0
        for i in range(ex_c):
            i += 1
            ex_id = tr_id + "_e" + str(i)
            gen_s = id2s_dic[ex_id]
            gen_e = id2e_dic[ex_id]
            ex_len = gen_e - gen_s
            tr_s = new_s
            tr_e = new_s + ex_len
            OUTTBED.write("%s\t%i\t%i\t%s\t0\t+\n" % (tr_id,tr_s,tr_e,ex_id))
            new_s = tr_e
    OUTTBED.close()

    # Overlap transcript sites with transcript exon regions.
    params = "-wb"
    check_cmd = "intersectBed -a " + tmp1_bed + " -b " + tmp2_bed + " " + params
    output = subprocess.getoutput(check_cmd)

    # Read in transcript site overlaps with transcript exon regions.
    site2c_dic = {}
    # Dictionaries for later outputting unique + split hits separately.
    siteid2pol_dic = {}
    siteid2sc_dic = {}
    partid2chrse_dic = {}
    siteid2parts_dic = {}
    part2siteids_dic = {}

    for line in output.split('\n'):
        cols = line.strip().split("\t")
        tr_id = cols[0]
        part_s = int(cols[1])
        part_e = int(cols[2])
        site_id = cols[3]
        site_sc = cols[4]
        ex_s = int(cols[7])
        ex_e = int(cols[8])
        ex_id = cols[9]
        ex_pol = id2pol_dic[ex_id]
        siteid2pol_dic[site_id] = ex_pol
        siteid2sc_dic[site_id] = site_sc
        if site_id in site2c_dic:
            site2c_dic[site_id] += 1
        else:
            site2c_dic[site_id] = 1
        # Hit part number.
        hit_c = site2c_dic[site_id]
        # Calculate genomic hit coordinates.
        # Plus strand case.
        gen_s = id2s_dic[ex_id] + part_s - ex_s
        gen_e = id2s_dic[ex_id] + part_e - ex_s
        # Minus strand case.
        if ex_pol == "-":
            gen_s = id2e_dic[ex_id] - part_e + ex_s
            gen_e = id2e_dic[ex_id] - part_s + ex_s
        # Store site_id parts.
        chrsepol = "%s,%i,%i,%s" %(id2chr_dic[ex_id],gen_s,gen_e, ex_pol)
        if chrsepol in part2siteids_dic:
            part2siteids_dic[chrsepol].append(site_id)
        else:
            part2siteids_dic[chrsepol] = [site_id]
        if site_id in siteid2parts_dic:
            siteid2parts_dic[site_id].append(chrsepol)
        else:
            siteid2parts_dic[site_id] = [chrsepol]

    # Sort parts and convert to string.
    siteid2pstr_dic = {}
    for site_id in siteid2parts_dic:
        siteid2parts_dic[site_id].sort()
        siteid2pstr_dic[site_id] = ""
        for pstr in siteid2parts_dic[site_id]:
            siteid2pstr_dic[site_id] += pstr
    pstr2siteids_dic = {}
    for site_id in siteid2pstr_dic:
        pstr = siteid2pstr_dic[site_id]
        if pstr in pstr2siteids_dic:
            pstr2siteids_dic[pstr].append(site_id)
        else:
            pstr2siteids_dic[pstr] = [site_id]

    #print("siteid2pstr_dic:", siteid2pstr_dic)
    #print("pstr2siteids_dic:", pstr2siteids_dic)

    # Get deduplicated site IDs.
    return_ids_dic = {}
    for pstr in pstr2siteids_dic:
        site_ids = pstr2siteids_dic[pstr]
        return_ids_dic[site_ids[0]] = 1

    if hits_tr_con_gen_reg_bed:
        GOUTBED = open(hits_tr_con_gen_reg_bed, "w")
    if hits_tr_con_gen_reg_split_bed:
        SPLITGOUTBED = open(hits_tr_con_gen_reg_split_bed, "w")

    return_tr_reg_dic = {}
    for site_id in tr_reg_dic:
        if site_id in return_ids_dic:
            return_tr_reg_dic[site_id] = tr_reg_dic[site_id]
            if hits_tr_con_gen_reg_bed:
                for split_reg in siteid2parts_dic[site_id]:
                    reg_info = split_reg.strip().split(",")
                    GOUTBED.write("%s\t%s\t%s\t%s\t0\t%s\n" %(reg_info[0], reg_info[1], reg_info[2], site_id, reg_info[3]))
            if hits_tr_con_gen_reg_split_bed:
                if len(siteid2parts_dic[site_id]) > 1:
                    for split_reg in siteid2parts_dic[site_id]:
                        reg_info = split_reg.strip().split(",")
                        SPLITGOUTBED.write("%s\t%s\t%s\t%s\t0\t%s\n" %(reg_info[0], reg_info[1], reg_info[2], site_id, reg_info[3]))

    if hits_tr_con_gen_reg_bed:
        GOUTBED.close()
    if hits_tr_con_gen_reg_split_bed:
        SPLITGOUTBED.close()

    assert return_tr_reg_dic, "return_tr_reg_dic empty (no regions to return)"
    assert len(return_ids_dic) == len(return_tr_reg_dic), "len(return_ids_dic) != len(return_tr_reg_dic)"

    if os.path.exists(tmp1_bed):
        os.remove(tmp1_bed)
    if os.path.exists(tmp2_bed):
        os.remove(tmp2_bed)

    return return_tr_reg_dic


################################################################################

def get_uniq_tr_size(tr_sites_bed, gen_exon_bed):
    """
    Get unique transcript space size, which the transcript sites inside
    tr_sites_bed cover. For this map the sites to genome (to gen_exon_bed,
    with genomic exon regions inside), and merge overlapping regions.
    Then sum up and return the region size.

    tr_sites_bed:
        Transcript sites on transcripts.

    gen_exon_bed:
        genomic exon regions containing the transcript sites.
        Exon ID with format: transcriptid_e[1...n] with n == number of
        exons for transcript.

    test_uniq_tr_size.tr_sites.bed:
    t1	990	1020	s1	0	+
    t2	10	40	s2	0	+
    t3	990	1020	s3	0	+

    test_uniq_tr_size.gen_ex.bed:
    chr1	1000	2000	t1_e1	0	+
    chr1	3000	4000	t1_e2	0	+
    chr1	3000	4000	t2_e1	0	+
    chr1	5000	6000	t3_e2	0	-
    chr1	7000	8000	t3_e1	0	-

    >>> tr_sites_bed = "test_data/test_uniq_tr_size.tr_sites.bed"
    >>> gen_exon_bed = "test_data/test_uniq_tr_size.gen_ex.bed"
    >>> get_uniq_tr_size(tr_sites_bed, gen_exon_bed)
    80

    """

    # Generate .tmp files.
    random_id = uuid.uuid1()
    tmp1_bed = str(random_id) + ".tmp1.bed"
    tmp2_bed = str(random_id) + ".tmp2.bed"

    # Read in exon region stats.
    id2chr_dic = {}
    id2s_dic = {}
    id2e_dic = {}
    id2pol_dic = {}
    exid2trid_dic = {}
    tr_exc_dic = {}  # count exon numbers.
    with open(gen_exon_bed) as f:
        for line in f:
            row = line.strip()
            cols = line.strip().split("\t")
            chr_id = cols[0]
            site_s = int(cols[1])
            site_e = int(cols[2])
            site_id = cols[3]
            site_pol = cols[5]
            id2chr_dic[site_id] = chr_id
            id2s_dic[site_id] = site_s
            id2e_dic[site_id] = site_e
            id2pol_dic[site_id] = site_pol
            if re.search(".+_e\d", site_id):
                m = re.search("(.+)_e\d", site_id)
                tr_id = m.group(1)
                exid2trid_dic[site_id] = tr_id
                if tr_id not in tr_exc_dic:
                    tr_exc_dic[tr_id] = 1
                else:
                    tr_exc_dic[tr_id] += 1
            else:
                assert False, "site ID \"%s\" missing added _e exon number" %(site_id)
    f.close()
    assert exid2trid_dic, "exid2trid_dic empty (nothing read in from %s)" %(gen_exon_bed)

    # Output exon regions with transcript coordinates.
    OUTTBED = open(tmp1_bed, "w")
    for tr_id in tr_exc_dic:
        ex_c = tr_exc_dic[tr_id]
        new_s = 0
        for i in range(ex_c):
            i += 1
            ex_id = tr_id + "_e" + str(i)
            gen_s = id2s_dic[ex_id]
            gen_e = id2e_dic[ex_id]
            ex_len = gen_e - gen_s
            tr_s = new_s
            tr_e = new_s + ex_len
            OUTTBED.write("%s\t%i\t%i\t%s\t0\t+\n" % (tr_id,tr_s,tr_e,ex_id))
            new_s = tr_e
    OUTTBED.close()

    # Overlap transcript sites with transcript exon regions.
    params = "-wb"
    check_cmd = "intersectBed -a " + tr_sites_bed + " -b " + tmp1_bed + " " + params
    output = subprocess.getoutput(check_cmd)

    # Read in transcript site overlaps with transcript exon regions.
    site2c_dic = {}
    # Dictionaries for later outputting unique + split hits separately.
    siteid2pol_dic = {}
    siteid2sc_dic = {}
    partid2chrse_dic = {}

    for line in output.split('\n'):
        cols = line.strip().split("\t")
        tr_id = cols[0]
        part_s = int(cols[1])
        part_e = int(cols[2])
        site_id = cols[3]
        site_sc = cols[4]
        ex_s = int(cols[7])
        ex_e = int(cols[8])
        ex_id = cols[9]
        ex_pol = id2pol_dic[ex_id]
        siteid2pol_dic[site_id] = ex_pol
        siteid2sc_dic[site_id] = site_sc
        if site_id in site2c_dic:
            site2c_dic[site_id] += 1
        else:
            site2c_dic[site_id] = 1
        # Hit part number.
        hit_c = site2c_dic[site_id]
        # Calculate genomic hit coordinates.
        # Plus strand case.
        gen_s = id2s_dic[ex_id] + part_s - ex_s
        gen_e = id2s_dic[ex_id] + part_e - ex_s
        # Minus strand case.
        if ex_pol == "-":
            gen_s = id2e_dic[ex_id] - part_e + ex_s
            gen_e = id2e_dic[ex_id] - part_s + ex_s
        # part ID.
        part_id = site_id + "_p" + str(hit_c)
        # Store chrse for each part ID.
        chrse = "%s\t%i\t%i" %(id2chr_dic[ex_id],gen_s,gen_e)
        partid2chrse_dic[part_id] = "%s\t%i\t%i" %(id2chr_dic[ex_id],gen_s,gen_e)

    # Output all hits (full and split).
    ALLHITSBED = open(tmp2_bed, "w")
    for site_id in site2c_dic:
        hit_c = site2c_dic[site_id]
        site_pol = siteid2pol_dic[site_id]
        site_sc = siteid2sc_dic[site_id]
        # For unique hit use site ID, for split hits use part IDs.
        if hit_c == 1:
            # Unique hits.
            part_id = site_id + "_p1"
            ALLHITSBED.write("%s\t%s\t%s\t%s\n" %(partid2chrse_dic[part_id],site_id,site_sc,site_pol))
        else:
            # Split hits.
            for i in range(hit_c):
                i += 1
                part_id = site_id + "_p" + str(i)
                ALLHITSBED.write("%s\t%s\t%s\t%s\n" %(partid2chrse_dic[part_id],part_id,site_sc,site_pol))

    ALLHITSBED.close()

    reg_len_sum = get_uniq_gen_size(tmp2_bed)

    if os.path.exists(tmp1_bed):
        os.remove(tmp1_bed)
    if os.path.exists(tmp2_bed):
        os.remove(tmp2_bed)

    return reg_len_sum


################################################################################

def get_uniq_gen_size(gen_sites_bed):
    """
    Get unique genomic space size, which the genomic sites inside
    gen_sites_bed cover.

    >>> gen_sites_bed = "test_data/test_gen_size.bed"
    >>> get_uniq_gen_size(gen_sites_bed)
    2500

    """

    params_str = '-s -c 4 -o distinct -delim ";"'
    check_cmd = "sort -k1,1 -k2,2n " + gen_sites_bed + " | mergeBed -i stdin " + params_str
    output = subprocess.getoutput(check_cmd)

    reg_len_sum = 0
    for line in output.split('\n'):
        cols = line.strip().split("\t")
        reg_s = int(cols[1])
        reg_e = int(cols[2])
        reg_len = reg_e - reg_s
        reg_len_sum += reg_len

    assert reg_len_sum, "no merged regions obtained from \"%s\"" %(check_cmd)
    return reg_len_sum


################################################################################

def bed_sort_file(in_bed, out_bed,
                  custom_params_str=False):
    """
    Use command line sort to sort the in_bed .bed file. Output sorted .bed
    file to out_bed.

    """
    # Parameter string.
    params_str = '-k1,1 -k2,2n'
    if custom_params_str:
        params_str = custom_params_str
    check_cmd = "sort " + params_str + " " + in_bed + " > " + out_bed
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "sort is complaining:\n%s\n%s" %(check_cmd, output)


################################################################################

def gtf_extract_exon_bed(in_gtf, out_bed,
                         out_intron_bed=False,
                         add_exon_id=False,
                         chr_id_style=1,
                         correct_min_ex_order=False,
                         tr2exc_dic=False,
                         tr_ids_dic=False):
    """
    Given a .gtf file with exon features, extract exon regions and store in
    .bed file. Optionally, a dictionary of transcript IDs can be provided,
    meaning that only exon regions from the given transcripts will be extracted.
    If out_intron_bed is set, an intronic regions .bed file will also be
    extracted, based on the exonic regions .bed information.

    Output .bed will look like this (note column 4 ID format with transcript
    ID followed by _e+exon_number):
    chr1	1000	2000	ENST001_e1	0	+
    chr1	3000	4000	ENST001_e2	0	+
    chr1	8000	9000	ENST002_e1	0	-
    chr1	6000	7000	ENST002_e2	0	-
    ...

    NOTE that function has been tested with .gtf files from Ensembl. .gtf files
    from different sources sometimes have a slightly different format, which
    could lead to incompatibilities / errors. See test files for format that
    works.

    Some tested Ensembl GTF files:
    Homo_sapiens.GRCh38.97.gtf.gz
    Mus_musculus.GRCm38.81.gtf.gz

    correct_min_ex_order:
        If set reverse number of exons for minus strand. This is necessary
        for some GTF files, which for minus strand transcripts assign lowest
        number to most upstream exon (same as for plus strand).

    >>> in_gtf = "test_data/map_test_in.gtf"
    >>> exp_out_bed = "test_data/gtf_exon_out_exp.bed"
    >>> exp_out_intron_bed = "test_data/gtf_intron_out_exp.bed"
    >>> out_bed = "test_data/gtf_exon_out.bed"
    >>> out_intron_bed = "test_data/gtf_intron_out.bed"
    >>> gtf_extract_exon_bed(in_gtf, out_bed, out_intron_bed=out_intron_bed)
    >>> diff_two_files_identical(out_bed, exp_out_bed)
    True
    >>> diff_two_files_identical(out_intron_bed, exp_out_intron_bed)
    True

    """

    # For reverse ording we need to have exon numbers.
    if correct_min_ex_order:
        assert tr2exc_dic, "tr2exc_dic needed if correct_min_ex_order True"

    # Output genomic exon regions.
    OUTBED = open(out_bed, "w")

    # Read in exon features from GTF file.
    c_gtf_ex_feat = 0
    # Start end coordinates of exons.
    exon_e_dic = {}
    exon_s_dic = {}
    # Transcript stats.
    tr2pol_dic = {}
    tr2chr_dic = {}
    # dic for sanity checking exon number order.
    tr2exon_nr_dic = {}

    proc_tr_dic = {}

    # Open GTF either as .gz or as text file.
    if re.search(".+\.gz$", in_gtf):
        f = gzip.open(in_gtf, 'rt')
    else:
        f = open(in_gtf, "r")
    for line in f:
        # Skip header.
        if re.search("^#", line):
            continue
        cols = line.strip().split("\t")
        chr_id = cols[0]
        feature = cols[2]
        feat_s = int(cols[3])
        feat_e = int(cols[4])
        feat_pol = cols[6]
        infos = cols[8]
        if not feature == "exon":
            continue

        # Restrict to standard chromosomes.
        new_chr_id = check_convert_chr_id(chr_id,
                                    id_style=chr_id_style)
        if not new_chr_id:
            continue
        else:
            chr_id = new_chr_id

        # Make start coordinate 0-base (BED standard).
        feat_s = feat_s - 1

        # Extract transcript ID.
        m = re.search('transcript_id "(.+?)"', infos)
        assert m, "transcript_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        transcript_id = m.group(1)
        # Extract exon number.
        m = re.search('exon_number "(\d+?)"', infos)
        # Try Gencode encoding.
        if not m:
            m = re.search('exon_number (\d+?);', infos)
        assert m, "exon_number entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        exon_nr = int(m.group(1))

        # Check if transcript ID is in transcript dic.
        if tr_ids_dic:
            if transcript_id not in tr_ids_dic:
                if len(tr_ids_dic) == len(proc_tr_dic):
                    break
                else:
                    continue
        proc_tr_dic[transcript_id] = 1

        # Store transcript stats.
        tr2pol_dic[transcript_id] = feat_pol
        tr2chr_dic[transcript_id] = chr_id

        # Check whether exon numbers are incrementing for each transcript ID.
        if not transcript_id in tr2exon_nr_dic:
            tr2exon_nr_dic[transcript_id] = exon_nr
        else:
            assert tr2exon_nr_dic[transcript_id] < exon_nr, "transcript ID \"%s\" without increasing exon number order in GTF file \"%s\"" %(transcript_id, in_gtf)
            tr2exon_nr_dic[transcript_id] = exon_nr

        # Reverse ordering for minus strand.
        if correct_min_ex_order and feat_pol == "-":
            assert transcript_id in tr2exc_dic, "transcript ID %s not in tr2exc_dic" %(transcript_id)
            exon_nr = tr2exc_dic[transcript_id] - exon_nr + 1
            assert exon_nr >= 1, "exon number < 1 assigned (%i) for transcript ID %s (# exons: %i)" %(exon_nr, transcript_id, tr2exc_dic[transcript_id])

        # Construct exon ID.
        exon_id = transcript_id + "_e" + str(exon_nr)

        # Count exon entry.
        c_gtf_ex_feat += 1

        # Store infos.
        exon_s_dic[exon_id] = feat_s
        exon_e_dic[exon_id] = feat_e

        # Output genomic exon region.
        OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,feat_s,feat_e,exon_id,feat_pol))

    OUTBED.close()
    f.close()

    # Check for read-in features.
    assert c_gtf_ex_feat, "no exon features read in from \"%s\"" %(in_gtf)

    # Output intron .bed.
    if out_intron_bed:
        tr2intron_nr_dic = {}
        OUTBED = open(out_intron_bed, "w")
        for tr_id in tr2pol_dic:
            tr_pol = tr2pol_dic[tr_id]
            chr_id = tr2chr_dic[tr_id]
            tr_c = tr2exon_nr_dic[tr_id]
            intron_c = 0
            tr2intron_nr_dic[tr_id] = 0
            # 1-exon transcripts, no introns.
            if tr_c == 1:
                continue
            ex_list = []
            for i in range(tr_c):
                ex_nr = i + 1
                ex_id = tr_id + "_e" + str(ex_nr)
                ex_list.append(ex_id)
            for i in range(len(ex_list)):
                ex1i = i
                ex2i = i + 1
                # At last exon, no more introns to add.
                if ex2i == len(ex_list):
                    break
                ex1id = ex_list[ex1i]
                ex2id = ex_list[ex2i]
                ex1s = exon_s_dic[ex1id]
                ex2s = exon_s_dic[ex2id]
                ex1e = exon_e_dic[ex1id]
                ex2e = exon_e_dic[ex2id]
                # Plus case.
                intron_s = ex1e
                intron_e = ex2s
                if tr_pol == "-":
                    intron_s = ex2e
                    intron_e = ex1s
                intron_id = tr_id + "_i" + str(ex2i)
                intron_c += 1

                OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,intron_s,intron_e,intron_id,tr_pol))
            tr2intron_nr_dic[tr_id] = intron_c
        OUTBED.close()
        # Sanity check exon + intron numbers.
        for tr_id in tr2exon_nr_dic:
            exon_nr = tr2exon_nr_dic[tr_id]
            intron_nr = tr2intron_nr_dic[tr_id]
            assert (exon_nr-1) == intron_nr, "intron number != exon number - 1 for \"%s\" (%i != %i - 1)" %(tr_id, intron_nr, exon_nr)


################################################################################

def gtf_extract_exon_numbers(in_gtf,
                             tr_ids_dic=False,
                             chr_id_style=1):
    """
    Given a .gtf file with exon features, return dictionary with transcript
    ID and exon number.

    tr_ids_dic:
    Give tr_ids_dic dictionary with transcript IDs to keep.

    >>> in_gtf = "test_data/test_border_annot.gtf"
    >>> tr_ids_dic = {'ENST1': 1, 'ENST2': 1, 'ENST3': 1}
    >>> gtf_extract_exon_numbers(in_gtf, tr_ids_dic=tr_ids_dic)
    {'ENST1': 1, 'ENST2': 2, 'ENST3': 2}

    """

    # Transcript ID to exon count dic.
    tr2exc_dic = {}
    # dic for sanity checking exon number order.
    tr2exon_nr_dic = {}

    # Open GTF either as .gz or as text file.
    if re.search(".+\.gz$", in_gtf):
        f = gzip.open(in_gtf, 'rt')
    else:
        f = open(in_gtf, "r")
    for line in f:
        # Skip header.
        if re.search("^#", line):
            continue
        cols = line.strip().split("\t")
        chr_id = cols[0]
        feature = cols[2]
        infos = cols[8]
        if not feature == "exon":
            continue

        # Restrict to standard chromosomes.
        new_chr_id = check_convert_chr_id(chr_id,
                                    id_style=chr_id_style)
        if not new_chr_id:
            continue
        else:
            chr_id = new_chr_id

        # Extract transcript ID.
        m = re.search('transcript_id "(.+?)"', infos)
        assert m, "transcript_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        transcript_id = m.group(1)
        # Extract exon number.
        m = re.search('exon_number "(\d+?)"', infos)
        # Try Gencode encoding.
        if not m:
            m = re.search('exon_number (\d+?);', infos)
        assert m, "exon_number entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        exon_nr = int(m.group(1))

        if tr_ids_dic:
            if transcript_id not in tr_ids_dic:
                continue

        # Count exon numbers.
        if not transcript_id in tr2exc_dic:
            tr2exc_dic[transcript_id] = 1
        else:
            tr2exc_dic[transcript_id] += 1

        # Check whether exon numbers are incrementing for each transcript ID.
        if not transcript_id in tr2exon_nr_dic:
            tr2exon_nr_dic[transcript_id] = exon_nr
        else:
            assert tr2exon_nr_dic[transcript_id] < exon_nr, "transcript ID \"%s\" without increasing exon number order in GTF file \"%s\"" %(transcript_id, in_gtf)
            tr2exon_nr_dic[transcript_id] = exon_nr
    f.close()

    # Check for read-in content.
    assert tr2exc_dic, "no exon features read in from \"%s\"" %(in_gtf)
    # Return to the castle.
    return tr2exc_dic


################################################################################

def gtf_check_exon_order(in_gtf):
    """
    Check exon_number ordering. Return True if ordering for minus strand
    and plus strand is different (i.e. for minus exon_number 1 is most downstream).
    Return False, if upstream to downstream order numbering is used for both
    plus and minus strand transcripts.

    >>> test_true_gtf = "test_data/test_order_true.gtf"
    >>> test_false_gtf = "test_data/test_order_false.gtf"
    >>> gtf_check_exon_order(test_true_gtf)
    True
    >>> gtf_check_exon_order(test_false_gtf)
    False

    """
    tr2exon_nr_dic = {}
    tr2exon_s_dic = {}

    check = 6666

    if re.search(".+\.gz$", in_gtf):
        f = gzip.open(in_gtf, 'rt')
    else:
        f = open(in_gtf, "r")
    for line in f:
        # Skip header.
        if re.search("^#", line):
            continue
        cols = line.strip().split("\t")
        chr_id = cols[0]
        feature = cols[2]
        feat_s = int(cols[3])
        feat_e = int(cols[4])
        feat_pol = cols[6]
        infos = cols[8]

        if feature != "exon":
            continue

        # Transcript ID.
        m = re.search('transcript_id "(.+?)"', infos)
        assert m, "transcript_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        transcript_id = m.group(1)
        # Exon number.
        m = re.search('exon_number "(\d+?)"', infos)
        # Try Gencode encoding.
        if not m:
            m = re.search('exon_number (\d+?);', infos)
        assert m, "exon_number entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        exon_nr = int(m.group(1))

        # Check whether exon numbers are incrementing for each transcript ID.
        if transcript_id not in tr2exon_nr_dic:
            tr2exon_nr_dic[transcript_id] = exon_nr
        else:
            assert tr2exon_nr_dic[transcript_id] < exon_nr, "transcript ID \"%s\" without increasing exon number order in GTF file \"%s\"" %(transcript_id, in_gtf)
            tr2exon_nr_dic[transcript_id] = exon_nr

        # Check ordering of exons for minus strand transcripts.
        if transcript_id not in tr2exon_s_dic:
            tr2exon_s_dic[transcript_id] = feat_s
        else:
            if feat_pol == "-":
                if tr2exon_s_dic[transcript_id] > feat_s:
                    check = True
                else:
                    check = False
                break
            elif feat_pol == "+":
                assert tr2exon_s_dic[transcript_id] < feat_s, "transcript ID \"%s\" on plus strand but exon region coordinates are decreasing" %(transcript_id)
    f.close()

    assert check != 6666, "no minus strand exon regions found in GTF file %s" %(in_gtf)
    return check


################################################################################

def gtf_extract_unique_exon_bed(in_gtf, out_bed,
                                no_tbt_filter=False,
                                tr_ids_dic=False,
                                biotype_filter=True,
                                correct_min_ex_order=False,
                                tr2exc_dic=False,
                                reject_tr_bt_dic=False,
                                gene_feat_check=False,
                                chr_id_style=1,
                                next2exids_dic=None,
                                exid2trid_dic=None,
                                trid2exc_dic=None,
                                trid2tsl_dic=None,
                                trid2tbt_dic=None,
                                trid2gid_dic=None,
                                trid2gna_dic=None,
                                trid2gbt_dic=None,
                                trid2len_dic=None,
                                next2reg_dic=None):
    """
    Given a .gtf file with exon features, extract exon unique (!) regions.
    Since the Ensembl exon_id regions are not unique regarding their genomic
    coordinates, create own IDs each representing one unique genomic region
    (unique start+end+strand info).

    Output .bed will look like this (column 4 ID == new exon ID):
    chr1	1000	2000	NEXT1	0	+
    chr1	3000	4000	NEXT2	0	+
    chr1	8000	9000	NEXT3	0	-
    chr1	6000	7000	NEXT4	0	-
    ...

    Note there are differences in GTF format between Gencode, Ensembl,
    and NCBI ...

    GTF content looking for:
    gene_id "id";
    transcript_id "id";
    transcript_biotype "id";
    ...

    correct_min_ex_order:
        If set reverse number of exons for minus strand. This is necessary
        for some GTF files, which for minus strand transcripts assign lowest
        number to most upstream exon (same as for plus strand).

    Ensembl GTF example:
    1       havana  gene    11869   14409   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene";
    1       havana  transcript      11869   14409   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; tag "basic"; transcript_support_level "1";
    1       havana  exon    11869   12227   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "1"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; exon_id "ENSE00002234944"; exon_version "1"; tag "basic"; transcript_support_level "1";
    1       havana  exon    12613   12721   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "2"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; exon_id "ENSE00003582793"; exon_version "1"; tag "basic"; transcript_support_level "1";
    1       havana  exon    13221   14409   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "3"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; exon_id "ENSE00002312635"; exon_version "1"; tag "basic"; transcript_support_level "1";

    StringTie custom GTF example:
    chr1	StringTie	transcript	149054033	149082430	.	-.	transcript_id "ENST00000613595.4"; gene_id "MSTRG.1496"; gene_name "NBPF9"; xloc "XLOC_004117"; ref_gene_id "ENSG00000269713.7"; cmp_ref "ENST00000613595.4"; class_code "="; tss_id "TSS10003";
    chr1	StringTie	exon	149054033	149055899	.	-	.transcript_id "ENST00000613595.4"; gene_id "MSTRG.1496"; exon_number "1";
    chr1	StringTie	exon	149056512	149056620	.	-	.transcript_id "ENST00000613595.4"; gene_id "MSTRG.1496"; exon_number "2";
    chr1	StringTie	exon	149060523	149060695	.	-	.transcript_id "ENST00000613595.4"; gene_id "MSTRG.1496"; exon_number "3";
    chr1	StringTie	exon	149061332	149061383	.	-	.transcript_id "ENST00000613595.4"; gene_id "MSTRG.1496"; exon_number "4";
    chr19	StringTie	transcript	3359583	3469217	.	+	.	transcript_id "MSTRG.11612.3"; gene_id "MSTRG.11612"; gene_name "NFIC"; xloc "XLOC_025438"; cmp_ref "ENST00000641145.1"; class_code "j"; tss_id "TSS63599";
    chr19	StringTie	exon	3359583	3359685	.	+	.	transcript_id "MSTRG.11612.3"; gene_id "MSTRG.11612"; exon_number "1";
    chr19	StringTie	exon	3381712	3382243	.	+	.	transcript_id "MSTRG.11612.3"; gene_id "MSTRG.11612"; exon_number "2";
    chr19	StringTie	exon	3425106	3425177	.	+	.	transcript_id "MSTRG.11612.3"; gene_id "MSTRG.11612"; exon_number "3";
    chr19	StringTie	exon	3433518	3433592	.	+	.	transcript_id "MSTRG.11612.3"; gene_id "MSTRG.11612"; exon_number "4";

        reject_tr_bt_dic = {
            "nonsense_mediated_decay" : 1,
            "retained_intron" : 1,
            "non_stop_decay" : 1,
            "processed_transcript" : 1
        }

    """
    # For reverse ording we need to have exon numbers.
    if correct_min_ex_order:
        assert tr2exc_dic, "tr2exc_dic needed if correct_min_ex_order True"

    # dic for sanity checking exon number order.
    tr2exon_nr_dic = {}

    # Reject transcripts with these biotypes.
    if not reject_tr_bt_dic:
        reject_tr_bt_dic = {
            "nonsense_mediated_decay" : 1,
            "retained_intron" : 1,
            "non_stop_decay" : 1,
            "processed_transcript" : 1
        }
    if no_tbt_filter:
        reject_tr_bt_dic = {}
    # Remember rejected transcript IDs.
    reject_tr_ids_dic = {}

    # Seen transcript IDs.
    proc_tr_dic = {}

    ok_features_dic = {'transcript': 1, 'exon': 1}

    # Store exon ID region data.
    reg_str_dic = {}
    reg_str_exon_ids_dic = {}

    # Open GTF either as .gz or as text file.
    if re.search(".+\.gz$", in_gtf):
        f = gzip.open(in_gtf, 'rt')
    else:
        f = open(in_gtf, "r")
    for line in f:
        # Skip header.
        if re.search("^#", line):
            continue
        cols = line.strip().split("\t")
        chr_id = cols[0]
        feature = cols[2]
        feat_s = int(cols[3])
        feat_e = int(cols[4])
        feat_pol = cols[6]
        infos = cols[8]

        if feature not in ok_features_dic:
            continue

        # Gene ID.
        m = re.search('gene_id "(.+?)"', infos)
        assert m, "gene_id entry missing for \"exon\" feature in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        gene_id = m.group(1)

        # Transcript ID.
        m = re.search('transcript_id "(.+?)"', infos)
        assert m, "transcript_id entry missing for \"exon\" feature in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        transcript_id = m.group(1)

        if feature == "transcript":

            if tr_ids_dic:
                if transcript_id not in tr_ids_dic:
                    continue

            # Only for GTF files with gene feature present (GENCODE, Ensembl .. ).
            gene_name = "-"
            gene_biotype = "-"
            tr_biotype = "-"
            m = re.search('gene_name "(.+?)"', infos)
            if m:
                gene_name = m.group(1)
            else:
                if gene_feat_check:
                    assert False, "gene_name entry missing for \"exon\" feature in GTF file \"%s\", line \"%s\"" %(in_gtf, line)

            m = re.search('gene_biotype "(.+?)"', infos)
            # Try Gencode encoding.
            if not m:
                m = re.search('gene_type "(.+?)"', infos)
            if m:
                gene_biotype = m.group(1)
            else:
                if gene_feat_check:
                    assert False, "gene_biotype or gene_type entry missing for \"exon\" feature in GTF file \"%s\", line \"%s\"" %(in_gtf, line)

            m = re.search('transcript_biotype "(.+?)"', infos)
            # Try Gencode encoding.
            if not m:
                m = re.search('transcript_type "(.+?)"', infos)
            if m:
                tr_biotype = m.group(1)
            else:
                if gene_feat_check:
                    assert False, "transcript_biotype or transcript_type entry missing for \"exon\" feature in GTF file \"%s\", line \"%s\"" %(in_gtf, line)

            # Transcript biotype filter.
            if biotype_filter:
                if tr_biotype in reject_tr_bt_dic:
                    reject_tr_ids_dic[transcript_id] = 1
                    continue

            # Get transcript support level (TSL).
            m = re.search('transcript_support_level "(.+?)"', infos)
            tr_tsl = "NA"
            if m:
                tr_tsl = m.group(1)
            # Gencode basic tag there?
            gc_basic = False
            if re.search('tag "basic"', infos):
                gc_basic = True
            ccds = False
            if re.search('tag "CCDS"', infos):
                ccds = True

            if trid2gid_dic is not None:
                trid2gid_dic[transcript_id] = gene_id
            if trid2gna_dic is not None:
                trid2gna_dic[transcript_id] = gene_name
            if trid2gbt_dic is not None:
                trid2gbt_dic[transcript_id] = gene_biotype
            if trid2tsl_dic is not None:
                trid2tsl_dic[transcript_id] = [tr_tsl, gc_basic, ccds]
            if trid2tbt_dic is not None:
                trid2tbt_dic[transcript_id] = tr_biotype


        elif feature == "exon":

            # If specified, only extract exons from these transcripts.
            if tr_ids_dic:
                if transcript_id not in tr_ids_dic:
                    if len(tr_ids_dic) == len(proc_tr_dic):
                        break
                    else:
                        continue

            proc_tr_dic[transcript_id] = 1

            # Transcript biotype filter.
            if transcript_id in reject_tr_ids_dic:
                continue

            # Restrict to standard chromosomes.
            new_chr_id = check_convert_chr_id(chr_id,
                                    id_style=chr_id_style)
            if not new_chr_id:
                continue
            else:
                chr_id = new_chr_id

            # Make start coordinate 0-base (BED standard).
            feat_s = feat_s - 1
            # Feature length.
            feat_l = feat_e - feat_s

            # Transcript length.
            if trid2len_dic is not None:
                if transcript_id in trid2len_dic:
                    trid2len_dic[transcript_id] += feat_l
                else:
                    trid2len_dic[transcript_id] = feat_l

            # Extract exon number.
            m = re.search('exon_number "(\d+?)"', infos)
            # Try Gencode encoding.
            if not m:
                m = re.search('exon_number (\d+?);', infos)
            assert m, "exon_number entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
            exon_nr = int(m.group(1))

            # Check whether exon numbers are incrementing for each transcript ID.
            if not transcript_id in tr2exon_nr_dic:
                tr2exon_nr_dic[transcript_id] = exon_nr
            else:
                assert tr2exon_nr_dic[transcript_id] < exon_nr, "transcript ID \"%s\" without increasing exon number order in GTF file \"%s\"" %(transcript_id, in_gtf)
                tr2exon_nr_dic[transcript_id] = exon_nr

            # Reverse ordering for minus strand.
            if correct_min_ex_order and feat_pol == "-":
                assert transcript_id in tr2exc_dic, "transcript ID %s not in tr2exc_dic" %(transcript_id)
                exon_nr = tr2exc_dic[transcript_id] - exon_nr + 1
                assert exon_nr >= 1, "exon number < 1 assigned (%i) for transcript ID %s (# exons: %i)" %(exon_nr, transcript_id, tr2exc_dic[transcript_id])
            # Construct exon ID.
            exon_id = transcript_id + "_e" + str(exon_nr)

            # Store exon data.
            check_reg_str = "%s,%i,%i,%s" %(chr_id,feat_s,feat_e,feat_pol)
            reg_str_dic[check_reg_str] = 1
            # Store exon IDs for this particular region.
            if check_reg_str in reg_str_exon_ids_dic:
                reg_str_exon_ids_dic[check_reg_str].append(exon_id)
            else:
                reg_str_exon_ids_dic[check_reg_str] = [exon_id]
            # Exon ID to transcript ID mapping.
            if exid2trid_dic is not None:
                exid2trid_dic[exon_id] = transcript_id

    f.close()

    assert reg_str_dic, "no exon regions read in"

    if reject_tr_ids_dic:
        print("# rejected transcript IDs (biotype filter):  %i" %(len(reject_tr_ids_dic)))

    # Store transcript exon numbers.
    if trid2exc_dic is not None:
        for tr_id in tr2exon_nr_dic:
            trid2exc_dic[tr_id] = tr2exon_nr_dic[tr_id]

    # Output genomic exon regions.
    OUTBED = open(out_bed, "w")

    c_ex = 0
    for reg_str in reg_str_dic:
        cols = reg_str.split(",")
        c_ex += 1
        ex_id = "NEXT" + str(c_ex)
        if next2exids_dic is not None:
            next2exids_dic[ex_id] = reg_str_exon_ids_dic[reg_str]
        if next2reg_dic is not None:
            next2reg_dic[ex_id] = [cols[0], int(cols[1]), int(cols[2]), cols[3]]
        OUTBED.write("%s\t%s\t%s\t%s\t0\t%s\n" % (cols[0], cols[1], cols[2], ex_id, cols[3]))
    OUTBED.close()


################################################################################

def get_exons_fully_ol_with_isr_introns(isr_intron_reg_dic, next2reg_dic,
                                        next_ids_dic=None,
                                        reg2cov_dic=None,
                                        max_read2isr_ratio=8,
                                        next2top_isrn_dic=False,
                                        tmp_out_folder=False,
                                        min_isrc=5):
    """
    Get exons fully overlapping with ISR-containing introns
    (i.e., the span of these introns).

    isr_intron_reg_dic:
        Intronic region string "chrid,s,e,pol" -> ISR count
    next2reg_dic:
        NEXT exon ID -> genomic region
        next2reg_dic = [chrid, s, e, pol]
    next_ids_dic:
        Dictionary of exon (NEXT) IDs to include in overlap calculations.
    min_isrc:
        Minimum ISR count an intronic region needs to be included in
        overlap calculation.
    reg2cov_dic:
        ["chr,s,e,pol"] -> read overlap count / coverage.
        Used to get read coverage of intron regions.
    max_read2isr_ratio:
        Maximum ratio of # total reads / # ISR reads, for filtering out
        fully overlapping exons. ISR introns with smaller ratios are not
        considered for filtering.
    next2top_isrn_dic:
        NEXT ID to top ISRN count mapping. If given, a fully overlapping
        exon with # ISRN reads >= # ISR reads of the overlapping intron
        does not get returned.
    tmp_out_folder:
        Provide output folder to store tmp files in.

    """
    assert isr_intron_reg_dic, "isr_intron_reg_dic empty"
    assert next2reg_dic, "next2reg_dic empty"

    random_id = uuid.uuid1()
    next_tmp_bed = str(random_id) + ".next_regions.tmp.bed"
    random_id = uuid.uuid1()
    isr_intron_tmp_bed = str(random_id) + ".intron_regions.tmp.bed"
    if tmp_out_folder:
        next_tmp_bed = tmp_out_folder + "/" + next_tmp_bed
        isr_intron_tmp_bed = tmp_out_folder + "/" + isr_intron_tmp_bed

    bed_write_reg_list_to_file(next2reg_dic, next_tmp_bed,
                               id2out_dic=next_ids_dic)

    ISROUTBED = open(isr_intron_tmp_bed, "w")
    reg_out_c = 0

    for intron_reg in isr_intron_reg_dic:
        cols = intron_reg.split(",")
        chr_id = cols[0]
        reg_s = cols[1]
        reg_e = cols[2]
        reg_pol = cols[3]
        # number of ISR reads in the intron_reg.
        reg_sc = isr_intron_reg_dic[intron_reg]
        if reg_sc >= min_isrc:
            if reg2cov_dic is not None and intron_reg in reg2cov_dic:
                # assert intron_reg in reg2cov_dic, "intron region \"%s\" not in reg2cov_dic" %(intron_reg)
                c_intron_reads = reg2cov_dic[intron_reg]
                read2isr_ratio = c_intron_reads / reg_sc

                if read2isr_ratio >= max_read2isr_ratio:
                    continue

            reg_out_c += 1
            reg_id = "intron_%i" %(reg_out_c)

            ISROUTBED.write("%s\t%s\t%s\t%s\t%i\t%s\n" %(chr_id, reg_s, reg_e, reg_id, reg_sc, reg_pol))

    if not reg_out_c:
        return {}
    # assert reg_out_c, "no introns output for ISR count threshold %i" %(min_isrc)

    ISROUTBED.close()

    # Overlap calculation to get exons fully containing introns.
    params = " -s -F 1.0 -wb"
    check_cmd = "intersectBed -a " + next_tmp_bed + " -b " + isr_intron_tmp_bed + " " + params
    output = subprocess.getoutput(check_cmd)

    """
    -F 1.0
    Only full overlaps (fraction of B).

    Example:

    $ intersectBed -a exons.bed -b introns.bed -s -F 1.0 -wb
    chr1	3400	3600	e2	0	+	chr1	3400	3600	i1	5	+

    """
    rem_nexts_dic = {}

    # If there are full overlaps.
    if output:
        for line in output.split('\n'):
            cols = line.strip().split("\t")
            next_id = cols[3]
            intron_id = cols[9]
            isrc = int(cols[10])
            if next2top_isrn_dic:
                assert next_id in next2top_isrn_dic, "NEXT ID %s not in next2top_isrn_dic" %(next_id)
                next_isrn = next2top_isrn_dic[next_id]

                if next_isrn >= isrc:
                    continue

            rem_nexts_dic[next_id] = 1

    if os.path.exists(isr_intron_tmp_bed):
        os.remove(isr_intron_tmp_bed)
    if os.path.exists(next_tmp_bed):
        os.remove(next_tmp_bed)

    return rem_nexts_dic


################################################################################

def output_ref_lengths(ref_len_dic, ref_len_out):
    """
    Output reference lengths.

    Output file format:
    chr1	210
    chr2	200
    ...

    """
    assert ref_len_dic, "ref_len_dic empty"

    OUTREFLEN = open(ref_len_out, "w")

    for ref_id in ref_len_dic:
        ref_len = ref_len_dic[ref_id]
        OUTREFLEN.write("%s\t%i\n" %(ref_id, ref_len))
    OUTREFLEN.close()


################################################################################

def read_in_id_val_cols(in_file,
                        id_col=3,
                        val_col=4,
                        val_type=1):
    """
    Read in ID value combination into dic. id_col + val_col 0-based.

    val_type:
        1 : float
        2 : integer

    """
    id2val_dic = {}

    with open(in_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            id_str = cols[id_col]
            val = cols[val_col]
            if val_type == 1:
                val = float(val)
            elif val_type == 2:
                val = int(val)
            assert id_str not in id2val_dic, "ID %s (col: %i) appears > 1 in file %s" %(id_str, id_col, in_file)
            id2val_dic[id_str] = val
    f.closed
    assert id2val_dic, "id2val_dic empty (nothing read in)"
    return id2val_dic


################################################################################

def read_in_ref_lengths(ref_len_in,
                        ref_len_dic=False):
    """
    Read in reference lengths from file.

    Input file format:
    chr1	210
    chr2	200
    ...

    """
    if not ref_len_dic:
        ref_len_dic = {}

    with open(ref_len_in) as f:
        for line in f:
            cols = line.strip().split("\t")
            ref_id = cols[0]
            ref_len = int(cols[1])
            if ref_id in ref_len_dic:
                assert ref_len == ref_len_dic[ref_id], "reference ID %s with differing lengths read in from %s (new != existing length, %i != %i)" %(ref_id, ref_len_in, ref_len, ref_len_dic[ref_id])
            else:
                ref_len_dic[ref_id] = ref_len
    f.closed
    assert ref_len_dic, "ref_len_dic empty (nothing read in)"
    return ref_len_dic


################################################################################

def read_in_sel_tr_regs(exon_sites_sel_tr_bed, id2selreg_dic,
                        id2dataset_dic=None,
                        dataset_id=1):
    """
    Read in selected transcript regions for each exonic site ID.
    Store in id2selreg_dic, format:
    site_id -> [tr_id, s, e, score]

    """

    with open(exon_sites_sel_tr_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            tr_id = cols[0]
            tr_s = int(cols[1])
            tr_e = int(cols[2])
            sitetrsc_id = cols[3]
            site_sc = cols[4]
            sitetrsc_cols = sitetrsc_id.split(",")
            site_id = sitetrsc_cols[0]
            tr_id_check = sitetrsc_cols[1]
            assert tr_id == tr_id_check, "tr_id != tr_id_check for site ID %s (%s != %s, input file: %s)" %(sitetrsc_id, tr_id, tr_id_check, exon_sites_all_tr_bed)
            if id2dataset_dic is not None:
                id2dataset_dic[site_id] = dataset_id
            assert site_id not in id2selreg_dic, "site ID %s already read in. Site IDs need to be unique for merging (also in between --in datasets!)" %(site_id)
            id2selreg_dic[site_id] = [tr_id, tr_s, tr_e, site_sc]
    f.closed
    assert id2selreg_dic, "id2selreg_dic empty (nothing read in)"


################################################################################

def get_site_to_pair_id_mapping(all_sites_igv_tsv, id2pairid_dic,
                                id2regtype_dic=None):
    """
    Get site ID to pair ID mapping for site IDs that form exon border pair.

    all_sites_igv_tsv content:
    site_id	assigned_region_type	new_merged_id	exon_gc_filter	igv_region	strand
    PUM2_K562_IDR_0001	exon_tc	-	-	chr1:42,176,874-42,176,938	-
    PUM2_K562_IDR_0002	exon_tc	-	-	chr5:72,914,320-72,914,373	+
    ...

    id2regtype_dic:
        site ID -> region type mapping.
        Region types: exon_tc, exon_gc, intron, intergen

    """

    with open(all_sites_igv_tsv) as f:
        for line in f:
            row = line.strip()
            cols = line.strip().split("\t")
            if cols[0] == "site_id":
                continue
            site_id = cols[0]
            reg_type = cols[1]
            pair_id = cols[2]
            if pair_id != "-":
                id2pairid_dic[site_id] = pair_id
            if id2regtype_dic is not None:
                id2regtype_dic[site_id] = reg_type
    f.close()


################################################################################

def read_in_genomic_regs(exon_sites_gen_regs_bed,
                         check_chr_ids_dic=False,
                         id2genreg_dic=False):
    """
    Read in genomic regions for each exonic site ID.
    Store in id2genreg_dic, format:
    site_id -> [chr_id, s, e, pol]

    >>> in_bed = "test_data/test3.bed"
    >>> id2genreg_dic = {}
    >>> read_in_genomic_regs(in_bed, id2genreg_dic=id2genreg_dic)
    {'CLIP1': ['chr1', 10, 20, '+'], 'CLIP2': ['chr1', 30, 45, '-']}

    """
    if not id2genreg_dic:
        id2genreg_dic = {}

    with open(exon_sites_gen_regs_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            chr_id = cols[0]
            gen_s = int(cols[1])
            gen_e = int(cols[2])
            site_id = cols[3]
            gen_pol = cols[5]
            assert site_id not in id2genreg_dic, "site ID %s already read in. Site IDs need to be unique for merging (also in between --in datasets!)" %(site_id)
            # Sanity checking if input dataset chromosome IDs found inside --gtf.
            if check_chr_ids_dic:
                assert chr_id in check_chr_ids_dic, "chromosome ID \"%s\" from \"%s\" not found in --gtf file" %(chr_id, exon_sites_gen_regs_bed)
            id2genreg_dic[site_id] = [chr_id, gen_s, gen_e, gen_pol]
    f.closed
    assert id2genreg_dic, "id2genreg_dic empty (nothing read in)"
    return id2genreg_dic


################################################################################

def read_in_all_tr_regs(exon_sites_all_tr_bed, id2allreg_dic,
                        sitetrid2sc_dic, id2bed_sc_dic):

    """
    Read in all transcript regions for each exonic site ID.
    Store in id2selreg_dic, format:
    site_id -> ["tr_id,s,e", "tr_id,s,e", ... ]

    test_all_tr_regs.bed:
    T1	100	110	s1,T1,5	0.1	+
    T2	100	110	s1,T2,6	0.1	+
    T1	200	210	s2,T1,7	0.2	+

    >>> in_bed = "test_data/test_all_tr_regs.bed"
    >>> id2allreg_dic = {}
    >>> sitetrid2sc_dic = {}
    >>> id2bed_sc_dic = {}
    >>> read_in_all_tr_regs(in_bed, id2allreg_dic, sitetrid2sc_dic, id2bed_sc_dic)
    >>> id2allreg_dic
    {'s1': ['T1,100,110', 'T2,100,110'], 's2': ['T1,200,210']}
    >>> sitetrid2sc_dic
    {'s1,T1': 5, 's1,T2': 6, 's2,T1': 7}
    >>> id2bed_sc_dic
    {'s1': '0.1', 's2': '0.2'}

    """

    with open(exon_sites_all_tr_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            tr_id = cols[0]
            tr_s = cols[1]
            tr_e = cols[2]
            sitetrsc_id = cols[3]
            site_sc = cols[4]
            sitetrsc_cols = sitetrsc_id.split(",")
            site_id = sitetrsc_cols[0]
            tr_id_check = sitetrsc_cols[1]
            comb_sc = int(sitetrsc_cols[2])
            assert tr_id == tr_id_check, "tr_id != tr_id_check for site ID %s (%s != %s, input file: %s)" %(sitetrsc_id, tr_id, tr_id_check, exon_sites_all_tr_bed)
            sitetrid = "%s,%s" %(site_id, tr_id)
            sitetrid2sc_dic[sitetrid] = comb_sc
            app_str = "%s,%s,%s" %(tr_id, tr_s, tr_e)
            if site_id in id2allreg_dic:
                id2allreg_dic[site_id].append(app_str)
            else:
                id2allreg_dic[site_id] = [app_str]
            id2bed_sc_dic[site_id] = site_sc

    f.closed
    assert id2allreg_dic, "id2allreg_dic empty (nothing read in)"


################################################################################

def check_convert_chr_id(chr_id,
                         id_style=1):
    """
    Check and convert chromosome IDs to format:
    chr1, chr2, chrX, ...
    If chromosome IDs like 1,2,X, .. given, convert to chr1, chr2, chrX ..
    Return False if given chr_id not standard and not convertable.

    Filter out scaffold IDs like:
    GL000009.2, KI270442.1, chr14_GL000009v2_random
    chrUn_KI270442v1 ...

    id_style:
        Defines to which style chromosome IDs should be converted to.
        0: Do not convert at all, just return chr_id.
        1: to chr1,chr2,...,chrM style.
        2: to 1,2,...,MT style.

    >>> chr_id = "chrX"
    >>> check_convert_chr_id(chr_id)
    'chrX'
    >>> chr_id = "4"
    >>> check_convert_chr_id(chr_id)
    'chr4'
    >>> chr_id = "MT"
    >>> check_convert_chr_id(chr_id)
    'chrM'
    >>> chr_id = "GL000009.2"
    >>> check_convert_chr_id(chr_id)
    False
    >>> chr_id = "chrUn_KI270442v1"
    >>> check_convert_chr_id(chr_id)
    False
    >>> chr_id = "chr2R"
    >>> check_convert_chr_id(chr_id)
    'chr2R'
    >>> chr_id = "3L"
    >>> check_convert_chr_id(chr_id)
    'chr3L'
    >>> chr_id = "4L"
    >>> check_convert_chr_id(chr_id)
    False
    >>> chr_id = "chrM"
    >>> check_convert_chr_id(chr_id, id_style=2)
    'MT'
    >>> chr_id = "chr2R"
    >>> check_convert_chr_id(chr_id, id_style=2)
    '2R'
    >>> chr_id = "5"
    >>> check_convert_chr_id(chr_id, id_style=2)
    '5'
    >>> chr_id = "chrA"
    >>> check_convert_chr_id(chr_id, id_style=2)
    False
    >>> chr_id = "chrA"
    >>> check_convert_chr_id(chr_id, id_style=0)
    'chrA'


    """
    assert chr_id, "given chr_id empty"

    if not id_style: # If id_style == 0 or False.
        return chr_id

    elif id_style == 1:
        if re.search("^chr", chr_id):
            if chr_id in add_chr_names_dic or re.search("^chr[\dMXY]+$", chr_id):
                return chr_id
            else:
                return False
        else:
            # Convert to "chr" IDs.
            if chr_id == "MT": # special case MT -> chrM.
                return "chrM"
            if chr_id in add_chr_names_dic or re.search("^[\dXY]+$", chr_id):
                return "chr" + chr_id
            else:
                return False

    elif id_style == 2:

        if re.search("^chr", chr_id):
            if chr_id == "chrM": # special case chrM -> MT.
                return "MT"
            if chr_id in add_chr_names_dic or re.search("^chr[\dXY]+$", chr_id):
                # Cut out chr suffix.
                m = re.search("chr(.+)", chr_id)
                assert m, "no match for regex search"
                chr_suffix = m.group(1)
                return chr_suffix
            else:
                return False

        else:
            if chr_id == "MT": # special case MT.
                return chr_id
            if chr_id in add_chr_names_dic or re.search("^[\dXY]+$", chr_id):
                return chr_id
            else:
                return False
    else:
        assert False, "invalid id_style set"


################################################################################

def bed_check_six_col_format(bed_file):
    """
    Check whether given .bed file has 6 columns.

    >>> test_bed = "test_data/test1.bed"
    >>> bed_check_six_col_format(test_bed)
    True
    >>> test_bed = "test_data/empty_file"
    >>> bed_check_six_col_format(test_bed)
    False

    """

    six_col_format = False
    with open(bed_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) == 6:
                six_col_format = True
            break
    f.closed
    return six_col_format


################################################################################

def bed_get_region_ids(bed_file,
                       check=True):
    """
    Read in .bed file, return region/site IDs (column 5 IDs).

    >>> test_file = "test_data/test3.bed"
    >>> bed_get_region_ids(test_file)
    {'CLIP1': 1, 'CLIP2': 1}

    """
    ids_dic = {}
    with open(bed_file) as f:
        for line in f:
            row = line.strip()
            cols = line.strip().split("\t")
            site_id = cols[3]
            assert site_id not in ids_dic, "column 4 IDs not unique in given .bed file \"%s\"" %(bed_file)
            ids_dic[site_id] = 1
    f.closed
    if check:
        assert ids_dic, "No IDs read into dictionary (input file \"%s\" empty or malformatted?)" % (bed_file)
    return ids_dic


################################################################################

def bed_check_unique_col4_ids(bed_file):
    """
    Check whether .bed file (6 column format with IDs in column 4)
    has unique column 4 IDs.

    >>> test_bed = "test_data/test1.bed"
    >>> bed_check_unique_ids(test_bed)
    True
    >>> test_bed = "test_data/test2.bed"
    >>> bed_check_unique_ids(test_bed)
    False

    """
    ids_dic = {}
    check = True
    with open(ids_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            site_id = cols[3]
            if site_id not in ids_dic:
                ids_dic[site_id] = 1
            else:
                check = False
    f.closed
    assert ids_dic, "IDs dictionary ids_dic empty"
    return check


################################################################################

def count_file_rows(in_file,
                    nr_cols=False):
    """
    Count number of file rows. If nr_cols set, demand certain (nr_cols) number
    of columns (separated by tab), in order for row to be counted.

    >>> test_file = "test_data/test1.bed"
    >>> count_file_rows(test_file)
    7
    >>> test_file = "test_data/empty_file"
    >>> count_file_rows(test_file)
    0

    """
    c = 0
    with open(in_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            if nr_cols:
                if len(cols) == nr_cols:
                    c += 1
            else:
                c += 1
    f.closed
    return c


################################################################################

def samtools_extract_reads_from_bed_regions(in_bed, in_bam, out_bam):
    """
    Get BAM reads from BED regions, store in new BAM file.

    OLD: samtools view -L ignores strandness info in BED files, will output
    overlapping reads on both strands.

    """
    assert os.path.exists(in_bed), "in_bed does not exist"
    assert os.path.exists(in_bam), "in_bam does not exist"

    check_cmd = "samtools view -L " + in_bed + " " + in_bam + " -o " + out_bam
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "samtools has problems with your input:\n%s\n%s" %(check_cmd, output)


################################################################################

def bam_extract_reads_from_bed_regions(in_bed, in_bam, out_bam,
                                       no_name_check=False,
                                       reverse_strand=False,
                                       sort_bed=False):
    """
    Get BAM reads from BED regions, store in new BAM file.

    Using intersectBed instead of samtools view -L, which allows -s for strand
    specific filtering.

    no_name_check:
        Activate intersectBed -nonamecheck
    reverse_strand:
        If gene reads are on reverse strand.

    """
    assert os.path.exists(in_bed), "in_bed does not exist"
    assert os.path.exists(in_bam), "in_bam does not exist"

    strand_set = "-s"
    if reverse_strand:
        strand_set = "-S"
    if no_name_check:
        strand_set += " -nonamecheck"

    if sort_bed:
        check_cmd = "sort -k1,1 -k2,2n " + in_bed + " | " + "intersectBed -abam " + in_bam + " -b stdin " + strand_set + " -sorted > " + out_bam
    else:
        check_cmd = "intersectBed -abam " + in_bam + " -b " + in_bed + " " + strand_set + " > " + out_bam
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "intersectBed has problems with your input:\n%s\n%s" %(check_cmd, output)


################################################################################

def bam_to_bed_get_isr_stats(in_bam, next_ol_bed,
                             isr_bed=False,
                             isr_ext_mode=1,
                             isr_max_reg_len=10,
                             reverse_strand=False,
                             nexts_cisrc_dic=False,
                             isr_intron_reg_dic=None,
                             new_is_ids=True):
    """

    Get intron-spanning reads (ISR), overlap their matched regions
    (isr_ext_mode=1) or ends (isr_ext_mode=2) with NEXT regions.
    Also extract intronic regions with ISR.

    in_bam:
        BAM file with gene region reads containing exonic sites.
    next_ol_bed:
        Exons of transcripts with NEXT overlapping exonic sites BED.
    isr_ext_mode:
        Extraction mode for IS read regions.
        1: take whole match regions.
        2: take end / start positions of IS read matches.
    isr_max_reg_len:
        Maximum length of IS read start end regions. If match is >
        isr_max_reg_len, use isr_max_reg_len as length of the region.
        If isr_max_reg_len=False, use full length of match as region
        length (if isr_ext_mode == 1).
    reverse_strand:
        Reads mapping to reverse strand (e.g. for certain RNA-seq datasets).
    isr_intron_reg_dic:
        Intronic region -> ISR count.
    isr_bed:
        Define BED file to store ISR reads in.

    Removed:
        next2isrc_dic:
            NEXT to overlapping IS read count.

    """
    assert os.path.exists(in_bam), "in_bam does not exist"
    assert os.path.exists(next_ol_bed), "next_ol_bed does not exist"

    if isr_bed:
        tmp_bed = isr_bed
    else:
        random_id = uuid.uuid1()
        tmp_bed = str(random_id) + ".bam_next_ol.tmp.bed"

    if not nexts_cisrc_dic:
        nexts_cisrc_dic = {} # NEXTs connecting IS read count.

    check_cmd = "bamToBed -i " + in_bam + " -bed12"
    output = subprocess.getoutput(check_cmd)

    c_ir_reads = 0
    ISRBEDOUT = open(tmp_bed, "w")
    for line in output.split('\n'):
        cols = line.strip().split("\t")
        # Only split reads.
        if cols[9] == "1":
            continue
        chr_id = cols[0]
        reg_s = int(cols[1])
        reg_e = int(cols[2])
        read_id = cols[3]
        read_pol = cols[5]
        if reverse_strand:
            read_pol = get_rev_strand(read_pol)
        # assert re.search(",", cols[10]), "-bed 12 column 11 of split read is missing \",\" in line \"%s\"" %(line)
        # assert re.search(",", cols[11]), "-bed 12 column 12 of split read is missing \",\" in line \"%s\"" %(line)

        """
        Output end positions (length 1) of the intron-spanning reads (i.e.,
        the last positions overlapping with exon(s)).

        chr1	1542166	1542167	isr_1	0	-
        chr1	1543868	1543869	isr_1	0	-
        chr1	1542166	1542167	isr_2	0	-
        chr1	1543868	1543869	isr_2	0	-
        ...

        In case of split reads (2 or more splits), output each intron split
        as separate IR read. Usually we should have cols[9] == 2, but for
        RNA-seq data (depending on type) we can also have more.

        """
        parts_list = cols[10].split(",")
        offsets_list = cols[11].split(",")

        """
        chr1	10168245	10171147	GACGG:HWI-D00611:119:C6K7PANXX:5:1316:4890:22830/2	255	+	10168245	10171147	255,0,0	2	25,10	0,2892
        chr1	10168259	10171148	AAAAC:HWI-D00611:119:C6K7PANXX:5:1103:13861:44307/2	255	+	10168259	10171148	255,0,0	2	11,11	0,2878
        chr1	10168270	10168313	ATTGA:HWI-D00611:119:C6K7PANXX:5:1315:5248:33314/2	255	+	10168270	10168313	255,0,0	1	43	0
        l: 43

        """

        for i in range(len(parts_list)-1):
            l_p1 = int(parts_list[i])
            l_p2 = int(parts_list[i+1])
            os1 = int(offsets_list[i])
            os2 = int(offsets_list[i+1])
            p1_e = reg_s + l_p1 + os1
            # p1_s = p1_e - 1
            p1_s = p1_e - 1
            p2_s = reg_s + os2
            # p2_e = p2_s + 1
            p2_e = p2_s + 1

            """
            Case: take full match for coverage calculations.
            This also influences exon border site discovery, as sites
            do not need to end exactly at exon ends anymore.

            """
            if isr_ext_mode == 1:
                p1_s = p1_e - l_p1
                p2_e = p2_s + l_p2
                if isr_max_reg_len:
                    if l_p1 > isr_max_reg_len:
                        p1_s = p1_e - isr_max_reg_len
                    if l_p2 > isr_max_reg_len:
                        p2_e = p2_s + isr_max_reg_len

            c_ir_reads += 1
            new_read_id = read_id
            if new_is_ids:
                new_read_id = "isr_%i" %(c_ir_reads)

            #isr_119637
            #chr11:7995358-7995367
            #chr11	7995357	7995367	isr_119637	0	+
            #chr11	7995944	7995947	isr_119637	0	+

            if isr_intron_reg_dic is not None:
                intron_s = p1_e
                intron_e = p2_s
                intron_reg = "%s,%i,%i,%s" %(chr_id, intron_s, intron_e, read_pol)
                if intron_reg in isr_intron_reg_dic:
                    isr_intron_reg_dic[intron_reg] += 1
                else:
                    isr_intron_reg_dic[intron_reg] = 1

            ISRBEDOUT.write("%s\t%i\t%i\t%s\t0\t%s\n" %(chr_id, p1_s, p1_e, new_read_id, read_pol))
            ISRBEDOUT.write("%s\t%i\t%i\t%s\t0\t%s\n" %(chr_id, p2_s, p2_e, new_read_id, read_pol))

    ISRBEDOUT.close()

    #
    #     l_p1 = int(parts[0])
    #     l_p2 = int(parts[1])
    #     # Get IS read mapped start+end positions.
    #     if is_ext_mode == 1:
    #         p1e = reg_s + l_p1
    #         p1s = p1e - 1
    #         p2s = reg_e - l_p2
    #         p2e = p2s + 1
    #     elif is_ext_mode == 2:
    #         p1s = reg_s
    #         p1e = reg_s + l_p1
    #         p2s = reg_e - l_p2
    #         p2e = reg_e
    #     else:
    #         assert False, "invalid is_ext_mode given"
    #     new_read_id = read_id
    #     if new_is_ids:
    #         new_read_id = "isr_%i" %(c_ir_reads)
    #     ISRBEDOUT.write("%s\t%i\t%i\t%s\t0\t%s\n" %(chr_id, p1s, p1e, new_read_id, read_pol))
    #     ISRBEDOUT.write("%s\t%i\t%i\t%s\t0\t%s\n" %(chr_id, p2s, p2e, new_read_id, read_pol))
    # ISRBEDOUT.close()

    # If no IS reads, no need to look at overlaps with NEXT exon ends.
    if not c_ir_reads:
        return nexts_cisrc_dic

    check_cmd = "intersectBed -a " + next_ol_bed + " -b " + tmp_bed + " -s -wb"

    output = subprocess.getoutput(check_cmd)

    if not isr_bed:
        if os.path.exists(tmp_bed):
            os.remove(tmp_bed)

    # Again if no overlap, return empty dics.
    if not output:
        return nexts_cisrc_dic

    isr2next_list_dic = {}
    for line in output.split('\n'):
        cols = line.strip().split("\t")
        next_id = cols[3]
        isr_id = cols[9]
        # if next_id in next2isrc_dic:
        #     next2isrc_dic[next_id] += 1
        # else:
        #     next2isrc_dic[next_id] = 1
        if isr_id in isr2next_list_dic:
            isr2next_list_dic[isr_id].append(next_id)
        else:
            isr2next_list_dic[isr_id] = [next_id]

    for isr_id in isr2next_list_dic:
        next_pairs_seen_dic = {}
        for next1 in isr2next_list_dic[isr_id]:
            for next2 in isr2next_list_dic[isr_id]:
                if next1 == next2:
                    continue
                nexts1 = "%s,%s" %(next1, next2)
                nexts2 = "%s,%s" %(next2, next1)
                if nexts1 in next_pairs_seen_dic:
                    continue
                if nexts2 in next_pairs_seen_dic:
                    continue

                con_id = "%s,%s" %(next1, next2)
                if con_id in nexts_cisrc_dic:
                    nexts_cisrc_dic[con_id] += 1
                else:
                    nexts_cisrc_dic[con_id] = 1
                next_pairs_seen_dic[nexts1] = 1
                next_pairs_seen_dic[nexts2] = 1

    return nexts_cisrc_dic


################################################################################

def get_isolated_transcripts(single_ex_tr_dic, tr2reg_dic,
                             tmp_out_folder=False):

    """
    Get isolated transcripts (no overlap with other transcripts).
    Return isolated transcript IDs dictionary.

    single_ex_tr_dic:
        Single exon transcript IDs dictionary
    tr2reg_dic:
        transcript ID -> genomic region [chr_id, tr_s, tr_e, gene_pol]

    >>> tr2reg_dic = {'t1': ['chr1', 1000, 2000, '+'], 't2': ['chr1', 1800, 2800, '+'], 't3': ['chr1', 3000, 4000, '+'], 't4': ['chr1', 5000, 6000, '+']}
    >>> single_ex_tr_dic = {'t1': 1, 't2': 1, 't3': 1}
    >>> get_isolated_transcripts(single_ex_tr_dic, tr2reg_dic)
    {'t3': 1}

    """

    assert single_ex_tr_dic, "single_ex_tr_dic empty"

    random_id = uuid.uuid1()
    tmp_bed = str(random_id) + ".gene_regions.tmp.bed"
    if tmp_out_folder:
        tmp_bed = tmp_out_folder + "/" + tmp_bed

    TREGOUT = open(tmp_bed, "w")
    for tr_id in single_ex_tr_dic:
        chr_id = tr2reg_dic[tr_id][0]
        tr_s = tr2reg_dic[tr_id][1]
        tr_e = tr2reg_dic[tr_id][2]
        tr_pol = tr2reg_dic[tr_id][3]
        TREGOUT.write("%s\t%i\t%i\t%s\t0\t%s\n" %(chr_id, tr_s, tr_e, tr_id, tr_pol))
    TREGOUT.close()

    # Overlap calculation with itself.
    params = " -s"
    check_cmd = "intersectBed -a " + tmp_bed + " -b " + tmp_bed + " " + params
    output = subprocess.getoutput(check_cmd)

    """
    $ cat genes.bed
    chr1	1000	2000	g1	0	+
    chr1	3000	4000	g2	0	+
    chr1	3800	4800	g3	0	+
    chr1	4600	5600	g4	0	+
    chr1	6000	7000	g5	0	+

    $ intersectBed -a genes.bed -b genes.bed -s
    chr1	1000	2000	g1	0	+
    chr1	3000	4000	g2	0	+
    chr1	3800	4000	g2	0	+
    chr1	3800	4000	g3	0	+
    chr1	3800	4800	g3	0	+
    chr1	4600	4800	g3	0	+
    chr1	4600	4800	g4	0	+
    chr1	4600	5600	g4	0	+
    chr1	6000	7000	g5	0	+

    So count column 4 ID appearances. Appearances == 1 means isolated (only
    overlapping with itself).

    """
    ol_tr_ids_dic = {}

    for line in output.split('\n'):
        cols = line.strip().split("\t")
        tr_id = cols[3]
        if tr_id in ol_tr_ids_dic:
            ol_tr_ids_dic[tr_id] += 1
        else:
            ol_tr_ids_dic[tr_id] = 1

    isolated_tr_ids_dic = {}
    for tr_id in ol_tr_ids_dic:
        if ol_tr_ids_dic[tr_id] == 1:
            isolated_tr_ids_dic[tr_id] = 1

    if os.path.exists(tmp_bed):
        os.remove(tmp_bed)

    return isolated_tr_ids_dic


################################################################################

def remove_overlapping_genes(gid2sel_tr_dic, gid2isrc_dic, tr2reg_dic,
                             remove_single_ex_genes=False,
                             trid2exc_dic=False,
                             tmp_out_folder=False,
                             min_isrc=2):
    """
    Overlap gene regions (using longest transcript of each gene) with
    each other, and remove gene A, if it fully overlaps with gene B
    and:
        1) gene A does not have intron-spanning reads, while
        gene B has. If both do not have, keep both. Minimum ISR count
        == min_isrc.
        This will also remove single exon genes, if the longer gene
        fully covers it and has > min_isrc ISR reads.
        Alternatively, set remove_single_ex_genes to True, to remove
        fully overlapping single exon genes in any case.

    gid2sel_tr_dic:
        Gene ID -> selected transcript ID (default: longest one).
    gid2isrc_dic:
        Intron-spanning read count of all transcripts belonging to
        gene, so:
        gene ID -> ISR count (ISR sum of all gene transcripts)
    tr2reg_dic:
        transcript ID -> genomic region [chr_id, tr_s, tr_e, gene_pol]
    min_isrc:
        Minimum ISR count used for filtering.

    >>> gid2sel_tr_dic = {'g1': 't1', 'g2': 't2', 'g3': 't3', 'g4': 't4'}
    >>> gid2isrc_dic = {'g1': 10, 'g2': 0, 'g3': 2, 'g4': 0}
    >>> tr2reg_dic = {'t1': ['chr1', 1000, 2000, '+'], 't2': ['chr1', 1200, 1600, '+'], 't3': ['chr1', 1400, 1800, '+'], 't4': ['chr1', 1000, 2000, '-']}
    >>> remove_overlapping_genes(gid2sel_tr_dic, gid2isrc_dic, tr2reg_dic)
    {'g2': 1}
    >>> gid2sel_tr_dic = {'g5': 't5', 'g6': 't6', 'g7': 't7'}
    >>> gid2isrc_dic = {'g5': 0, 'g6': 15, 'g7': 0}
    >>> tr2reg_dic = {'t5': ['chr1', 5000, 6000, '+'], 't6': ['chr1', 5000, 6000, '+'], 't7': ['chr1', 7000, 8000, '+']}
    >>> remove_overlapping_genes(gid2sel_tr_dic, gid2isrc_dic, tr2reg_dic)
    {'g5': 1}
    >>> gid2sel_tr_dic = {'g1': 't1', 'g2': 't2'}
    >>> gid2isrc_dic = {'g1': 0, 'g2': 0}
    >>> tr2reg_dic = {'t1': ['chr1', 1000, 2000, '+'], 't2': ['chr1', 1600, 1800, '+']}
    >>> trid2exc_dic = {'t1': 5, 't2': 1}
    >>> remove_overlapping_genes(gid2sel_tr_dic, gid2isrc_dic, tr2reg_dic)
    {}
    >>> remove_overlapping_genes(gid2sel_tr_dic, gid2isrc_dic, tr2reg_dic, remove_single_ex_genes=True, trid2exc_dic=trid2exc_dic)
    {'g2': 1}

    """

    random_id = uuid.uuid1()
    tmp_bed = str(random_id) + ".gene_regions.tmp.bed"
    if tmp_out_folder:
        tmp_bed = tmp_out_folder + "/" + tmp_bed

    # Write gene regions BED for overlap calculation.
    GREGOUT = open(tmp_bed, "w")
    gid2len_dic = {}
    for gid in gid2sel_tr_dic:
        tr_id = gid2sel_tr_dic[gid]
        chr_id = tr2reg_dic[tr_id][0]
        tr_s = tr2reg_dic[tr_id][1]
        tr_e = tr2reg_dic[tr_id][2]
        tr_pol = tr2reg_dic[tr_id][3]
        tr_sc = gid2isrc_dic[gid]
        gid2len_dic[gid] = tr_e - tr_s
        GREGOUT.write("%s\t%i\t%i\t%s\t%i\t%s\n" %(chr_id, tr_s, tr_e, gid, tr_sc, tr_pol))
    GREGOUT.close()

    # Overlap calculation with itself.
    params = " -s -F 1.0 -wao"
    check_cmd = "intersectBed -a " + tmp_bed + " -b " + tmp_bed + " " + params
    output = subprocess.getoutput(check_cmd)

    """
    -F 1.0
    Only full overlaps (fraction of B), i.e. shorter genes get reported
    (and self overlaps).

    Example:
    $ intersectBed -a a.bed -b a.bed -s -F 1.0 -wao
    chr1	1000	2000	g1	10	+	chr1	1000	2000	g1	10	+	1000
    chr1	1000	2000	g1	10	+	chr1	1200	1600	g2	0	+	400
    chr1	1000	2000	g1	10	+	chr1	1400	1800	g3	2	+	400
    chr1	1200	1600	g2	0	+	chr1	1200	1600	g2	0	+	400
    chr1	1400	1800	g3	2	+	chr1	1400	1800	g3	2	+	400

    """
    gid2remove_dic = {}

    for line in output.split('\n'):
        cols = line.strip().split("\t")
        gid1 = cols[3]
        gid2 = cols[9]
        # overlap_len = int(12)
        if gid1 == gid2:
            continue
        # Gene ID 1 should be longer gene, check.
        gid1_len = gid2len_dic[gid1]
        gid2_len = gid2len_dic[gid2]
        assert gid1_len >= gid2_len, "reported gene ID 1 (%s) length < gene ID 2 (%s) length (%i < %i)" %(gid1, gid2, gid1_len, gid2_len)
        # Compare and remove if criterion satisfied.
        gid1_isrc = gid2isrc_dic[gid1]
        gid2_isrc = gid2isrc_dic[gid2]
        if gid2_isrc == 0 and gid1_isrc >= min_isrc:
            gid2remove_dic[gid2] = 1
        # Remove overlapping single exon genes no matter what.
        if remove_single_ex_genes:
            tid2 = gid2sel_tr_dic[gid2]
            assert trid2exc_dic, "remove_single_ex_genes=True requires trid2exc_dic"
            assert trid2exc_dic[tid2], "transcript ID %s not in trid2exc_dic" %(tid2)
            tid2_exc = trid2exc_dic[tid2]
            if tid2_exc == 1:
                gid2remove_dic[gid2] = 1

    if os.path.exists(tmp_bed):
        os.remove(tmp_bed)

    return gid2remove_dic


################################################################################

def get_trid_isrc_full_con(tr_id, tr_exc, exid2next_dic, nexts_cisrc_dic):
    """

    Get intron-spanning read count for transcript with ID tr_id. Also return
    if transcript exons are fully connected by intron-spanning reads
    (True, False).

    tr_id:
        Transcript ID
    tr_exc:
        Transcript exon count.
    exid2next_dic:
        exon ID to NEXT ID mapping
    nexts_cisrc_dic:
        Connected NEXT IDs with format "NEXT1,NEXT2" and mapping:
        "NEXT1,NEXT2" -> connecting IS read count

    >>> exid2next_dic = {"t1_e1" : "NEXT1", "t1_e2" : "NEXT2", "t1_e3": "NEXT3", "t2_e1": "NEXT4"}
    >>> nexts_cisrc_dic = {"NEXT1,NEXT2": 4, "NEXT2,NEXT3": 2}
    >>> get_trid_isrc_full_con("t1", 3, exid2next_dic, nexts_cisrc_dic)
    (6, True)
    >>> nexts_cisrc_dic = {"NEXT2,NEXT3": 5}
    >>> get_trid_isrc_full_con("t1", 3, exid2next_dic, nexts_cisrc_dic)
    (5, False)
    >>> get_trid_isrc_full_con("t2", 1, exid2next_dic, nexts_cisrc_dic)
    (0, False)

    """

    if tr_exc == 1:
        return 0, False

    isr_c = 0
    fully_con = True

    # print("tr_id:", tr_id)

    for i in range(tr_exc-1):
        ex_nr1 = i + 1
        ex_nr2 = i + 2
        ex_id1 = tr_id + "_e%i" %(ex_nr1)
        ex_id2 = tr_id + "_e%i" %(ex_nr2)
        # print(ex_id1, ex_id2)
        next1 = exid2next_dic[ex_id1]
        next2 = exid2next_dic[ex_id2]
        nexts1 = next1 + "," + next2
        nexts2 = next2 + "," + next1
        nexts1_yes = False
        nexts2_yes = False
        isr_c_pair = 0
        if nexts1 in nexts_cisrc_dic:
            nexts1_yes = True
        if nexts2 in nexts_cisrc_dic:
            nexts2_yes = True
        if nexts1_yes and nexts2_yes:
            assert False, "NEXT ID combination appears twice in nexts_cisrc_dic (%s, %s)" %(nexts1, nexts2)
        if nexts1_yes:
            isr_c_pair = nexts_cisrc_dic[nexts1]
        elif nexts2_yes:
            isr_c_pair = nexts_cisrc_dic[nexts2]
        isr_c += isr_c_pair
        if not nexts1_yes and not nexts2_yes:
            fully_con = False
        # print("isrc:", isr_c_pair)

    return isr_c, fully_con


################################################################################

def get_rev_strand(strand):
    """
    Get reverse strand.

    >>> get_rev_strand("-")
    '+'

    """
    if strand == "+":
        return "-"
    elif strand == "-":
        return "+"
    else:
        assert False, "invalid strand information given (%s)" %(strand)


################################################################################

def check_tr_id_full_coverage(tr_id, trid2exc_dic, regid2nc_dic,
                              pseudo_counts=False):
    """
    Check if each exon of a given transcript ID is covered by > 0 reads.
    If so return True, else False.

    >>> trid2exc_dic = {"t1" : 2, "t2" : 2, "t3" : 1}
    >>> regid2nc_dic = {"t1_e1": 0.4, "t1_e2": 1.2, "t2_e1": 0.4, "t2_e2": 0.0, "t3_e1": 1.6}
    >>> check_tr_id_full_coverage("t1", trid2exc_dic, regid2nc_dic)
    True
    >>> check_tr_id_full_coverage("t2", trid2exc_dic, regid2nc_dic)
    False
    >>> check_tr_id_full_coverage("t3", trid2exc_dic, regid2nc_dic)
    True

    """
    min_ex_cov = 0
    if pseudo_counts:
        min_ex_cov = 1
    for i in range(trid2exc_dic[tr_id]):
        ex_nr = i + 1
        ex_id = tr_id + "_e" + str(ex_nr)
        ex_cov = regid2nc_dic[ex_id]
        if ex_cov == min_ex_cov:
            return False
    return True


################################################################################

def get_sid_trid_combination_score(site_id, tr_ids_list, idfilt2best_trids_dic):
    """
    Get site ID - transcript ID combination score, based on selected
    transcripts for each of the 10 different filter settings.

    10 transcript quality filter settings:
    EIR
    EXB
    TSC
    ISRN
    ISR
    ISRFC
    SEO
    FUCO
    TCOV
    TSL

    idfilt2best_trids_dic:
        "site_id,filter_id"
        -> top transcript ID(s) after applying filter on exon IDs > min_eir

    >>> site_id = "s1"
    >>> idfilt2best_trids_dic = {"s1,EIR" : ["t1"], "s1,EXB" : ["t1"], "s1,TSC" : ["t1"], "s1,ISRN" : ["t1"], "s1,ISR" : ["t1"], "s1,ISRFC" : ["t1"], "s1,SEO" : ["t1"], "s1,FUCO" : ["t1"], "s1,TCOV" : ["t1"], "s1,TSL" : ["t1"]}
    >>> tr_ids_list = ["t1"]
    >>> get_sid_trid_combination_score(site_id, tr_ids_list, idfilt2best_trids_dic)
    {'t1': 10}
    >>> idfilt2best_trids_dic = {"s1,EIR" : ["t1", "t2"], "s1,EXB" : ["t1", "t2"], "s1,TSC" : ["t1"], "s1,ISRN" : ["t2"], "s1,ISR" : ["t1"], "s1,ISRFC" : ["t1"], "s1,SEO" : ["t1"], "s1,FUCO" : ["t1", "t2"], "s1,TCOV" : ["t1"], "s1,TSL" : ["t2"]}
    >>> tr_ids_list = ["t1", "t2", "t3"]
    >>> get_sid_trid_combination_score(site_id, tr_ids_list, idfilt2best_trids_dic)
    {'t1': 8, 't2': 5, 't3': 0}

    """
    assert tr_ids_list, "tr_ids_list empty"
    filter_ids = ["EIR", "EXB", "TSC", "ISRN", "ISR", "ISRFC", "SEO", "FUCO", "TCOV", "TSL"]
    trid2comb_sc_dic = {}
    for tr_id in tr_ids_list:
        trid2comb_sc_dic[tr_id] = 0
    for tr_id in tr_ids_list:
        for fid in filter_ids:
            sitefiltid = "%s,%s" %(site_id, fid)
            if tr_id in idfilt2best_trids_dic[sitefiltid]:
                trid2comb_sc_dic[tr_id] += 1
    return trid2comb_sc_dic


################################################################################

def list_found_duplicates(in_list):
    """
    Check list for duplicate entries. Return True if duplicates found,
    and False if not duplicates found.

    >>> in_list = ["hallo", "hello"]
    >>> list_found_duplicates(in_list)
    False
    >>> in_list = ["hallo", "hello", "hollo", "hello"]
    >>> list_found_duplicates(in_list)
    True

    """
    if len(set(in_list)) == len(in_list):
        return False
    else:
        return True


################################################################################

def get_exid_isr_hood_count(exon_id, exid2trid_dic,
                                exid2next_dic, nexts_cisrc_dic):
    """
    Given an exon ID, get intron-spanning read count that connects it to
    its neighboring exons. Return count.

    >>> exid2trid_dic = {"t1_e1" : "t1", "t1_e2" : "t1", "t2_e1" : "t2", "t2_e2" : "t2", "t2_e3" : "t2", "t3_e1" : "t3"}
    >>> exid2next_dic = {"t1_e1" : "n1", "t1_e2" : "n2", "t2_e1" : "n3", "t2_e2" : "n2", "t2_e3" : "n4", "t3_e1" : "n5"}
    >>> nexts_cisrc_dic = {"n1,n2" : 10, "n2,n3" : 10, "n2,n3" : 8, "n2,n4" : 5}
    >>> get_exid_isr_hood_count("t1_e2", exid2trid_dic, exid2next_dic, nexts_cisrc_dic)
    10
    >>> get_exid_isr_hood_count("t2_e2", exid2trid_dic, exid2next_dic, nexts_cisrc_dic)
    13
    >>> get_exid_isr_hood_count("t3_e1", exid2trid_dic, exid2next_dic, nexts_cisrc_dic)
    0

    """
    isrn_c = 0
    assert re.search(".+_e\d", exon_id), "exon ID %s has invalid format"
    m = re.search("(.+)_e(\d+)", exon_id)
    tr_id = m.group(1)
    exon_nr = int(m.group(2))
    next = exid2next_dic[exon_id]
    # Neighbor exons.
    us_exon_id = tr_id + "_e" + str(exon_nr-1)
    ds_exon_id = tr_id + "_e" + str(exon_nr+1)
    if us_exon_id in exid2next_dic:
        us_next = exid2next_dic[us_exon_id]
        nexts1 = "%s,%s" %(next,us_next)
        nexts2 = "%s,%s" %(us_next,next)
        if nexts1 in nexts_cisrc_dic:
            isrn_c += nexts_cisrc_dic[nexts1]
        elif nexts2 in nexts_cisrc_dic:
            isrn_c += nexts_cisrc_dic[nexts2]
    if ds_exon_id in exid2next_dic:
        ds_next = exid2next_dic[ds_exon_id]
        nexts1 = "%s,%s" %(next,ds_next)
        nexts2 = "%s,%s" %(ds_next,next)
        if nexts1 in nexts_cisrc_dic:
            isrn_c += nexts_cisrc_dic[nexts1]
        elif nexts2 in nexts_cisrc_dic:
            isrn_c += nexts_cisrc_dic[nexts2]
    return isrn_c


################################################################################

def read_in_12col_bed(in_bed, bed_row_dic=False):
    """
    Read in 12-column BED (in_bed), store in bed_row_dic.

    Region ID in column BED 4:
    chr14	23321288	23326185	ENST00000216727.8	0	+	23321288 ...

    """

    if not bed_row_dic:
        bed_row_dic = {}

    with open(in_bed) as f:
        for line in f:
            row = line.strip()
            cols = line.strip().split("\t")
            reg_id = cols[3]
            bed_row_dic[reg_id] = row
    f.closed
    return bed_row_dic


################################################################################

def write_12col_bed(bed_row_dic, out_bed):
    """
    Write BED row dictionary to BED file.

    """

    assert bed_row_dic, "given bed_row_dic empty"

    OUTBED = open(out_bed, "w")
    c_out = 0
    for reg_id in bed_row_dic:
        bed_row = bed_row_dic[reg_id]
        OUTBED.write("%s\n" %(bed_row))
        c_out += 1
    OUTBED.close()
    assert c_out, "nothing output by write_12col_bed()"


################################################################################

def get_exon_border_site_pairs(sites_bed, isr_bed,
                               max_site_dist=10,
                               id2gen_se_dic=False,
                               id2next_list_dic=False):
    """
    Get exon border site pairs, based on intron-spanning reads connecting
    the two sites. If for a given site > 1 connection is supported by
    intron-spanning reads, also save this. Return two dictionaries:
    site_id --> [connected_site_id1, connected_site_id2, ...]
    "site_id1,site_id2" --> # of connecting reads.

    id2gen_se_dic:
        site_id -> [gen_start, gen_end]
    id2next_list_dic:
        Use site ID -> NEXT list mapping, to remove border pairs on same
        NEXT exon regions (can happen if --isr-ext-mode 2)
        Also use idnext2ol_se_dic to only remove sites nearby from
        id2ids_dic.

    test_exb_sites.bed
    chr1	1090	2000	s1	0	+
    chr1	3000	3010	s2	0	+
    chr1	4000	4020	s3	0	+
    chr1	5000	5020	s4	0	+

    test_exb_isr.bed
    chr1	1999	2000	isr1	0	+
    chr1	1999	2000	isr2	0	+
    chr1	1999	2000	isr3	0	+
    chr1	1999	2000	isr4	0	+
    chr1	3000	3001	isr1	0	+
    chr1	3000	3001	isr2	0	+
    chr1	3000	3001	isr3	0	+
    chr1	4000	4001	isr4	0	+

    Use "-EB-" to separate site IDs and form new ID:
    id1-EB-id2

    Resulting isr2site_id_list_dic:
    {'isr1': ['s1', 's2'], 'isr2': ['s1', 's2'], 'isr3': ['s1', 's2'],
    'isr4': ['s1', 's3']}

    >>> sites_bed = "test_data/test_exb_sites.bed"
    >>> isr_bed = "test_data/test_exb_isr.bed"
    >>> get_exon_border_site_pairs(sites_bed, isr_bed)
    ({'s1': ['s2', 's3'], 's2': ['s1'], 's3': ['s1']}, {'s1-EB-s2': 3, 's2-EB-s1': 3, 's1-EB-s3': 1, 's3-EB-s1': 1})
    >>> sites_bed = "test_data/test_exb_sites2.bed"
    >>> get_exon_border_site_pairs(sites_bed, isr_bed)
    ({}, {})

    """

    id2ids_dic = {}
    ids2isrc_dic = {}

    assert os.path.exists(sites_bed), "%s transcript context genomic BED file not found" %(sites_bed)
    assert os.path.exists(isr_bed), "%s intron-spanning reads BED file not found" %(isr_bed)

    check_cmd = "intersectBed -a " + sites_bed + " -b " + isr_bed + " -s -wb"
    output = subprocess.getoutput(check_cmd)

    # Again if no overlap, return empty dics.
    if not output:
        return id2ids_dic, ids2isrc_dic

    # Connect intron-spanning read ID to site IDs.
    isr2site_id_list_dic = {}
    for line in output.split('\n'):
        cols = line.strip().split("\t")
        site_id = cols[3]
        isr_id = cols[9]
        if isr_id in isr2site_id_list_dic:
            isr2site_id_list_dic[isr_id].append(site_id)
        else:
            isr2site_id_list_dic[isr_id] = [site_id]

    """
    Resulting isr2site_id_list_dic:
    {'isr1': ['s1', 's2'], 'isr2': ['s1', 's2'], 'isr3': ['s1', 's2'],
    'isr4': ['s1', 's3']}
    """

    # Site IDs to connecting intron-spanning read count.
    for isr_id in isr2site_id_list_dic:
        for id1 in isr2site_id_list_dic[isr_id]:
            for id2 in isr2site_id_list_dic[isr_id]:
                if id1 == id2:
                    continue

                if id1 in id2ids_dic:
                    id2ids_dic[id1].append(id2)
                else:
                    id2ids_dic[id1] = [id2]

                ids = "%s-EB-%s" %(id1, id2)
                if ids in ids2isrc_dic:
                    ids2isrc_dic[ids] += 1
                else:
                    ids2isrc_dic[ids] = 1

    rem_sids_dic = {}

    for site_id in id2ids_dic:

        # Make lists non-redundant.
        id2ids_dic[site_id] = list(set(id2ids_dic[site_id]))
        id2ids_dic[site_id].sort()

        # Remove connections on same NEXT exons + in close distance.
        if id2gen_se_dic and id2next_list_dic:
            sid_next_list = id2next_list_dic[site_id]
            sid_s = id2gen_se_dic[site_id][0] # site ID genomic start.
            sid_e = id2gen_se_dic[site_id][1] # site ID genomic end.
            rem_assoc_sids_dic = {}
            for sid2 in id2ids_dic[site_id]:
                sid2_next_list = id2next_list_dic[sid2]
                for sid2_next in sid2_next_list:
                    if sid2_next in sid_next_list:
                        sid2_s = id2gen_se_dic[sid2][0]
                        sid2_e = id2gen_se_dic[sid2][1]
                        gen_dist = get_site_ends_distance(sid_s, sid_e, sid2_s, sid2_e)
                        if gen_dist <= max_site_dist:
                            rem_assoc_sids_dic[sid2] = 1
                            break
            if rem_assoc_sids_dic:
                new_list = []
                for sid in id2ids_dic[site_id]:
                    if sid not in rem_assoc_sids_dic:
                        new_list.append(sid)
                if new_list:
                    new_list.sort()
                    id2ids_dic[site_id] = new_list
                else:
                    rem_sids_dic[site_id] = 1
    if rem_sids_dic:
        for rem_sid in rem_sids_dic:
            del id2ids_dic[rem_sid]

    return id2ids_dic, ids2isrc_dic


################################################################################

def remove_no_common_tr_exb_pairs(id2ids_dic, id2tr_list_dic):
    """
    Remove exon border pairs which have no common transcripts.

    >>> id2ids_dic = {'s1': ['s2', 's3'], 's2': ['s1'], 's3': ['s1']}
    >>> id2tr_list_dic = {'s1': ['t1', 't2'], 's2': ['t3'], 's3': ['t2']}
    >>> remove_no_common_tr_exb_pairs(id2ids_dic, id2tr_list_dic)
    {'s1': ['s3'], 's3': ['s1']}

    """
    if not id2ids_dic:
        return id2ids_dic
    assert id2tr_list_dic, "id2tr_list_dic empty"
    rem_sid_list = []
    for sid1 in id2ids_dic:
        if sid1 not in id2tr_list_dic:
            rem_sid_list.append(sid1)
            continue
        sid1_tr_list = id2tr_list_dic[sid1]
        new_list = []
        for sid2 in id2ids_dic[sid1]:
            sid2_tr_list = id2tr_list_dic[sid2]
            if two_lists_get_intersect(sid1_tr_list, sid2_tr_list):
                new_list.append(sid2)
        # If no connections left, remove site_id from id2ids_dic.
        if not new_list:
            rem_sid_list.append(sid1)
        else:
            new_list.sort()
            id2ids_dic[sid1] = new_list

    if rem_sid_list:
        for rem_sid in rem_sid_list:
            del id2ids_dic[rem_sid]

    return id2ids_dic


################################################################################

def get_connected_exb_sites(site_id, id2exb_pair_dic, seen_dic):
    """
    Recursively get connected exon border pair sites.

    >>> id2exb_pair_dic = {'id1': 'id2', 'id2': 'id3', 'id3': 'id4', 'id4': 'id3', 'id5': 'id6', 'id6': 'id5'}
    >>> seen_dic = {}
    >>> get_connected_exb_sites("id1", id2exb_pair_dic, seen_dic)
    >>> seen_dic
    {'id1': 1, 'id2': 1, 'id3': 1, 'id4': 1}
    >>> seen_dic = {}
    >>> get_connected_exb_sites("id5", id2exb_pair_dic, seen_dic)
    >>> seen_dic
    {'id5': 1, 'id6': 1}

    """
    seen_dic[site_id] = 1
    assert site_id in id2exb_pair_dic, "site ID %s not in id2exb_pair_dic" %(site_id)
    pair_id = id2exb_pair_dic[site_id]
    if pair_id in seen_dic:
        return 0
    else:
        get_connected_exb_sites(pair_id, id2exb_pair_dic, seen_dic)


################################################################################

def get_highest_isr_exb_pair(con_exb_sites_dic, ids2isrc_dic):
    """
    Given a dictionary of connected exon border site IDs, get pair with
    highest intron-spanning read count between them.

    If all have same ISR count, return the first one seen.

    >>> con_exb_sites_dic = {'id1': 1, 'id2': 1, 'id3': 1, 'id4': 1}
    >>> ids2isrc_dic = {"id1-EB-id2": 10, "id2-EB-id1": 10, "id2-EB-id3": 15, "id3-EB-id2": 15, "id3-EB-id4": 20, "id4-EB-id3": 20}
    >>> get_highest_isr_exb_pair(con_exb_sites_dic, ids2isrc_dic)
    (['id3', 'id4'], 20)

    """

    eb_ids_list = []
    seen_dic = {}
    best_pair = []
    best_isrc = 0
    sids_list = []
    for sid1 in con_exb_sites_dic:
        sids_list.append(sid1)
        for sid2 in con_exb_sites_dic:
            if sid1 == sid2:
                continue
            sids1 = sid1 + "-EB-" + sid2
            sids2 = sid2 + "-EB-" + sid1
            if sids1 in seen_dic:
                continue
            seen_dic[sids1] = 1
            seen_dic[sids2] = 1
            if sids1 in ids2isrc_dic:
                isrc = ids2isrc_dic[sids1]
                if isrc > best_isrc:
                    best_isrc = isrc
                    best_pair = [sid1, sid2]

    assert best_pair, "no highest ISRC pair extracted for connected exon border sites %s" %(",".join(sids_list))
    best_pair.sort()
    return best_pair, best_isrc


################################################################################

def get_exb_group_best_con(exb_group_ids_list, id2exb_pair_dic):
    """
    Given a list of site IDs connected at exon borders with intron-spanning
    reads. Only one of them should have a connection in both directions,
    while the other have different connections.

    >>> exb_group_ids_list = ['id1', 'id2', 'id3']
    >>> id2exb_pair_dic = {'id1': 'id2', 'id2': 'id3', 'id3': 'id2'}
    >>> get_exb_group_best_con(exb_group_ids_list, id2exb_pair_dic)
    ['id2', 'id3']

    """
    assert exb_group_ids_list, "exb_group_ids_list empty"
    sids = []
    for sid in exb_group_ids_list:
        pid = id2exb_pair_dic[sid]
        sid2 = id2exb_pair_dic[pid]
        if sid == sid2:
            sids = [sid, pid]
            sids.sort()
            break
    assert sids, "no exon border site pair with connections in both directions found for %s" %(','.join(exb_group_ids_list))
    return sids


################################################################################

def get_border_pair_exons(site_id, id2exids_dic, exid2trid_dic,
                          id2ids_dic, ids2isrc_dic,
                          min_isr_c=2):
    """
    Given a site ID and overlapping exon IDs (id2exids_dic) > min_eir,
    check if site ID is at exon border and is connected through this
    border via intron-spanning reads to another border site.
    This info is given by: id2ids_dic, ids2isrc_dic
    Return list of exon IDs which are compatible with this connection,
    i.e. exon IDs from transcripts both sites possibly share. If more
    than one connection is present, choose the connection with most
    connecting reads.

    id2exids_dic:
        site ID -> exon IDs > min_eir mapping.
    exid2trid_dic:
        exon ID -> transcript ID mapping.
    id2ids_dic:
        site ID -> site IDs connected by exon border list mapping.
    ids2isrc_dic:
        site ID pair -> connecting intron-spanning read count.

    Several (special) cases possible:
    1) > 1 connection between the respective exon border site and other sites.
    In this case choose the connection with highest ISR count.
    2) Support by ISR counts but no common transcript ID. In this case
    keep the existing exon IDs.
    3) Only one exon ID. In this case keep the ID, even if it does not support
    the connection.
    4) No connection. Here return all exon IDs / do not change the list.

    >>> id2exids_dic = {"id1": ["t1_e1", "t2_e2", "t4_e2"], "id2" : ["t1_e2", "t2_e3", "t3_e3"], "id3":  ["t3_e4", "t7_e1"], "id4": ["t4_e3", "t5_e1"], "id6": ["t6_e666"]}
    >>> exid2trid_dic = {"t1_e1": "t1", "t1_e2": "t1", "t2_e2": "t2", "t2_e3": "t2", "t3_e3": "t3", "t3_e4": "t3", "t4_e2": "t4", "t4_e3": "t4", "t6_e666": "id6", "t5_e1" : "t5", "t7_e1" : "t7"}
    >>> id2ids_dic = {"id1": ["id2", "id4"], "id2": ["id1"], "id4": ["id1"], "id6": ["id3"], "id3": ["id6"]}
    >>> ids2isrc_dic = {"id1-EB-id2": 10, "id2-EB-id1": 10, "id1-EB-id4": 2, "id4-EB-id1": 2, "id3-EB-id6": 5, "id6-EB-id3": 5}
    >>> get_border_pair_exons("id1", id2exids_dic, exid2trid_dic, id2ids_dic, ids2isrc_dic)
    (['t1_e1', 't2_e2'], 'id2')
    >>> get_border_pair_exons("id2", id2exids_dic, exid2trid_dic, id2ids_dic, ids2isrc_dic)
    (['t1_e2', 't2_e3'], 'id1')
    >>> get_border_pair_exons("id3", id2exids_dic, exid2trid_dic, id2ids_dic, ids2isrc_dic)
    (['t3_e4', 't7_e1'], '')
    >>> get_border_pair_exons("id4", id2exids_dic, exid2trid_dic, id2ids_dic, ids2isrc_dic)
    (['t4_e3'], 'id1')
    >>> get_border_pair_exons("id4", id2exids_dic, exid2trid_dic, id2ids_dic, ids2isrc_dic, min_isr_c=3)
    (['t4_e3', 't5_e1'], '')
    >>> get_border_pair_exons("id6", id2exids_dic, exid2trid_dic, id2ids_dic, ids2isrc_dic)
    (['t6_e666'], '')

    """
    if site_id not in id2ids_dic:
        return id2exids_dic[site_id], ""
    # if len(id2exids_dic[site_id]) == 1:
    #     return id2exids_dic[site_id], ""
    else:
        # Select best connection.
        best_con_c = 0
        best_con_id = ""
        for sid in id2ids_dic[site_id]:
            ids = "%s-EB-%s" %(site_id, sid)
            con_c = ids2isrc_dic[ids]
            if con_c > best_con_c:
                best_con_c = con_c
                best_con_id = sid
        assert best_con_c, "no best connection ID selected"
        # Look for common transcripts.
        ex_ids1 = id2exids_dic[site_id]
        ex_ids2 = id2exids_dic[best_con_id]
        if best_con_c < min_isr_c:
            return id2exids_dic[site_id], ""
        common_ex_ids = []
        for exid1 in ex_ids1:
            m = re.search(".+_e(\d+)", exid1)
            assert m, "invalid exon ID %s" %(exid1)
            exid1_nr = int(m.group(1))
            trid1 = exid2trid_dic[exid1]
            for exid2 in ex_ids2:
                m = re.search(".+_e(\d+)", exid2)
                assert m, "invalid exon ID %s" %(exid2)
                exid2_nr = int(m.group(1))
                trid2 = exid2trid_dic[exid2]
                """
                If both sites support same transcript and their exons are
                adjacent to each other!
                """
                if trid1 == trid2 and abs(exid1_nr-exid2_nr) == 1:
                    common_ex_ids.append(exid1)
        # If no common transcript.
        if common_ex_ids:
            return common_ex_ids, best_con_id
        else:
            return id2exids_dic[site_id], ""


################################################################################

def select_id_from_scores_dic(id1, id2, sc_dic,
                              get_worse=False,
                              rev_filter=False):
    """
    Based on ID to score mapping, return better (or worse) scoring ID.

    >>> id1 = "id1"
    >>> id2 = "id2"
    >>> id3 = "id3"
    >>> sc_dic = {'id1' : 5, 'id2': 3, 'id3': 3}
    >>> select_id_from_scores_dic(id1, id2, sc_dic)
    'id1'
    >>> select_id_from_scores_dic(id1, id2, sc_dic, get_worse=True)
    'id2'
    >>> select_id_from_scores_dic(id1, id2, sc_dic, rev_filter=True, get_worse=True)
    'id1'
    >>> select_id_from_scores_dic(id1, id2, sc_dic, rev_filter=True)
    'id2'
    >>> select_id_from_scores_dic(id2, id3, sc_dic)
    False

    """

    sc_id1 = sc_dic[id1]
    sc_id2 = sc_dic[id2]
    if sc_id1 > sc_id2:
        if rev_filter:
            if get_worse:
                return id1
            else:
                return id2
        else:
            if get_worse:
                return id2
            else:
                return id1
    elif sc_id1 < sc_id2:
        if rev_filter:
            if get_worse:
                return id2
            else:
                return id1
        else:
            if get_worse:
                return id1
            else:
                return id2
    else:
        return False


################################################################################

def bam_check_file_empty(in_bam):
    """
    Check for empty BAM file.

    >>> empty_bam = "test_data/empty.bam"
    >>> bam_check_file_empty(empty_bam)
    True

    """
    assert os.path.exists(in_bam), "in_bam does not exist"


    check_cmd = "samtools view -c " + in_bam
    output = int(subprocess.getoutput(check_cmd).strip())
    if output:
        return False
    else:
        return True


################################################################################

def bam_count_reads(in_bam):
    """
    Count and return reads inside in_bam.

    >>> empty_bam = "test_data/empty.bam"
    >>> bam_count_reads(empty_bam)
    0

    """
    assert os.path.exists(in_bam), "in_bam does not exist"

    check_cmd = "samtools view -c " + in_bam
    count = int(subprocess.getoutput(check_cmd).strip())
    return count


################################################################################

def filter_bam_file(in_bam, out_bam,
                    pp_mode=2):
    """
    Filter input bam file in_bam to keep only R2 reads (pp_mode=2),
    or keep only R1 reads (pp_mode=3).

    """
    # Filter.
    filter_flag = "130"
    if pp_mode == 3:
        filter_flag = "0x40"
    check_cmd = "samtools view -hb -f " + filter_flag + " " + in_bam + " -o " + out_bam
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "samtools view has problems with your input:\n%s\n%s" %(check_cmd, output)


################################################################################

def get_extended_gen_seqs(args, id2row_dic,
                          ref_len_dic=False,
                          id2out_dic=False,
                          tmp_out_folder=False,
                          rr_ratios_dic=None):
    """
    Extend genomic regions and return extended genomic site sequences.

    """
    random_id = uuid.uuid1()
    zero_sc_tmp_bed = str(random_id) + ".zero_sc.tmp.bed"
    random_id = uuid.uuid1()
    zero_sc_tmp_fa = str(random_id) + ".zero_sc.tmp.fa"
    if tmp_out_folder:
        zero_sc_tmp_bed = tmp_out_folder + "/" + zero_sc_tmp_bed
        zero_sc_tmp_fa = tmp_out_folder + "/" + zero_sc_tmp_fa

    # Output BED regions with zero scores.
    bed_write_row_dic_into_file(id2row_dic, zero_sc_tmp_bed,
                                ext_mode=args.seq_ext_mode,
                                ext_lr=args.seq_ext,
                                zero_scores=True,
                                id2out_dic=id2out_dic,
                                chr_len_dic=ref_len_dic)

    # Get genomic sequences.
    bed_extract_sequences_from_2bit(zero_sc_tmp_bed,
                                    zero_sc_tmp_fa,
                                    args.in_2bit,
                                    lc_repeats=True)
    site_seqs_dic = read_fasta_into_dic(zero_sc_tmp_fa,
                                        dna=False,
                                        skip_n_seqs=False)

    # Calculate repeat region ratios for each site.
    if rr_ratios_dic is not None:
        gen_rr_ratios_dic = get_seqs_dic_repeat_region_ratios(site_seqs_dic)
        for site_id in gen_rr_ratios_dic:
            rr_ratios_dic[site_id] = gen_rr_ratios_dic[site_id]

    if os.path.exists(zero_sc_tmp_bed):
        os.remove(zero_sc_tmp_bed)
    if os.path.exists(zero_sc_tmp_fa):
        os.remove(zero_sc_tmp_fa)

    return site_seqs_dic


################################################################################

def pm_ext_merge_bed_regions(id2row_dic,
                             ref_len_dic,
                             id2sc_dic=None,
                             id2len_dic=None,
                             id2gen_se_dic=None,
                             new_stem_id=False,
                             tmp_out_folder=False,
                             merge_ext=0):
    """
    peakhood extract --pre-merge
    Extend and merge BED regions, return merged regions (not best like
    in other merge operations).
    Return new ID to row dictionary.

    """
    random_id = uuid.uuid1()
    m1_tmp_bed = str(random_id) + ".pre_merge1.tmp.bed"
    random_id = uuid.uuid1()
    m2_tmp_bed = str(random_id) + ".pre_merge2.tmp.bed"
    if tmp_out_folder:
        m1_tmp_bed = tmp_out_folder + "/" + m1_tmp_bed
        m2_tmp_bed = tmp_out_folder + "/" + m2_tmp_bed

    bed_write_row_dic_into_file(id2row_dic, m1_tmp_bed,
                                ext_mode=1,
                                ext_lr=merge_ext,
                                zero_scores=False,
                                chr_len_dic=ref_len_dic)

    bed_sort_merge_output_ol_regions(m1_tmp_bed, m2_tmp_bed,
                                     tmp_out_folder=tmp_out_folder,
                                     new_stem_id=new_stem_id)

    pm_id2row_dic = bed_read_rows_into_dic(m2_tmp_bed,
                                           id2sc_dic=id2sc_dic,
                                           id2len_dic=id2len_dic,
                                           check_chr_id_format=False,
                                           id2gen_se_dic=id2gen_se_dic)

    if os.path.exists(m1_tmp_bed):
        os.remove(m1_tmp_bed)
    if os.path.exists(m2_tmp_bed):
        os.remove(m2_tmp_bed)

    return pm_id2row_dic


################################################################################

def merge_filter_bam_files(list_bam, out_bam,
                           tmp_out_folder=False,
                           pp_mode=1):
    """
    Preprocess and filter input --bam files.

    pp_mode:
        BAM preprocessing mode.
        1: no filtering after merging.
        2: filter to keep only R2 reads.
        3: filter to keep only R1 reads.

    """
    bam_files = ""
    for bam_file in list_bam:
        bam_files += " %s" %(bam_file)

    if pp_mode != 1:
        random_id = uuid.uuid1()
        tmp_bam = str(random_id) + ".tmp.bam"
        if tmp_out_folder:
            tmp_bam = tmp_out_folder + "/" + tmp_bam

        # Merge.
        check_cmd = "samtools merge -f " + tmp_bam + " " + bam_files
        output = subprocess.getoutput(check_cmd)
        error = False
        if output:
            error = True
        assert error == False, "samtools merge has problems with your input:\n%s\n%s" %(check_cmd, output)

        # Filter.
        filter_flag = "130"
        if pp_mode == 3:
            filter_flag = "0x40"
        # -hb include header and output BAM.
        check_cmd = "samtools view -hb -f " + filter_flag + " " + tmp_bam + " -o " + out_bam
        output = subprocess.getoutput(check_cmd)
        error = False
        if output:
            error = True
        assert error == False, "samtools view has problems with your input:\n%s\n%s" %(check_cmd, output)

        # Delete tmp files.
        if os.path.exists(tmp_bam):
            os.remove(tmp_bam)
    else:
        # Merge.
        check_cmd = "samtools merge -f " + out_bam + " " + bam_files
        output = subprocess.getoutput(check_cmd)
        error = False
        if output:
            error = True
        assert error == False, "samtools merge has problems with your input:\n%s\n%s" %(check_cmd, output)


################################################################################

def intersect_bed_files(a_file, b_file, params, out_file,
                        a2b_list_dic=None,
                        ab2ol_se_dic=None,
                        sorted_out=False):
    """
    Intersect two .bed files, using intersectBed.

    a2b_list_dic:
        If -wb is chosen, this dictionary is used to store -a ID to -b IDs
        stored as list.
    ab2ol_se_dic:
        a_id,b_id -> [a_overlap_s, a_overlap_e] (zero-based, one-based)

    """

    check_cmd = "intersectBed -a " + a_file + " -b " + b_file + " " + params + " > " + out_file
    if sorted_out:
        check_cmd = "intersectBed -a " + a_file + " -b " + b_file + " " + params + " | " + "sort -k1,1 -k2,2n > " + out_file
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "intersectBed has problems with your input:\n%s\n%s" %(check_cmd, output)

    if a2b_list_dic is not None:
        f = open(out_file, "r")
        for line in f:
            cols = line.strip().split("\t")
            a_s = int(cols[1])
            a_e = int(cols[2])
            a_id = cols[3]
            b_id = cols[9]
            if ab2ol_se_dic is not None:
                ab_id = "%s,%s" %(a_id, b_id)
                ab2ol_se_dic[ab_id] = [a_s, a_e]
            if a_id in a2b_list_dic:
                a2b_list_dic[a_id].append(b_id)
            else:
                a2b_list_dic[a_id] = [b_id]
        f.close()
        # Make list entries unique.
        for a_id in a2b_list_dic:
            a2b_list_dic[a_id] = list(set(a2b_list_dic[a_id]))


################################################################################

def intersect_genes_with_sites(genes_bed, sites_bed, tmp_out):
    """
    Intersect gene regions with sites, and return site ID to genes mapping,
    and gene regions dictionary to output later as BED.

    """

    params = "-s -wa -wb"
    check_cmd = "intersectBed -a " + genes_bed + " -b " + sites_bed + " " + params + " > " + tmp_out
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "intersectBed has problems with your input:\n%s\n%s" %(check_cmd, output)

    id2gene_list_dic = {}
    gene2reg_dic = {}
    f = open(tmp_out, "r")
    for line in f:
        cols = line.strip().split("\t")
        chr_id = cols[0]
        gene_s = int(cols[1])
        gene_e = int(cols[2])
        gene_id = cols[3]
        gene_pol = cols[5]
        site_id = cols[9]
        if site_id in id2gene_list_dic:
            id2gene_list_dic[site_id].append(gene_id)
        else:
            id2gene_list_dic[site_id] = [gene_id]
        gene2reg_dic[gene_id] = [chr_id, gene_s, gene_e, gene_pol]
    f.close()
    return id2gene_list_dic, gene2reg_dic


################################################################################

def bed_write_reg_list_to_file(id2reg_dic, out_bed,
                               id2out_dic=None):
    """
    Write dictionary of region lists to BED file.
    Example dictionary:
    {'gene1': ["chr1", 10, 20, "+"], ...}

    id2out_dic:
        IDs dictionary for which to output regions.

    >>> id2reg_dic = {"site_666" : ["chr6", 66, 666, "-"], "camping_site" : ["chr10", 10, 100, "+"]}
    >>> out_exp_bed = "test_data/reg_list_out.exp.bed"
    >>> out_tmp_bed = "test_data/reg_list_out.tmp.bed"
    >>> bed_write_reg_list_to_file(id2reg_dic, out_tmp_bed)
    >>> diff_two_files_identical(out_exp_bed, out_tmp_bed)
    True

    """
    assert id2reg_dic, "given id2reg_dic empty"
    OUTBED = open(out_bed, "w")
    c_out = 0
    for reg_id in id2reg_dic:
        reg_list = id2reg_dic[reg_id]
        if id2out_dic is not None:
            if not reg_id in id2out_dic:
                continue
        c_out += 1
        OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" %(reg_list[0], reg_list[1], reg_list[2], reg_id, reg_list[3]))
    OUTBED.close()
    assert c_out, "nothing was output"


################################################################################

def get_site_len_list(len_in):
    """
    Read in site lengths (one length per row from file len_in).
    Return list of lengths.

    """
    site_len_list = []

    with open(len_in) as f:
        for line in f:
            site_len = line.strip()
            site_len_list.append(int(site_len))
    f.closed
    assert site_len_list, "site_len_list empty (no lengths read in from %s)" %(len_in)
    return site_len_list


################################################################################

def gtf_check_gene_feat(in_gtf,
                        n_rows_check=10000):
    """
    Extract gene regions from in_gtf GTF file, and output to out_bed BED
    file.

    n_rows_check:
        Number of rows to check.

    >>> true_gtf = "test_data/gene_test_in.gtf"
    >>> false_gtf = "test_data/test_order_true.gtf"
    >>> gtf_check_gene_feat(true_gtf)
    True
    >>> gtf_check_gene_feat(false_gtf)
    False

    """

    c_in = 0
    check = False
    if re.search(".+\.gz$", in_gtf):
        f = gzip.open(in_gtf, 'rt')
    else:
        f = open(in_gtf, "r")
    for line in f:
        # Skip header.
        if re.search("^#", line):
            continue
        cols = line.strip().split("\t")
        feature = cols[2]
        c_in += 1
        if c_in > n_rows_check:
            break
        if feature == "gene":
            check = True
            break
    f.close()
    return check


################################################################################

def bed_intersect_sites_genes_get_infos(sites_bed, genes_bed, id2gids_dic,
                                        tmp_out_folder=False):
    """
    Intersect gene regions with sites, and return site_id -> overlapping
    gene IDs mapping.

    >>> sites_bed = "test_data/test_intersect.sites.bed"
    >>> genes_bed = "test_data/test_intersect.genes.bed"
    >>> id2gids_dic = {}
    >>> bed_intersect_sites_genes_get_infos(sites_bed, genes_bed, id2gids_dic)
    >>> id2gids_dic
    {'site1': ['ENSG1', 'ENSG2'], 'site2': ['ENSG1'], 'site4': ['ENSG3']}

    """

    params = "-wb -s -f 0.8"

    # Generate .tmp files.
    random_id = uuid.uuid1()
    tmp_out = str(random_id) + ".intersect.tmp.out"
    if tmp_out_folder:
        tmp_out = tmp_out_folder + "/" + tmp_out

    check_cmd = "intersectBed -a " + sites_bed + " -b " + genes_bed + " " + params
    output = subprocess.getoutput(check_cmd)

    for line in output.split('\n'):
        cols = line.strip().split("\t")
        site_id = cols[3]
        gene_id = cols[9]
        if site_id in id2gids_dic:
            id2gids_dic[site_id].append(gene_id)
        else:
            id2gids_dic[site_id] = [gene_id]


################################################################################

def gtf_extract_gene_bed(in_gtf, out_bed,
                         gene_ids_dic=False,
                         gid2gn_dic=None,
                         chr_id_style=1,
                         chr_ids_dic=None,
                         gid2gbt_dic=None):
    """
    Extract gene regions from in_gtf GTF file, and output to out_bed BED
    file.

    gene_ids_dic:
    Dictionary with gene IDs for filtering (keeping dic IDs).

    >>> in_gtf = "test_data/gene_test_in.gtf"
    >>> exp_out_bed = "test_data/gtf_gene_out.exp.bed"
    >>> tmp_out_bed = "test_data/gtf_gene_out.tmp.bed"
    >>> gtf_extract_gene_bed(in_gtf, tmp_out_bed)
    >>> diff_two_files_identical(tmp_out_bed, exp_out_bed)
    True

    """

    # Output gene regions.
    OUTBED = open(out_bed, "w")
    c_out = 0
    # Open GTF either as .gz or as text file.
    if re.search(".+\.gz$", in_gtf):
        f = gzip.open(in_gtf, 'rt')
    else:
        f = open(in_gtf, "r")
    for line in f:
        # Skip header.
        if re.search("^#", line):
            continue
        cols = line.strip().split("\t")
        chr_id = cols[0]
        feature = cols[2]
        feat_s = int(cols[3])
        feat_e = int(cols[4])
        feat_pol = cols[6]
        infos = cols[8]
        if not feature == "gene":
            continue

        # Restrict to standard chromosomes.
        new_chr_id = check_convert_chr_id(chr_id,
                                id_style=chr_id_style)
        if not new_chr_id:
            continue
        else:
            chr_id = new_chr_id

        # Make start coordinate 0-base (BED standard).
        feat_s = feat_s - 1

        # Extract gene ID and from infos.
        m = re.search('gene_id "(.+?)"', infos)
        assert m, "gene_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        gene_id = m.group(1)

        # Check if gene ID is in gene dic.
        if gene_ids_dic:
            if not gene_id in gene_ids_dic:
                continue

        # Only for GTF files with gene feature present (GENCODE, Ensembl .. ).
        gene_name = "-"
        gene_biotype = "-"
        m = re.search('gene_name "(.+?)"', infos)
        if m:
            gene_name = m.group(1)

        m = re.search('gene_biotype "(.+?)"', infos)
        # Try Gencode encoding.
        if not m:
            m = re.search('gene_type "(.+?)"', infos)
        if m:
            gene_biotype = m.group(1)

        if gid2gn_dic is not None:
            gid2gn_dic[gene_id] = gene_name
        if gid2gbt_dic is not None:
            gid2gbt_dic[gene_id] = gene_biotype

        if chr_ids_dic is not None:
            chr_ids_dic[chr_id] = 1

        # Output gene region.
        c_out += 1
        OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,feat_s,feat_e,gene_id,feat_pol))

    OUTBED.close()
    f.close()

    assert c_out, "no regions output to out_bed. Invalid in_gtf (e.g. chromosome IDs inside --gtf should be either of type \"1\", or \"chr1\"), or too restrictive gene_ids_dic filtering?"


################################################################################

def read_ids_into_dic(ids_file,
                      check_dic=True,
                      ids_dic=False):
    """
    Read in IDs file, where each line stores one ID.

    >>> test_ids_file = "test_data/test.ids"
    >>> ids_dic = read_ids_into_dic(test_ids_file)
    >>> print(ids_dic)
    {'clip1': 1, 'clip2': 1, 'clip3': 1}

    """
    if not ids_dic:
        ids_dic = {}
    # Read in file content.
    with open(ids_file) as f:
        for line in f:
            row_id = line.strip()
            ids_dic[row_id] = 1
    f.closed
    if check_dic:
        assert ids_dic, "IDs dictionary ids_dic empty"
    return ids_dic


################################################################################

def gtf_get_transcript_ids(in_gtf,
                           chr_ids_dic=None):
    """
    Get transcript IDs from in_gtf GTF file.

    >>> in_gtf = "test_data/gene_test_in.gtf"
    >>> gtf_get_transcript_ids(in_gtf)
    {'ENST01': 1, 'ENST02': 1}

    """
    # Transcript IDs dictionary.
    tr_ids_dic = {}
    # Open GTF either as .gz or as text file.
    if re.search(".+\.gz$", in_gtf):
        f = gzip.open(in_gtf, 'rt')
    else:
        f = open(in_gtf, "r")
    for line in f:
        # Skip header.
        if re.search("^#", line):
            continue
        cols = line.strip().split("\t")
        chr_id = cols[0]
        feature = cols[2]
        infos = cols[8]
        if not feature == "transcript":
            continue

        if chr_ids_dic is not None:
            chr_ids_dic[chr_id] = 1

        # Extract transcript ID.
        m = re.search('transcript_id "(.+?)"', infos)
        assert m, "transcript_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        transcript_id = m.group(1)

        # Store transcript ID.
        tr_ids_dic[transcript_id] = 1
    f.close()
    # Check and return to barracks.
    assert tr_ids_dic, "no transcript IDs read in"
    return tr_ids_dic


################################################################################

def bed_get_region_pols(in_bed):
    """
    Read in .bed file, and store polarities for each region in dictionary
    (unique column 4 ID has to be present).
    Return dictionary with mappings region ID -> region polarity

    >>> test_bed = "test_data/test4.bed"
    >>> bed_get_region_pols(test_bed)
    {'tr1_e1': '+', 'tr2_e1': '-'}

    """
    id2pol_dic = {}
    with open(in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            site_id = cols[3]
            site_pol = cols[5]
            id2pol_dic[site_id] = site_pol
    f.closed
    assert id2pol_dic, "nothing read in for in_bed \"%s\"" %(in_bed)
    return id2pol_dic


################################################################################

def bed_merge_transcript_regions(in_bed, out_bed,
                                 new_ids=False,
                                 id2pol_dic=False):
    """
    Merge transcript regions BED (strand specific).

    id2pol_dic:
        Region ID to genomic polarity mapping.

    >>> in_bed = "test_data/test.merge_tr_reg_in.bed"
    >>> out_bed = "test_data/test.merge_tr_reg_out.tmp.bed"
    >>> exp_bed = "test_data/test.merge_tr_reg_out.exp.bed"
    >>> bed_merge_transcript_regions(in_bed, out_bed)
    >>> diff_two_files_identical(out_bed, exp_bed)
    True

    """
    assert os.path.isfile(in_bed), "cannot open in_bed \"%s\"" % (in_bed)

    # Get region polarities.
    if not id2pol_dic:
        id2pol_dic = bed_get_region_pols(in_bed)

    # Sort and merge transcript regions.
    check_cmd = 'sort -k1,1 -k2,2n ' + in_bed + ' | mergeBed -i stdin -s -c 4 -o distinct -delim ";"'
    output = subprocess.getoutput(check_cmd)

    MRGOUT = open(out_bed, "w")
    c_read = 0
    for line in output.split('\n'):
        cols = line.strip().split("\t")
        c_read += 1
        chr_id = cols[0]
        reg_s = cols[1]
        reg_e = cols[2]
        merged_ids_list = cols[3].split(";")
        reg_pol = id2pol_dic[merged_ids_list[0]]
        new_id = cols[3]
        if new_ids:
            new_id = "tr_reg_%i" %(c_read)
        MRGOUT.write("%s\t%s\t%s\t%s\t0\t%s\n" %(chr_id, reg_s, reg_e, new_id, reg_pol))
    MRGOUT.close()
    assert c_read, "no merged transcript regions output"


################################################################################

def bed_sort_merge_output_top_entries(in_bed, out_bed,
                                      alpha_merge=False,
                                      check_chr_id_format=False,
                                      rev_filter=False):
    """
    Sort in_bed file, use mergeBed from bedtools to merge overlapping entries,
    then select for each overlapping set the entry with highest score and
    output it to out_bed.

    alpha_merge:
        Alphabetical merge (if no site scores available).

    >>> in_bed = "test_data/test5.bed"
    >>> out_bed = "test_data/test5.tmp.bed"
    >>> exp_bed = "test_data/test5.exp.bed"
    >>> bed_sort_merge_output_top_entries(in_bed, out_bed)
    >>> diff_two_files_identical(out_bed, exp_bed)
    True

    """
    assert os.path.isfile(in_bed), "cannot open in_bed \"%s\"" % (in_bed)
    # Generate .tmp files.
    random_id = uuid.uuid1()
    tmp_bed = str(random_id) + ".tmp.bed"
    # Read in_bed rows into dictionary.
    id2row_dic = bed_read_rows_into_dic(in_bed,
                                        check_chr_id_format=check_chr_id_format)
    # Get region scores.
    id2sc_dic = bed_get_region_id_scores(in_bed)

    # Sort file.
    bed_sort_file(in_bed, out_bed)
    # Merge .bed.
    bed_merge_file(out_bed, tmp_bed)
    # Output file.
    OUTBED = open(out_bed,"w")
    # Open merged .bed file, and select top entry for each overlap set.
    with open(tmp_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            ids = cols[3].split(";")
            if alpha_merge:
                ids.sort()
                best_id = ids[0]
            else:
                best_id = "-"
                best_sc = -666666
                if rev_filter:
                    best_sc = 666666
                for site_id in ids:
                    assert site_id in id2sc_dic, "site ID \"%s\" not found in id2sc_dic" % (site_id)
                    site_sc = id2sc_dic[site_id]
                    if rev_filter:
                        if site_sc < best_sc:
                            best_sc = site_sc
                            best_id = site_id
                    else:
                        if site_sc > best_sc:
                            best_sc = site_sc
                            best_id = site_id
            assert best_id in id2row_dic, "site ID \"%s\" not found in id2row_dic" % (best_id)
            OUTBED.write(id2row_dic[best_id] + "\n")
    f.closed
    OUTBED.close()
    if os.path.exists(tmp_bed):
        os.remove(tmp_bed)


################################################################################

def bed_sort_merge_output_ol_regions(in_bed, out_bed,
                                     tmp_out_folder=False,
                                     new_stem_id=False):
    """
    Sort in_bed file, use mergeBed from bedtools to merge overlapping entries,
    then select for each overlapping set the entry with highest score and
    output it to out_bed.

    new_stem_id:
        New stem ID for assigning new IDs to merged regions.
        By default, merge site IDs in regions to new IDs.

    test5.bed
    chr1	3000	4000	CLIP2	2.57	+
    chr2	1000	2500	CLIP1	1.58	+
    chr1	3500	5000	CLIP3	3.11	+

    test5_2.exp.bed
    chr1	3000	5000	CLIP2_CLIP3	2.84	+
    chr2	1000	2500	CLIP1	1.58	+

    >>> in_bed = "test_data/test5.bed"
    >>> out_bed = "test_data/test5_2.tmp.bed"
    >>> exp_bed = "test_data/test5_2.exp.bed"
    >>> bed_sort_merge_output_ol_regions(in_bed, out_bed)
    >>> diff_two_files_identical(out_bed, exp_bed)
    True

    """
    assert os.path.isfile(in_bed), "cannot open in_bed \"%s\"" % (in_bed)
    # Generate .tmp files.
    random_id = uuid.uuid1()
    tmp_bed = str(random_id) + ".tmp.bed"
    if tmp_out_folder:
        tmp_bed = tmp_out_folder + "/" + tmp_bed

    # Get region scores.
    id2pol_dic = {}
    id2sc_dic = bed_get_region_id_scores(in_bed,
                                         id2pol_dic=id2pol_dic)

    # Sort file.
    bed_sort_file(in_bed, out_bed)
    # Merge .bed.
    bed_merge_file(out_bed, tmp_bed)

    """
    chr1	1000	2000	r1	1	+
    chr1	2000	3000	r2	2	+
    chr1	2500	2700	r4	4	+
    chr1	5000	6000	r7	3	+

    merged:
    chr1	1000	3000	r1;r2;r4
    chr1	5000	6000	r7
    """

    OUTBED = open(out_bed,"w")
    c_out = 0

    with open(tmp_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            chr_id = cols[0]
            gen_s = cols[1]
            gen_e = cols[2]
            ids = cols[3].split(";")
            c_out += 1
            # New score == average score.
            sc_list = []
            for sid in ids:
                sc_list.append(id2sc_dic[sid])
            if len(sc_list) > 1:
                new_sc = statistics.mean(sc_list)
            else:
                new_sc = sc_list[0]
            # New site ID.
            if new_stem_id:
                new_id = new_stem_id + "_" + str(c_out)
            else:
                new_id = '_'.join(ids)
            # Region polarity.
            site_pol = id2pol_dic[ids[0]]

            OUTBED.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(chr_id, gen_s, gen_e, new_id, new_sc, site_pol))

    f.closed
    OUTBED.close()
    if os.path.exists(tmp_bed):
        os.remove(tmp_bed)


################################################################################

def bed_get_region_id_scores(in_bed, no_float=False,
                             id2pol_dic=None):
    """
    Read in .bed file, and store scores for each region in dictionary
    (unique column 4 ID and column 5 score have to be present).
    Return dictionary with mappings region ID -> region score

    >>> test_bed = "test_data/test5.bed"
    >>> bed_get_region_id_scores(test_bed)
    {'CLIP2': 2.57, 'CLIP1': 1.58, 'CLIP3': 3.11}

    """
    id2sc_dic = {}
    # Open input .bed file.
    with open(in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            site_id = cols[3]
            site_sc = float(cols[4])
            site_pol = cols[5]
            if no_float:
                site_sc = cols[4]
            id2sc_dic[site_id] = site_sc
            if id2pol_dic is not None:
                id2pol_dic[site_id] = site_pol
    f.closed
    assert id2sc_dic, "nothing read in for in_bed \"%s\"" %(in_bed)
    return id2sc_dic


################################################################################

def bed_merge_file(in_bed, out_bed,
                   custom_params_str=False):
    """
    Use mergeBed from bedtools to merge overlapping .bed entries, storing
    the region IDs to later pick one region for each set of overlapping
    regions.

    >>> in_bed = "test_data/test.sorted.bed"
    >>> out_bed = "test_data/test.sorted.merged.tmp.bed"
    >>> out_exp_bed = "test_data/test.sorted.merged.exp.bed"
    >>> bed_merge_file(in_bed, out_bed)
    >>> diff_two_files_identical(out_bed, out_exp_bed)
    True

    """
    # Check for bedtools.
    assert is_tool("bedtools"), "bedtools not in PATH"
    # Parameter string.
    params_str = '-s -c 4 -o distinct -delim ";"'
    if custom_params_str:
        params_str = custom_params_str
    check_cmd = "mergeBed -i " + in_bed + " " + params_str + " > " + out_bed
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "mergeBed is complaining:\n%s\n%s" %(check_cmd, output)


################################################################################

def bed_output_exon_12col_bed_file(tr_ids_dic, out_12col_bed,
                                   trid2reg_dic, trid2exc_dic, exid2reg_dic,
                                   exid2next_dic=False):

    """

    Output 12-column BED of transcript regions, for IGV viewing.

    tr_ids_dic:
        Transcript IDs to output dictionary.
    trid2reg_dic:
        Transcript ID to genomic region mapping
        transcript_id -> [chr_id, s, e, pol]
    trid2exc_dic:
        Transcript ID to exon count mapping
    exid2reg_dic:
        Exon ID to genomic region mapping
        exon_id -> [chr_id, s, e, pol]
    exid2next_dic:
        If exon ID to NEXT ID mapping is given, treat exid2reg_dic as
        NEXT -> exonic region mapping.

    """

    assert tr_ids_dic, "tr_ids_dic empty"

    OUT12BED = open(out_12col_bed, "w")

    for tr_id in tr_ids_dic:

        tr_chr_id = trid2reg_dic[tr_id][0]
        tr_gen_s = trid2reg_dic[tr_id][1] # 0-based.
        tr_gen_e = trid2reg_dic[tr_id][2]
        tr_gen_pol = trid2reg_dic[tr_id][3]
        tr_exc = trid2exc_dic[tr_id]

        ex_len_list = []
        ex_offset_list = []
        range_start = 0
        range_stop = tr_exc
        range_step = 1
        range_add = 1
        if tr_gen_pol == "-":
            range_start = tr_exc
            range_stop = 0
            range_step = -1
            range_add = 0

        for i in range(range_start, range_stop, range_step):
            ex_nr = i + range_add
            ex_id = tr_id + "_e" + str(ex_nr)
            ex_next = ex_id

            if exid2next_dic:
                assert ex_id in exid2next_dic, "exon ID %s not in exid2next_dic" %(ex_id)
                ex_next = exid2next_dic[ex_id]

            chr_id = exid2reg_dic[ex_next][0]
            assert tr_chr_id == chr_id, "transcript chromosome ID != NEXT region chromosome ID (%s != %s)" %(tr_chr_id, chr_id)
            gen_s = exid2reg_dic[ex_next][1] # 0-based.
            gen_e = exid2reg_dic[ex_next][2]
            gen_pol = exid2reg_dic[ex_next][3]
            assert tr_gen_pol == gen_pol, "transcript gene polarity != NEXT region polarity (%s != %s)" %(tr_gen_pol, gen_pol)

            ex_l = gen_e - gen_s
            ex_offset = gen_s - tr_gen_s
            ex_len_list.append(str(ex_l))
            ex_offset_list.append(str(ex_offset))

        # Output 12-col BED (IGV compatible).
        ex_len_str = ",".join(ex_len_list)
        ex_offset_str = ",".join(ex_offset_list)
        bed_out = "%s\t%i\t%i\t%s\t0\t%s\t%i\t%i\t100,100,100\t%i\t%s\t%s" %(tr_chr_id, tr_gen_s, tr_gen_e, tr_id, tr_gen_pol, tr_gen_s, tr_gen_e, tr_exc, ex_len_str, ex_offset_str)
        OUT12BED.write("%s\n" %(bed_out))
    OUT12BED.close()


################################################################################

def bed_get_score_to_count_dic(in_bed):
    """
    Given an .bed file in_bed, store scores and count how many times each
    score appears. Return dictionary with score -> count mapping.

    >>> in_bed = "test_data/test1.bed"
    >>> bed_get_score_to_count_dic(in_bed)
    {'1': 2, '0': 2, '2': 1, '3': 2}

    """
    assert os.path.isfile(in_bed), "cannot open in_bed \"%s\"" % (in_bed)
    # Read in IDs.
    sc2c_dic = {}
    with open(in_bed) as f:
        for line in f:
            row = line.strip()
            cols = line.strip().split("\t")
            site_sc = cols[4]
            if site_sc in sc2c_dic:
                sc2c_dic[site_sc] += 1
            else:
                sc2c_dic[site_sc] = 1
    f.closed
    return sc2c_dic


################################################################################

def gtf_extract_transcript_bed(in_gtf, out_bed,
                               trid2reg_dic=None,
                               chr_id_style=1,
                               tr_ids_dic=False):
    """
    Extract transcript regions from in_gtf GTF file, and output to out_bed BED
    file.

    tr_ids_dic:
        Dictionary with transcript IDs for filtering (keeping dic IDs).
    trid2reg_dic:
        Store genomic region of transcript [chr_id, s, e, pol]

    >>> in_gtf = "test_data/gene_test_in.gtf"
    >>> exp_out_bed = "test_data/gtf_transcript_out.exp.bed"
    >>> tmp_out_bed = "test_data/gtf_transcript_out.tmp.bed"
    >>> gtf_extract_transcript_bed(in_gtf, tmp_out_bed)
    >>> diff_two_files_identical(tmp_out_bed, exp_out_bed)
    True

    """

    # Output transcript regions.
    OUTBED = open(out_bed, "w")
    c_out = 0
    # Open GTF either as .gz or as text file.
    if re.search(".+\.gz$", in_gtf):
        f = gzip.open(in_gtf, 'rt')
    else:
        f = open(in_gtf, "r")
    for line in f:
        # Skip header.
        if re.search("^#", line):
            continue
        cols = line.strip().split("\t")
        chr_id = cols[0]
        feature = cols[2]
        feat_s = int(cols[3])
        feat_e = int(cols[4])
        feat_pol = cols[6]
        infos = cols[8]
        if not feature == "transcript":
            continue

        # Restrict to standard chromosomes.
        new_chr_id = check_convert_chr_id(chr_id,
                                id_style=chr_id_style)
        if not new_chr_id:
            continue
        else:
            chr_id = new_chr_id

        # Make start coordinate 0-base (BED standard).
        feat_s = feat_s - 1

        # Extract transcript ID.
        m = re.search('transcript_id "(.+?)"', infos)
        assert m, "transcript_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        transcript_id = m.group(1)

        # Check if transcript ID is in transcript dic.
        if tr_ids_dic:
            if not transcript_id in tr_ids_dic:
                continue

        if trid2reg_dic is not None:
            trid2reg_dic[transcript_id] = [chr_id, feat_s, feat_e, feat_pol]

        # Output genomic exon region.
        c_out += 1
        OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,feat_s,feat_e,transcript_id,feat_pol))

    OUTBED.close()
    f.close()

    assert c_out, "no regions output to out_bed. Invalid in_gtf or too restrictive tr_ids_dic filtering?"


################################################################################

def ph_extract_generate_html_report(out_folder, hoodlib_path,
                            eir_stats_dic=False,
                            ei_ratios_list=False,
                            eib_ratios_list=False,
                            intergen_site_ids_dic=False,
                            intron_site_ids_dic=False,
                            all_ex_site_ids_dic=False,
                            tc_ex_site_ids_dic=False,
                            gc_ex_site_ids_dic=False,
                            all_tr_site_seqs_dic=False,
                            sel_tr_site_seqs_dic=False,
                            gen_rr_ratios_dic=False,
                            tr_rr_ratios_dic=False,
                            site_lengths_list=False,
                            copy_logo=True,
                            html_report_out="report.peakhood_extract.html",
                            plots_subfolder="html_plots"):
    """
    HTML report for peakhood extract.


    eir_stats_dic:
        Various exon-intron ratio stats dictionary.
    ei_ratios_list:
        Exon-intron ratios list.
    eib_ratios_list:
        Exon-intron border region ratios list.
    all_ex_site_ids_dic:
        All exonic site IDs (ID -> sequence)
    intron_site_ids_dic:
        Intronic site IDs (ID -> sequence)
    intergen_site_ids_dic:
        Intergenic site IDs (ID -> sequence)
    tc_ex_site_ids_dic:
        Transcript context exonic site IDs (ID -> sequence)
    gc_ex_site_ids_dic:
        Genomic context exonic site IDs (ID -> sequence)
    all_tr_site_seqs_dic:
        All transcript context site_id,tr_id combinations (ID -> sequence)
    sel_tr_site_seqs_dic:
        Selected transcript context site_id,tr_id combinations (ID -> sequence)
    gen_rr_ratios_dic:
        Genomic site to repeat region ratio
        (repeat region nucleotides / all site nucleotides)
        site ID -> ratio
    tr_rr_ratios_dic:
        Transcript site to repeat region ratio
        (repeat region nucleotides / all site nucleotides)
        site_id,tr_id combination ID -> ratio
    copy_logo:
        Copy logo to results plots folder.

    """

    assert os.path.exists(out_folder), "out_folder does not exist"
    assert os.path.exists(hoodlib_path), "hoodlib_path does not exist"
    assert eir_stats_dic, "eir_stats_dic empty"
    assert ei_ratios_list, "ei_ratios_list empty"
    assert eib_ratios_list, "eib_ratios_list empty"
    assert site_lengths_list, "site_lengths_list empty"

    from markdown import markdown

    plots_folder = plots_subfolder
    plots_out_folder = out_folder + "/" + plots_folder
    if not os.path.exists(plots_out_folder):
        os.makedirs(plots_out_folder)
    html_out = out_folder + "/" + "report.peakhood_extract.html"
    if html_report_out:
        html_out = html_report_out

    # Plot files.
    ei_ratio_density_plot = "ei_ratio_density_plot.png"
    ei_ratio_density_plot_out = plots_out_folder + "/" + ei_ratio_density_plot
    eib_ratio_density_plot = "eib_ratio_density_plot.png"
    eib_ratio_density_plot_out = plots_out_folder + "/" + eib_ratio_density_plot
    lengths_plot = "set_lengths_plot.png"
    lengths_plot_out = plots_out_folder + "/" + lengths_plot

    # Paths to logo and ploty.
    logo1_path = hoodlib_path + "/content/logo1.png"
    logo2_path = hoodlib_path + "/content/logo2.png"
    # sorttable_js_path = hoodlib_path + "/content/sorttable.js"
    # plotly_js_path = hoodlib_path + "/content/plotly-latest.min.js"
    # assert os.path.exists(plotly_js_path), "plotly js %s not found" %(plotly_js_path)

    # Copy logo to plots folder.
    if copy_logo:
        logo_out = plots_out_folder + "/" + "logo.png"
        shutil_copy_file(logo1_path, logo_out)
        logo1_path = plots_folder + "/" + "logo.png"

    mdtext = """
<head>
<title>Peakhood - Context Extraction Report</title>
</head>

<img src="%s" alt="ph_logo"
	title="ph_logo" width="650" />

<p><body style="font-family:sans-serif" link="#007af4" vlink="#007af4" alink="#007af4"></p>

""" %(logo1_path)

    mdtext += """

# Context Extraction Report

List of available context extraction statistics generated
by Peakhood (peakhood extract):

- [Input site length statistics](#site-length-stats)
- [Input site length distribution](#site-length-plot)
- [Site region type statistics](#site-region-stats)
- [Exon-intron coverage ratio statistics](#ei-ratio-stats)
- [Exon-intron coverage ratio distribution](#ei-ratio-plot)
- [Exon-intron border coverage ratio statistics](#eib-ratio-stats)
- [Exon-intron border coverage ratio distribution](#eib-ratio-plot)
- [Repeat region statistics](#rep-reg-stats)"""


    """
    Site lengths statistics.

    """

    mdtext += """
## Input site length statistics ### {#site-length-stats}

**Table:** Input site length statistics
(min, max, mean, and median length) in nucleotides (nt).


"""

    mdtext += "| Attribute | &nbsp; Value &nbsp; | \n"
    mdtext += "| :-: | :-: |\n"
    mdtext += "| # sites | %i |\n" %(len(site_lengths_list))
    mdtext += "| min site length | %i |\n" %(min(site_lengths_list))
    mdtext += "| max site length | %i |\n" %(max(site_lengths_list))
    mdtext += "| mean site length | %.1f |\n" %(statistics.mean(site_lengths_list))
    mdtext += "| median site length | %i |\n" %(statistics.median(site_lengths_list))
    mdtext += "\n&nbsp;\n&nbsp;\n"

    # Make site length distribution box plot.
    create_site_lengths_plot(site_lengths_list, lengths_plot_out)
    lengths_plot_path = plots_folder + "/" + lengths_plot

    mdtext += """
## Input site length distribution ### {#site-length-plot}

Input site length distribution, after pre-merging of book-ended and
overlapping input sites (if set)
and pre-filtering (if set, e.g. by score or length).
Note that set --pre-merge leads to increased lengths if there are adjacent
or overlapping sites. Moreover, set --max-len (default 200) limits the
maximum site length, but this can increase again if --pre-merge is set
(since --pre-merge is applied after --max-len filtering).

"""
    mdtext += '<img src="' + lengths_plot_path + '" alt="Site length distribution"' + "\n"
    mdtext += 'title="Site length distribution" width="500" />' + "\n"
    mdtext += """

**Figure:** Input site length distribution (after pre-filtering and-pre-merging sites).

&nbsp;

"""

    """
    Site region type statistics

    To add ?
    - [Site region type distribution](#site-region-plot)

    """

    c_intergen_sites = len(intergen_site_ids_dic)
    c_intron_sites = len(intron_site_ids_dic)
    c_exon_sites_merged_tc = len(sel_tr_site_seqs_dic)
    c_exon_sites_tc = len(tc_ex_site_ids_dic)
    c_exon_sites_gc = len(gc_ex_site_ids_dic)

    mdtext += """
## Site region type statistics ### {#site-region-stats}

**Table:**
Assigned site region type statistics. Exonic sites can be either assigned to
transcript context (TC) or genomic context (GC), depending on the read information
in the input BAM file. In addition, transcript context site count with
merged exon border sites is given (MEXB).
Intronic sites are sites with insufficient
(see --min-exon-overlap, default >= 90 percent) or no overlap with any
exonic region from the input GTF file. Intergenic sites do not overlap
with any transcript regions from the input GTF file. Depending on which
pipeline was used to determine the input CLIP-seq peak regions, there might
be little or no intergenic sites due to pre-filtering for gene regions.

"""
    mdtext += "| &nbsp; Region type &nbsp; | &nbsp; &nbsp; Count &nbsp; &nbsp; | \n"
    mdtext += "| :-: | :-: |\n"
    mdtext += "| exonic (GC) | %i |\n" %(c_exon_sites_gc)
    mdtext += "| exonic (TC) |  %i |\n" %(c_exon_sites_tc)
    mdtext += "| exonic (TC) MEXB |  %i |\n" %(c_exon_sites_merged_tc)
    mdtext += "| intronic | %i |\n" %(c_intron_sites)
    mdtext += "| intergenic |  %i |\n" %(c_intergen_sites)
    mdtext += "\n&nbsp;\n&nbsp;\n"


    """
    Exon-intron coverage ratio statistics
    """

    # EIR stats.
    c_uniq_exons = eir_stats_dic["c_uniq_exons"]
    eir_mean = eir_stats_dic["eir_mean"]
    eir_median = eir_stats_dic["eir_median"]
    eir_stdev = eir_stats_dic["eir_stdev"]
    eir_perc5 = eir_stats_dic["eir_perc5"]
    eir_perc25 = eir_stats_dic["eir_perc25"]
    eir_perc50 = eir_stats_dic["eir_perc50"]
    eir_perc75 = eir_stats_dic["eir_perc75"]
    eir_perc95 = eir_stats_dic["eir_perc95"]
    eir_min = eir_stats_dic["eir_min"]
    eir_max = eir_stats_dic["eir_max"]

    mdtext += """
## Exon-intron coverage ratio statistics ### {#ei-ratio-stats}

**Table:**
Exon-intron coverage ratios statistics for unique exon regions (# unique exons: %i)
containing CLIP-seq sites. A unique exon region can include several annotated
exon regions, since GTF files usually contain exons with different IDs but
identical regions.
The unique exon region ratio is the average ratio of all exon regions with
the same coordinates as the unique exon region.
The ratio of an exon region is calculated
by dividing the exon coverage (reads / region length) through the coverage
of the neighboring intron(s). In case of two introns, the average coverage
of the two introns is used as the divisor. In case of no introns, a fixed
value above the threshold is assigned.

""" %(c_uniq_exons)

    mdtext += "| &nbsp; Attribute &nbsp; | &nbsp; &nbsp; Value &nbsp; &nbsp; | \n"
    mdtext += "| :-: | :-: |\n"
    mdtext += "| # unique exons | %i |\n" %(c_uniq_exons)
    mdtext += "| min ratio | %.4f |\n" %(eir_min)
    mdtext += "| max ratio | %.4f |\n" %(eir_max)
    mdtext += "| mean ratio | %.4f |\n" %(eir_mean)
    mdtext += "| stdev ratio | %.4f |\n" %(eir_stdev)
    mdtext += "| median ratio | %.4f |\n" %(eir_median)
    mdtext += "| 25th percentile ratio | %.4f |\n" %(eir_perc25)
    mdtext += "| 50th percentile ratio | %.4f |\n" %(eir_perc50)
    mdtext += "| 75th percentile ratio | %.4f |\n" %(eir_perc75)
    mdtext += "\n&nbsp;\n&nbsp;\n"

    kde_s = eir_min
    kde_e = eir_perc95
    kde_clip = [kde_s, kde_e]
    kde_bw_adjust = 0.4
    plot_ei_ratio_density(ei_ratios_list, ei_ratio_density_plot_out,
                          x_label="Exon-intron coverage ratio",
                          y_label="Density",
                          fig_width=7,
                          fig_height=3,
                          kde_bw_adjust=kde_bw_adjust,
                          kde_clip=kde_clip,
                          x_0_to_100=False)

    plot_path = plots_folder + "/" + ei_ratio_density_plot

    mdtext += """
## Exon-intron coverage ratio distribution ### {#ei-ratio-plot}

This plot shows the distribution of exon-intron coverage ratios for
unique exon regions containing CLIP-seq sites.

"""
    mdtext += '<img src="' + plot_path + '" alt="ei_raio_plot"' + "\n"
    mdtext += 'title="Exon-intron ratio distribution" width="650" />' + "\n"
    mdtext += """

**Figure:** Distribution of exon-intron ratios (exon coverage divided by
surrounding intron coverage) for unique exon regions containing CLIP-seq sites.
To prevent suboptimal scaling due to outliers, ratios are plotted only up
to the 95th percentile ratio.

&nbsp;

"""

    """
    Exon-intron border coverage ratio statistics
    """

    # EIBR stats.
    c_eibr_regions = eir_stats_dic["c_eibr_regions"]
    eibr_mean = eir_stats_dic["eibr_mean"]
    eibr_median = eir_stats_dic["eibr_median"]
    eibr_stdev = eir_stats_dic["eibr_stdev"]
    eibr_perc5 = eir_stats_dic["eibr_perc5"]
    eibr_perc25 = eir_stats_dic["eibr_perc25"]
    eibr_perc50 = eir_stats_dic["eibr_perc50"]
    eibr_perc75 = eir_stats_dic["eibr_perc75"]
    eibr_perc95 = eir_stats_dic["eibr_perc95"]
    eibr_min = eir_stats_dic["eibr_min"]
    eibr_max = eir_stats_dic["eibr_max"]

    mdtext += """
## Exon-intron border coverage ratio statistics ### {#eib-ratio-stats}

**Table:**
Exon-intron border coverage ratio statistics for unique exon regions
containing CLIP-seq sites. A unique exon region can include several annotated
exon regions, since GTF files usually contain exons with different IDs but
identical regions. Note that not all unique exon regions might be considered,
since exon-intron borders with small read coverages are not considered.
The ratio is calculated for each exon-intron border, taking a small border
region on the intron as well as on the exon, and calculating the coverage ratio
between the two. This is done for both exon ends, and the average or the
ratio with more reads is returned for each exon region. These ratios are then
merged to one ratio for each unique exon region.

"""

    mdtext += "| &nbsp; Attribute &nbsp; | &nbsp; &nbsp; Value &nbsp; &nbsp; | \n"
    mdtext += "| :-: | :-: |\n"
    mdtext += "| # considered exons | %i |\n" %(c_eibr_regions)
    mdtext += "| min ratio | %.4f |\n" %(eibr_min)
    mdtext += "| max ratio | %.4f |\n" %(eibr_max)
    mdtext += "| mean ratio | %.4f |\n" %(eibr_mean)
    mdtext += "| stdev ratio | %.4f |\n" %(eibr_stdev)
    mdtext += "| median ratio | %.4f |\n" %(eibr_median)
    mdtext += "| 25th percentile ratio | %.4f |\n" %(eibr_perc25)
    mdtext += "| 50th percentile ratio | %.4f |\n" %(eibr_perc50)
    mdtext += "| 75th percentile ratio | %.4f |\n" %(eibr_perc75)
    mdtext += "\n&nbsp;\n&nbsp;\n"

    kde_s = eibr_min
    kde_e = eibr_perc95
    kde_clip = [kde_s, kde_e]
    kde_bw_adjust = 0.4
    plot_ei_ratio_density(eib_ratios_list, eib_ratio_density_plot_out,
                          x_label="Exon-intron border coverage ratio",
                          y_label="Density",
                          fig_width=7,
                          fig_height=3,
                          kde_bw_adjust=kde_bw_adjust,
                          kde_clip=kde_clip,
                          x_0_to_100=False)

    plot_path = plots_folder + "/" + eib_ratio_density_plot

    mdtext += """
## Exon-intron border coverage ratio distribution ### {#ei-ratio-plot}

This plot shows the distribution of exon-intron border coverage ratios for
 unique exon regions containing CLIP-seq sites.

"""
    mdtext += '<img src="' + plot_path + '" alt="ei_raio_plot"' + "\n"
    mdtext += 'title="Exon-intron border ratio distribution" width="650" />' + "\n"
    mdtext += """

**Figure:** Distribution of exon-intron border ratios (exon border coverages divided by
adjacent intron border coverages) for unique exon regions containing CLIP-seq sites.
To prevent suboptimal scaling due to outliers, ratios are plotted only up
to the 95th percentile ratio.

&nbsp;

"""

    """
    Repeat region coverage statistics

    """

    # Exonic sites, selected transcript context, transcript.
    sel_exs_tc_tr_rrc = 0.0
    sel_tc_tr_c = len(sel_tr_site_seqs_dic)
    for sitetrid in sel_tr_site_seqs_dic:
        sel_exs_tc_tr_rrc += tr_rr_ratios_dic[sitetrid]
    sel_exs_tc_tr_perc = 0.0
    if sel_tc_tr_c:
        sel_exs_tc_tr_perc = (sel_exs_tc_tr_rrc / sel_tc_tr_c) * 100

    # Exonic sites, genomic context, genomic.
    exs_gc_gen_rrc = 0.0
    exs_gc_gen_c = len(gc_ex_site_ids_dic)
    for site_id in gc_ex_site_ids_dic:
        exs_gc_gen_rrc += gen_rr_ratios_dic[site_id]
    exs_gc_gen_perc = 0.0
    if exs_gc_gen_c:
        exs_gc_gen_perc = (exs_gc_gen_rrc / exs_gc_gen_c) * 100

    # Intronic sites.
    intronic_rrc = 0.0
    intronic_c = len(intron_site_ids_dic)
    for site_id in intron_site_ids_dic:
        intronic_rrc += gen_rr_ratios_dic[site_id]
    intronic_perc = 0.0
    if intronic_c:
        intronic_perc = (intronic_rrc / intronic_c) * 100

    # Intergenic sites.
    intergen_rrc = 0.0
    intergen_c = len(intergen_site_ids_dic)
    for site_id in intergen_site_ids_dic:
        intergen_rrc += gen_rr_ratios_dic[site_id]
    intergen_perc = 0.0
    if intergen_c:
        intergen_perc = (intergen_rrc / intergen_c) * 100

    mdtext += """
## Repeat region statistics ### {#rep-reg-stats}

**Table:**
Repeat region content statistics for different region types.
The percentage of repeat regions found in each region type set is given.

"""
    mdtext += "| &nbsp; Region type &nbsp; | &nbsp; Count &nbsp; | &nbsp; &nbsp; Percentage &nbsp; &nbsp; | \n"
    mdtext += "| :-: | :-: | :-: |\n"
    mdtext += "| exonic (GC) | %i | %.3f |\n" %(exs_gc_gen_c, exs_gc_gen_perc)
    mdtext += "| exonic (TC) | %i | %.3f |\n" %(sel_tc_tr_c, sel_exs_tc_tr_perc)
    mdtext += "| intronic    | %i | %.3f |\n" %(intronic_c, intronic_perc)
    mdtext += "| intergenic  | %i | %.3f |\n" %(intergen_c, intergen_perc)
    mdtext += "\n&nbsp;\n&nbsp;\n"

    # Convert mdtext to html.
    md2html = markdown(mdtext, extensions=['attr_list', 'tables'])
    OUTHTML = open(html_out,"w")
    OUTHTML.write("%s\n" %(md2html))
    OUTHTML.close()
    # change <table> to sortable.
    check_cmd = "sed -i 's/<table>/<table class=" + '"sortable"' + ">/g' " + html_out
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "sed command returned error:\n%s" %(output)


################################################################################

def create_site_lengths_plot(site_len_list, out_plot,
                             scale_zero_max=False):
    """
    Create a box plot, showing the distribution of site lengths.

    Peakhood colors:
    #b237fd  purple
    #007af4  blue
    #00b4f7  light blue
    #00d8ff  lighter blue

    """
    # Checker.
    assert site_len_list, "given list site_len_list empty"
    if scale_zero_max:
        # Get maximum length for scaling.
        max_l = max(site_len_list)
        # Get next highest number % 10.
        max_y = max_l
        while max_y % 10:
             max_y += 1
    # Make pandas dataframe.
    test_label = "Input sites"
    data = {'set': [], 'length': []}
    test_c = len(site_len_list)
    data['set'] += test_c*[test_label]
    data['length'] += site_len_list
    df = pd.DataFrame (data, columns = ['set','length'])

    # Make plot.
    sns.set(style="darkgrid")
    fig, ax = plt.subplots()
    sns.boxplot(x="set", y="length", data=df, palette=['#00b4f7'],
                width=0.7, linewidth = 1.5, boxprops=dict(alpha=.7))
    # Modify.
    ax.set_ylabel("Length (nt)",fontsize=18)
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=12)
    if scale_zero_max:
        ax.set_ylim([0,max_y])
    ax.set(xlabel=None)
    # Store plot.
    fig.savefig(out_plot, dpi=125, bbox_inches='tight')


################################################################################

def plot_ei_ratio_density(set_scores, out_plot,
                          set_label="Positives",
                          x_label="k-mer score",
                          y_label="Density",
                          fig_width=5,
                          fig_height=4,
                          kde_bw_adjust=1,
                          x_0_to_100=False,
                          kde_clip=False):
    """
    Exon-intron ratios distribution plot.


    PH graffiti colors:
    #b237fd : pink
    #007af4 : blue
    #00d7fa : light blue
    #0bf80b : slime green


    """
    assert set_scores, "set_scores empty"
    if not kde_clip:
        max_sc = max(set_scores)
        min_sc = min(set_scores)
        kde_clip = [min_sc, max_sc]

    data = {'score': []}
    data['score'] += set_scores
    df = pd.DataFrame (data, columns = ['score'])

    # Make plot.
    sns.set(style="darkgrid")
    fig, ax = plt.subplots()
    sns.kdeplot(x="score", data=df, color='#b237fd',
                clip=kde_clip, bw_adjust=kde_bw_adjust)
    fig.set_figwidth(fig_width)
    fig.set_figheight(fig_height)
    ax.set(xlabel=x_label)
    ax.set_ylabel(y_label)
    #ax.tick_params(axis='x', labelsize=18)
    #ax.tick_params(axis='y', labelsize=14)
    fig.savefig(out_plot, dpi=150, bbox_inches='tight')


################################################################################

def ph_merge_generate_html_report(out_folder, hoodlib_path,
                            ex_sites_perc_dic=False,
                            tc_sites_perc_dic=False,
                            exb_sites_perc_dic=False,
                            add_stats_dic=False,
                            set_stats_dd=False,
                            set2site_len_dic=False,
                            copy_logo=True,
                            html_report_out="report.peakhood_merge.html",
                            plots_subfolder="html_plots"):
    """
    HTML report for peakhood merge.

    """

    assert os.path.exists(out_folder), "out_folder does not exist"
    assert os.path.exists(hoodlib_path), "hoodlib_path does not exist"
    assert ex_sites_perc_dic, "ex_sites_perc_dic empty"
    assert tc_sites_perc_dic, "tc_sites_perc_dic empty"
    assert exb_sites_perc_dic, "exb_sites_perc_dic empty"
    assert add_stats_dic, "add_stats_dic empty"
    assert set_stats_dd, "set_stats_dd empty"
    assert set2site_len_dic, "set2site_len_dic empty"

    from markdown import markdown

    plots_folder = plots_subfolder
    plots_out_folder = out_folder + "/" + plots_folder
    if not os.path.exists(plots_out_folder):
        os.makedirs(plots_out_folder)
    html_out = out_folder + "/" + "report.peakhood_merge.html"
    if html_report_out:
        html_out = html_report_out

    # Plot files.
    exon_perc_plot = "exon_perc_plot.png"
    exon_perc_plot_out = plots_out_folder + "/" + exon_perc_plot
    lengths_plot = "set_lengths_plot.png"
    lengths_plot_out = plots_out_folder + "/" + lengths_plot

    # Paths to logo and ploty.
    logo1_path = hoodlib_path + "/content/logo1.png"
    logo2_path = hoodlib_path + "/content/logo2.png"
    # sorttable_js_path = hoodlib_path + "/content/sorttable.js"
    # plotly_js_path = hoodlib_path + "/content/plotly-latest.min.js"
    # assert os.path.exists(plotly_js_path), "plotly js %s not found" %(plotly_js_path)

    # Copy logo to plots folder.
    if copy_logo:
        logo_out = plots_out_folder + "/" + "logo.png"
        shutil_copy_file(logo1_path, logo_out)
        logo1_path = plots_folder + "/" + "logo.png"

    mdtext = """
<head>
<title>Peakhood - Merge Extracted Datasets Report</title>
</head>

<img src="%s" alt="ph_logo"
	title="ph_logo" width="675" />

<p><body style="font-family:sans-serif" link="#007af4" vlink="#007af4" alink="#007af4"></p>

""" %(logo1_path)


    mdtext += """

# Merge Extracted Datasets Report

List of available statistics generated
by Peakhood (peakhood merge):

- [Merged dataset statistics](#merge-stats)
- [Site region type statistics](#exon-perc-stats)
- [Input site length distribution](#site-length-plot)
- [Exonic site percentages distribution](#exon-perc-plot)"""


    """
    Merged dataset statistics.

    """

    c_all_tr = add_stats_dic["c_all_tr"]
    c_sel_tr = add_stats_dic["c_sel_tr"]
    c_pair_sites_all_tr = add_stats_dic["c_pair_sites_all_tr"]
    c_pair_sites_all_tr_diff_set = add_stats_dic["c_pair_sites_all_tr_diff_set"]
    c_pair_sites_sel_tr = add_stats_dic["c_pair_sites_sel_tr"]
    c_pair_sites_sel_tr_diff_set = add_stats_dic["c_pair_sites_sel_tr_diff_set"]

    mdtext += """
## Merged dataset statistics ### {#merge-stats}

**Table:**
Merged dataset statistics, listing over all datasets: the number of
transcripts containing sites, the number of selected transcripts (most likely
transcripts from peakhood extract) with sites, the number of site pairs on
all transcripts (same and different datasets / RBPs), the number of site pairs
from different datasets (different RBPs), the number of site pairs
on the selected transcripts (same and different RBPs), and the number of
site pairs on the selected transcripts (different RBPs).

"""

    mdtext += "| &nbsp; Description &nbsp; | &nbsp; &nbsp; Count &nbsp; &nbsp; | \n"
    mdtext += "| :-: | :-: |\n"
    mdtext += "| # all transcripts with sites | %i |\n" %(c_all_tr)
    mdtext += "| # selected transcripts with sites |  %i |\n" %(c_sel_tr)
    mdtext += "| # site pairs on all transcripts |  %i |\n" %(c_pair_sites_all_tr)
    mdtext += "| # site pairs (from different datasets) | %i |\n" %(c_pair_sites_all_tr_diff_set)
    mdtext += "| # site pairs on selected transcripts |  %i |\n" %(c_pair_sites_sel_tr)
    mdtext += "| # site pairs (from different datasets) | %i |\n" %(c_pair_sites_sel_tr_diff_set)
    mdtext += "\n&nbsp;\n&nbsp;\n"


    """
    Site region type statistics.

    """

    mdtext += """
## Site region type statistics ### {#exon-perc-stats}

**Table:**
Site region type statistics. For each input dataset,
different region types with corresponding site numbers are given:
number of all dataset sites,
number of intronic sites, number of intergenic sites,
number of exonic sites with assigned genomic context,
number of exonic sites with assigned transcript context (TC),
number of exonic sites with assigned transcript context after merging
adjacent exon border sites (TCM),
number of exonic sites at exon borders connected by intron-spanning reads (before merging),
percentage of exonic sites (exonic sites / all sites), and
percentage of exonic transcript context sites (TC sites / all exonic sites).

"""

    mdtext += "| &nbsp; Dataset &nbsp; | # all sites | # intronic | # intergenic | # exonic (GC) | # exonic (TC) | # exonic (TCM) | # exon border | % exon / all | % TC / exonic | \n"
    mdtext += "| :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | \n"

    for plot_data_id in set_stats_dd:

        c_all_sites = set_stats_dd[plot_data_id]["c_all_sites"]
        c_intron_sites = set_stats_dd[plot_data_id]["c_intron_sites"]
        c_intergen_sites = set_stats_dd[plot_data_id]["c_intergen_sites"]
        c_exon_sites_gc = set_stats_dd[plot_data_id]["c_exon_sites_gc"]
        c_exon_sites_tc = set_stats_dd[plot_data_id]["c_exon_sites_tc"]
        c_exon_sites_merged_tc = set_stats_dd[plot_data_id]["c_exon_sites_merged_tc"]
        c_exonic_sites = set_stats_dd[plot_data_id]["c_exonic_sites"]
        c_exb_sites = set_stats_dd[plot_data_id]["c_exb_sites"]
        perc_exonic_sites = 0.0
        perc_exon_sites_tc = 0.0
        if c_exonic_sites and c_all_sites:
            perc_exonic_sites = (c_exonic_sites / c_all_sites) * 100
        if c_exon_sites_tc and c_exonic_sites:
            perc_exon_sites_tc = (c_exon_sites_tc / c_exonic_sites) * 100

        mdtext += "| %s | %i | %i | %i | %i | %i | %i | %i | %.2f | %.2f | \n" %(plot_data_id, c_all_sites, c_intron_sites, c_intergen_sites, c_exon_sites_gc, c_exon_sites_tc, c_exon_sites_merged_tc, c_exb_sites, perc_exonic_sites, perc_exon_sites_tc)

    mdtext += "\n&nbsp;\n&nbsp;\n"

    create_merge_site_lengths_plot(set2site_len_dic, lengths_plot_out)
    lengths_plot_path = plots_folder + "/" + lengths_plot

    mdtext += """
## Input site length distribution ### {#site-length-plot}

Input site length distribution for every --in input dataset,
after pre-merging of book-ended and overlapping input sites (if set)
and pre-filtering (if set, e.g. by score or length).
Note that set --pre-merge leads to increased lengths if there are adjacent
or overlapping sites. Moreover, set --max-len (default 200) limits the
maximum site length, but this can increase again if --pre-merge is set
(since --pre-merge is applied after --max-len filtering).

"""
    mdtext += '<img src="' + lengths_plot_path + '" alt="Site length distribution"' + "\n"
    mdtext += 'title="Site length distribution" width="700" />' + "\n"
    mdtext += """

**Figure:** Input site length distribution for every --in input dataset.

&nbsp;

"""


    mdtext += """
## Exonic site percentages distribution ### {#exon-perc-plot}

Percentages of exonic sites (exonic sites / all sites), assigned
transcript context (TC) sites (TC / exonic sites), and
paired exonic sites at exon borders (EXBS) (EXBS / TC),
for all --in input datasets. Two sites at exon borders form a pair
if they are connected via intron-spanning reads,
and thus likely form one single site instead of two separate.
Pair sites at exon borders consequently
get merged by Peakhood, so usually only 1 of 2 sites remains.
Moreover, if there is > 1 connection for a set of exon border sites,
Peakhood only keeps the one featuring the most intron-spanning reads.

"""
    create_exon_perc_plot(ex_sites_perc_dic, tc_sites_perc_dic,
                          exb_sites_perc_dic, exon_perc_plot_out)
    exon_perc_plot_path = plots_folder + "/" + exon_perc_plot
    mdtext += '<img src="' + exon_perc_plot_path + '" alt="Exonic site percentages distribution"' + "\n"
    mdtext += 'title="Exonic site percentages distribution" width="800" />' + "\n"
    mdtext += """
**Figure:** Percentages of exonic sites (exonic sites / all sites),
assigned transcript context (TC) sites (TC / exonic sites),
and exonic sites at exon borders connected by intron-spanning reads
(EXBS) (EXBS / TC),
for all --in input datasets. The higher the first two percentages,
the more likely the RBP associated to the dataset binds to a
spliced context (introns removed).

&nbsp;

"""
    # Convert mdtext to html.
    md2html = markdown(mdtext, extensions=['attr_list', 'tables'])
    OUTHTML = open(html_out,"w")
    OUTHTML.write("%s\n" %(md2html))
    OUTHTML.close()
    # change <table> to sortable.
    # check_cmd = "sed -i 's/<table>/<table class=" + '"sortable"' + ">/g' " + html_out
    # output = subprocess.getoutput(check_cmd)
    # error = False
    # if output:
    #     error = True
    # assert error == False, "sed command returned error:\n%s" %(output)


################################################################################

def create_merge_site_lengths_plot(set2site_len_dic, out_plot):
    """
    Create a box plot, showing the site length distributions of
    input datasets.

    Peakhood colors:
    #b237fd  purple
    #007af4  blue
    #00b4f7  light blue
    #00d8ff  lighter blue

    """

    # Checker.
    assert set2site_len_dic, "given dictionary set2site_len_dic empty"

    data = {'set': [], 'length': []}
    for set_id in set2site_len_dic:
        c_sites = len(set2site_len_dic[set_id])
        data['set'] += c_sites*[set_id]
        data['length'] += set2site_len_dic[set_id]
    df = pd.DataFrame (data, columns = ['set','length'])

    # Scale height depending on # of sets.
    c_ids = len(set2site_len_dic)
    fheight = 1.5 * c_ids

    # Make plot.
    sns.set(style="darkgrid")

    # g = sns.boxplot(data=df, orient="h", palette=["#00b4f7"],
    #                 y="set", x="length",
    #                 width=0.7, linewidth = 1.5, boxprops=dict(alpha=.7))
    # g = sns.catplot(x="perc", y="feat_id", hue="set", data=df,
    #                 kind="bar", palette=["#00d8ff", "#007af4", "#b237fd"],
    #                 legend=False)

    g = sns.catplot(x="length", y="set", data=df,
                    kind="box", palette=["#00d8ff"],
                    legend=False)

    g.fig.set_figwidth(18)
    g.fig.set_figheight(fheight)
    # Modify axes.
    ax = g.axes
    ax[0,0].set_xlabel("Site length distribution",fontsize=20)
    ax[0,0].set(ylabel=None)
    ax[0,0].tick_params(axis='x', labelsize=16)
    ax[0,0].tick_params(axis='y', labelsize=20)
    # Add legend at specific position.
    #plt.legend(loc=(1.01, 0.4), fontsize=16)
    g.savefig(out_plot, dpi=100, bbox_inches='tight')


################################################################################

def create_exon_perc_plot(ex_sites_perc_dic, tc_sites_perc_dic,
                          exb_sites_perc_dic, out_plot):
    """
    Create a grouped bar plot, showing the percentages of exonic sites
    and transcript context sites for each dataset.

    Peakhood colors:
    #b237fd  purple
    #007af4  blue
    #00b4f7  light blue
    #00d8ff  lighter blue

    """

    # Checker.
    assert ex_sites_perc_dic, "given dictionary ex_sites_perc_dic empty"
    assert tc_sites_perc_dic, "given dictionary tc_sites_perc_dic empty"
    assert exb_sites_perc_dic, "given dictionary exb_sites_perc_dic empty"

    # Make pandas dataframe.
    ex_label = "exonic / all"
    tc_label = "TC / exonic"
    exb_label = "EXBS / TC"
    data = {'set': [], 'feat_id': [], 'perc': []}

    # feat_id: RBP / dataset name.
    # set: "TC", "exonic"
    for feat_id in ex_sites_perc_dic:
        data['set'].append(ex_label)
        data['feat_id'].append(feat_id)
        data['perc'].append(ex_sites_perc_dic[feat_id])
    for feat_id in tc_sites_perc_dic:
        data['set'].append(tc_label)
        data['feat_id'].append(feat_id)
        data['perc'].append(tc_sites_perc_dic[feat_id])
    for feat_id in exb_sites_perc_dic:
        data['set'].append(exb_label)
        data['feat_id'].append(feat_id)
        data['perc'].append(exb_sites_perc_dic[feat_id])
    df = pd.DataFrame (data, columns = ['set','feat_id', 'perc'])

    # Scale height depending on # of features.
    c_ids = len(ex_sites_perc_dic)
    fheight = 1.5 * c_ids

    # Make plot.
    sns.set(style="darkgrid")
    g = sns.catplot(x="perc", y="feat_id", hue="set", data=df,
                    kind="bar", palette=["#00d8ff", "#007af4", "#b237fd"],
                    legend=False)
    g.fig.set_figwidth(18)
    g.fig.set_figheight(fheight)
    # Modify axes.
    ax = g.axes
    ax[0,0].set_xlabel("Percentage of sites (%)",fontsize=20)
    ax[0,0].set(ylabel=None)
    ax[0,0].tick_params(axis='x', labelsize=16)
    ax[0,0].tick_params(axis='y', labelsize=20)
    # Add legend at specific position.
    plt.legend(loc=(1.01, 0.4), fontsize=16)
    g.savefig(out_plot, dpi=100, bbox_inches='tight')


################################################################################

def print_some_banner():
    """
    Print some banner.

    """
    banner = []

    a = """
                               $$\\
                               $$ |
  $$$$$$\   $$$$$$\   $$$$$$\  $$ |  $$\\
 $$  __$$\ $$  __$$\  \____$$\ $$ | $$  |
 $$ /  $$ |$$$$$$$$ | $$$$$$$ |$$$$$$  /
 $$ |  $$ |$$   ____|$$  __$$ |$$  _$$<
 $$$$$$$  |\$$$$$$$\ \$$$$$$$ |$$ | \$$\\
 $$  ____/  \_______| \_______|\__|  \__|
 $$ |                                $$\\
 $$ |                                $$ |
 $$$$$$$\   $$$$$$\   $$$$$$\   $$$$$$$ |
 $$  __$$\ $$  __$$\ $$  __$$\ $$  __$$ |
 $$ |  $$ |$$ /  $$ |$$ /  $$ |$$ /  $$ |
 $$ |  $$ |$$ |  $$ |$$ |  $$ |$$ |  $$ |
 $$ |  $$ |\$$$$$$  |\$$$$$$  |\$$$$$$$ |
 \__|  \__| \______/  \______/  \_______|

"""

    b = """

  ██▓███  ▓█████ ▄▄▄       ██ ▄█▀ ██░ ██  ▒█████░          ▒█████▄
 ▓██░  ██▒▓█   ▀▒████▄     ██▄█▒ ▓██░ ██▒▒██▒  ██▒ ▒██████▒▒██▀ ██▌
 ▓██░ ██▓▒▒███  ▒██  ▀█▄  ▓███▄░ ▒██▀▀██░▒██░  ██▒▒██▒  ██▒░██   █▌
 ▒██▄█▓▒ ▒▒▓█  ▄░██▄▄▄▄██ ▓██ █▄ ░▓█ ░██ ▒██   ██▒▒██░  ██▒░▓█▄   ▌
 ▒██▒░   ░░▒████▒▓█   ▓██▒▒██▒ █▄░▓█▒░██▓░ ████▓▒░▒██   ██░░▒████▓
 ▒██▒░          ▒▓█   ▓█▒░       ░▓█▒░█▓▒░          ▒███▓▒░

"""

    c = """

  ██▓ ███▄    █    ▓█████▄  ▄▄▄          ██░ ██  ▒█████           ▒█████▄
 ▓██▒ ██ ▀█   █    ▒██▀ ██▌▒████▄       ▓██░ ██▒▒██▒  ██▒ ▒██████▒▒██▀ ██▌
 ▒██▒▓██  ▀█ ██▒   ░██   █▌▒██  ▀█▄     ▒██▀▀██░▒██░  ██▒▒██▒  ██▒░██   █▌
 ░██░▓██▒  ▐▌██▒   ░▓█▄   ▌░██▄▄▄▄██    ░▓█ ░██ ▒██   ██▒▒██░  ██▒░▓█▄   ▌
 ░██░▒██░   ▓██░   ░▒████▓  ▓█   ▓██▒   ░▓█▒░██▓░ ████▓▒░▒██   ██░░▒████▓
 ░██░                       ▓█   ▓█▒░   ░▓█▒░█▓▒░          ▒███▓▒░

"""

    banner.append(a)
    # banner.append(a)
    # banner.append(a)
    # banner.append(a)
    # banner.append(b)
    return(random.choice(banner))


################################################################################

def extract_multicore_wrapper(extract_out_folder, extract_cmd, dataset_id):
    output = subprocess.getoutput(extract_cmd)

    # Save output.
    run_log_file = extract_out_folder + "/run.peakhood_extract.log"
    RUNLOGOUT = open(run_log_file, "w")
    RUNLOGOUT.write(output)
    RUNLOGOUT.close()

    # Check for errors in log file.
    error_found = check_string_in_file(run_log_file, "AssertionError")
    if error_found:
        print(output)
    assert not error_found, "An assertion error was raised during this peakhood extract run, check run log file %s for details" % (
        run_log_file)

    # Check for results.
    stats_out_file = extract_out_folder + "/extract_stats.out"
    assert os.path.exists(
        stats_out_file), "missing extract_stats.out file inside folder %s. Probably this peakhood extract run produced errors, check run log file %s" % (
    extract_out_folder, run_log_file)
    extr_stats_dic = read_settings_into_dic(stats_out_file,
                                                    check=False)
    assert extr_stats_dic, "no stats extracted from extract_stats.out file inside folder %s. Probably this peakhood extract run produced errors, check run log file %s" % (
    extract_out_folder, run_log_file)
    c_intergen_sites = int(extr_stats_dic["c_intergen_sites"])
    c_intron_sites = int(extr_stats_dic["c_intron_sites"])
    c_exon_sites_tc = int(extr_stats_dic["c_exon_sites_tc"])
    c_exon_sites_merged_tc = int(extr_stats_dic["c_exon_sites_merged_tc"])
    c_exon_sites_gc = int(extr_stats_dic["c_exon_sites_gc"])
    c_all_sites = c_intergen_sites + c_intron_sites + c_exon_sites_gc + c_exon_sites_tc
    c_exonic_sites = int(extr_stats_dic["c_exonic_sites"])
    c_exb_sites = int(extr_stats_dic["c_exb_sites"])

    # Percentage of exonic sites.
    perc_exonic_sites = "0.0 %"
    if c_exonic_sites and c_all_sites:
        perc_exonic_sites = "%.2f " % (
                    (c_exonic_sites / c_all_sites) * 100) + "%"
    # Percentage of spliced context sites.
    perc_exonic_tc_sites = "0.0 %"
    if c_exon_sites_tc and c_exonic_sites:
        perc_exonic_tc_sites = "%.2f " % (
                    (c_exon_sites_tc / c_exonic_sites) * 100) + "%"
    # Percentage of exon border sites (connected by ISR reads).
    perc_exb_sites = "0.0 %"
    if c_exb_sites and c_exon_sites_tc:
        perc_exb_sites = "%.2f " % (
                    (c_exb_sites / c_exon_sites_tc) * 100) + "%"

    dataset_print = ""
    dataset_print += "dataset:  %s\n" % (dataset_id)
    dataset_print += "# of all sites                                    %i\n" % (c_all_sites)
    dataset_print += "# of intronic sites:                              %i\n" %(c_intron_sites)
    dataset_print += "# of intergenic sites:                            %i\n" %(c_intergen_sites)
    dataset_print += "# of exonic sites (assigned genome context):      %i\n" %(c_exon_sites_gc)
    dataset_print += "# of exonic sites (assigned transcript context):  %i\n" %(c_exon_sites_tc)
    dataset_print += "# of sites after merging exon border sites:       %i\n" %(c_exon_sites_merged_tc)
    dataset_print += "Percentage (# exonic sites / # all input sites):\n%s\n" %(perc_exonic_sites)
    dataset_print += "Percentage (# transcript context sites / # exonic sites):\n%s\n" %(perc_exonic_tc_sites)
    dataset_print += "Percentage (# transcript context sites / # exonic sites):\n%s\n" %(perc_exonic_tc_sites)
    dataset_print += "Percentage (# exon border sites / # transcript context sites):\n%s\n\n" %(perc_exb_sites)
    return dataset_print
