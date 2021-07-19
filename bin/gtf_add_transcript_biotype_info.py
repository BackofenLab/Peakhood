#!/usr/bin/env python3

import subprocess
import argparse
import gzip
import sys
import os
import re

"""


GTF example rows:

StringTie:
chr1	StringTie	transcript	11869	14409	.	+	.	transcript_id "ENST00000456328.2"; gene_id "MSTRG.47"; gene_name "DDX11L1"; xloc "XLOC_000047"; ref_gene_id "ENSG00000223972.5"; cmp_ref "ENST00000456328.2"; class_code "="; tss_id "TSS47";

Ensembl GRCh38.103:
1	havana	transcript	11869	14409	.	+	.	gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; tag "basic"; transcript_support_level "1";

gencode.v37:
chr1	HAVANA	transcript	11869	14409	.	+	.	gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";


"""

################################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Extract transcript biotype info from --in-gtf,
    and add it to --add-gtf file. Write new GTF to --out-gtf.

    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="gtf_add_transcript_biotype_info.py",
                                description=help_description,
                                formatter_class=argparse.MetavarTypeHelpFormatter)

    # Argument groups.
    p_man = p.add_argument_group("REQUIRED ARGUMENTS")
    p_opt = p.add_argument_group("OPTIONAL ARGUMENTS")
    
    # Arguments.
    p_opt.add_argument("-h", "--help",
                   action="help",
                   help="Print help message")
    p_man.add_argument("--in-gtf",
                   dest="in_gtf",
                   type=str,
                   required = True,
                   help = "GTF file to take gene and transcript infos from")
    p_man.add_argument("--add-gtf",
                   dest="add_gtf",
                   type=str,
                   required = True,
                   help = "GTF file to add gene and transcript infos to")
    p_man.add_argument("--out-gtf",
                   dest="out_gtf",
                   type=str,
                   required = True,
                   help = "Output GTF file to write -add-gtf with added infos to")
    p_opt.add_argument("--only-std-chr",
                   dest="only_std_chr",
                   default = False,
                   action = "store_true",
                   help = "Output only standard chromosome entries (chr1, chr2 ... ) (default: False)")
    return p



################################################################################

def gtf_get_transcript_biotypes(in_gtf):
    """
    Get transcript biotype info for each transript ID.
    Return transcript ID -> transcript biotype mapping.

    """

    # Biotype to count dic.
    trid2tbt_dic = {}
 
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
        feature = cols[2]
        infos = cols[8]
        if feature != "transcript":
            continue

        # Extract transcript ID.
        m = re.search('transcript_id "(.+?)"', infos)
        assert m, "transcript_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        transcript_id = m.group(1)

        # Extract transcript biotype.
        m = re.search('transcript_biotype "(.+?)"', infos)
        if not m:
            m = re.search('transcript_type "(.+?)"', infos)
        assert m, "transcript_biotype or transcript_type entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        transcript_biotype = m.group(1)
        
        trid2tbt_dic[transcript_id] = transcript_biotype

    f.close()
    # Check and return to shack.
    assert trid2tbt_dic, "no transcript biotype information read in"
    return trid2tbt_dic


################################################################################

def check_convert_chr_id(chr_id):
    """
    Check and convert chromosome IDs to format:
    chr1, chr2, chrX, ...
    If chromosome IDs like 1,2,X, .. given, convert to chr1, chr2, chrX ..
    Return False if given chr_id not standard and not convertable.

    Filter out scaffold IDs like:
    GL000009.2, KI270442.1, chr14_GL000009v2_random
    chrUn_KI270442v1 ...

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

    """
    assert chr_id, "given chr_id empty"

    if re.search("^chr", chr_id):
        if not re.search("^chr[\dMXY]+$", chr_id):
            chr_id = False
    else:
        # Convert to "chr" IDs.
        if chr_id == "MT":
            chr_id = "M"
        if re.search("^[\dMXY]+$", chr_id):
            chr_id = "chr" + chr_id
        else:
            chr_id = False
    return chr_id


################################################################################

def gtf_add_infos(add_gtf, out_gtf, trid2tbt_dic,
                  only_std_chr=False):
    """
    Add infos to add_gtf. Output to out_gtf.

    """


    GTFOUT = open(out_gtf, "w")

    c_added_info = 0
    c_no_new_info = 0

    if re.search(".+\.gz$", add_gtf):
        f = gzip.open(add_gtf, 'rt')
    else:
        f = open(add_gtf, "r")
    for line in f:

        if re.search("^#", line):
            GTFOUT.write(line)
            continue

        cols = line.strip().split("\t")
        col1 = cols[0]
        col2 = cols[1]
        feature = cols[2]
        col4 = cols[3]
        col5 = cols[4]
        col6 = cols[5]
        col7 = cols[6]
        col8 = cols[7]
        infos = cols[8]

        if only_std_chr:
            new_chr_id = check_convert_chr_id(col1)
            if not new_chr_id:
                continue
            else:
                col1 = new_chr_id

        if feature == "transcript":

            m = re.search('transcript_id "(.+?)"', infos)
            assert m, "transcript_id entry missing in GTF file \"%s\", line \"%s\"" %(add_gtf, line)
            tr_id = m.group(1)
        
            tr_biotype_str = 'transcript_type "-";'
            if tr_id in trid2tbt_dic:
                c_added_info += 1
                tr_biotype = trid2tbt_dic[tr_id]
                tr_biotype_str = 'transcript_type "' + tr_biotype + '";'
            else:
                c_no_new_info += 1


            new_infos = infos + " " + tr_biotype_str

            GTFOUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(col1, col2, feature, col4, col5, col6, col7, col8, new_infos))

        else:
            GTFOUT.write(line)

    f.close()
    GTFOUT.close()

    print("# --add-gtf transcript IDs with added info from --in-gtf:  %s" %(c_added_info))
    print("# --add-gtf transcript IDs not appearing in --in-gtf:      %s" %(c_no_new_info))



################################################################################


if __name__ == '__main__':

    parser = setup_argument_parser()
    args = parser.parse_args()

    assert os.path.exists(args.in_gtf), "--in-gtf %s not found" %(args.in_gtf)
    assert os.path.exists(args.add_gtf), "--add-gtf %s not found" %(args.add_gtf)
    assert args.in_gtf != args.add_gtf, "--in-gtf == --add-gtf (!)"

    print("Get transcript biotype info from --in-gtf ... ")
    trid2tbt_dic = gtf_get_transcript_biotypes(args.in_gtf)
    
    c_trids = len(trid2tbt_dic)
    tbt2c_dic = {}
    for tr_id in trid2tbt_dic:
        tbt = trid2tbt_dic[tr_id]
        if tbt not in tbt2c_dic:
            tbt2c_dic[tbt] = 1
        else:
            tbt2c_dic[tbt] += 1
    c_tbts = len(tbt2c_dic)

    print("# of transcript IDs read in from --in-gtf:  %i" %(c_trids))
    print("# of associated transcript biotypes:        %i" %(c_tbts))
    print("")

    if tbt2c_dic:
        print("Encountered transcript biotypes and counts:")
        for tbt, tbt_c in sorted(tbt2c_dic.items(), key=lambda item: item[1], reverse=True):
            print("\"%s\" %i" %(tbt, tbt_c))
    print("")

    print("Add info to --add-gtf and store in new --out-gtf ...") 
    gtf_add_infos(args.add_gtf, args.out_gtf, trid2tbt_dic,
                  only_std_chr=args.only_std_chr)
    print("")

