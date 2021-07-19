#!/usr/bin/env python3


import argparse
import os


################################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Given a --in BED file (6-column), generate unique column 4 IDs, 
    output to STDOUT.

    """

    p = argparse.ArgumentParser(add_help=False,
                                prog="bed_generate_unique_ids.py",
                                description=help_description,
                                formatter_class=argparse.MetavarTypeHelpFormatter)

    p.add_argument("-h", "--help",
                   action="help",
                   help="Print help message")
    p.add_argument("--in",
                   dest="in_bed",
                   type=str,
                   metavar='str',
                   required = True,
                   help = "Input BED file (6-column BED) to generate new col-4 IDs for")
    p.add_argument("--id",
                   dest="stem_id",
                   type=str,
                   metavar='str',
                   required = True,
                   help = "Part of ID string common to all IDs. E.g. --id CLIP will generate IDs CLIP_1, CLIP_2, CLIP_3 etc. ")
    return p


################################################################################

if __name__ == '__main__':

    parser = setup_argument_parser()
    args = parser.parse_args()

    assert os.path.exists(args.in_bed), "--in %s BED file not found" %(args.in_bed)

    c_out = 0
    with open(args.in_bed) as f:
        for line in f:
            row = line.strip()
            cols = line.strip().split("\t")
            chr_id = cols[0]
            site_s = cols[1]
            site_e = cols[2]
            site_sc = cols[4]
            site_pol = cols[5]

            c_out += 1
            new_id = args.stem_id + "_" + str(c_out)
            print("%s\t%s\t%s\t%s\t%s\t%s" %(chr_id, site_s, site_e, new_id, site_sc, site_pol))
    f.closed


