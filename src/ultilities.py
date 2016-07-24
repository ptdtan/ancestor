#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module provides some basic FASTA I/O
"""

import logging

from string import maketrans
from Bio import SeqIO

logger = logging.getLogger()

class FastaError(Exception):
    pass

def write_fasta_dict(fasta_dict, filename):
    """
    Writes dictionary with fasta to file
    """
    with open(filename, "w") as f:
        for header in sorted(fasta_dict):
            f.write(">{0}\n".format(header))

            for i in range(0, len(fasta_dict[header]), 60):
                f.write(fasta_dict[header][i:i + 60] + "\n")


COMPL = maketrans("ATGCURYKMSWBVDHNXatgcurykmswbvdhnx",
                  "TACGAYRMKSWVBHDNXtacgayrmkswvbhdnx")
def reverse_complement(string):
    return string[::-1].translate(COMPL)

def _validate_seq(sequence):
    VALID_CHARS = "ACGTURYKMSWBDHVNXatgcurykmswbvdhnx"
    if len(sequence.translate(None, VALID_CHARS)):
        return False
    return True

def _calc_n50(scaffolds_lengths, assembly_len):
    n50 = 0
    sum_len = 0
    for l in sorted(scaffolds_lengths, reverse=True):
        sum_len += l
        if sum_len > assembly_len / 2:
            n50 = l
            break
    return n50
