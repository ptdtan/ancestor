#!/usr/bin/env python

import Ragoutqc as qc
import commands
import sys
PCR_INFERER = "./pcr_inferer"

def inferer(chr_name, preName, nextName, bam , info=None):
    cmd = " ".join([PCR_INFERER, bam, chr_name]+info)
    estimated_gaps = commands.getstatusoutput(cmd)
    print preName, nextName, estimated_gaps

def main(LINKS, BAM, meanINS):
    chr_scaffolds = qc.read_links(LINKS)
    for chr, contigs in chr_scaffolds.items():
        for pre, next in qc._pair_iter(contigs):
			try:
				preName, preStart, preLen, preGaps,_p = pre.strip().split()
				nextName, nextStart, nextLen, nextGaps,_n = next.strip().split()
				if int(preGaps) > 50:
					continue
				inferer(chr_name, preName, nextName, bam,
								info=[preStart, nextStart, nextLen, meanINS, preGaps])
			except KeyboardInterrupt:
				sys.exit(0)
if __name__ == "__main__":
    links, bams, mean = sys.argv[1:]
