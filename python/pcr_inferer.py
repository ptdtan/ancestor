#!/usr/bin/env python

import qc as qc
import commands
import sys
import multithreads as mt

PCR_INFERER = "./pcr_inferer"

def inferer(chr_name, preName, nextName, bam , info=None, threshold=50):
    cmd = " ".join([PCR_INFERER, bam, chr_name]+info+[threshold])
    signal, estimated_gaps = commands.getstatusoutput(cmd)
    if estimated_gaps:
        print preName, nextName, estimated_gaps

def main(LINKS, BAM, meanINS, threshold):
    pool = mt.ThreadPool(20)
    chr_scaffolds = qc.read_links(LINKS)
    for chr_name, contigs in chr_scaffolds.items():
        for pre, next in qc._pair_iter(contigs):
			try:
				preName, preStart, preLen, preGaps = pre.strip().split()[:4]
				nextName, nextStart, nextLen = next.strip().split()[:3]
				if int(preGaps) > 50:
					continue
				pool.add_task(inferer, chr_name, preName, nextName, bam,
								info=[preStart, nextStart, nextLen, meanINS, preGaps], threshold = threshold)
			except KeyboardInterrupt:
				sys.exit(0)
    pool.wait_completion()
if __name__ == "__main__":
    links, bam, mean, threshold  = sys.argv[1:]
    main(links, bam, mean, threshold )

