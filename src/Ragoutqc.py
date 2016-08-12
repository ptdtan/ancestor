#!/usr/bin/env python

import sys
import os
import networkx as nx
import subprocess
import commands
from copy import copy
import logging

def exeCommand(sCommand):
    ###Get all output data
    #print sCommand
    process = subprocess.Popen(sCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        close_fds=True)
    out, errData = process.communicate()

    ###Get all response data
    ###If there is error
    if ((errData != None) and (len(errData) > 0)):
        sys.stderr.write("Command has error:{0}".format(errData))
	if out:
		sys.stdout.write(out)
    pass

def run_command(chr, wall, pre, next, bam, info = None):
	start, end = wall-500, wall+500
	command = ' '.join(['sambamba', 'view', '-t', '32', '-F',  'proper_pair',  bam,  "".join([chr, ':',
						str(start),'-', str(end)]), '|', './check_span.py', str(wall), pre, next, chr, " "]) + " ".join(info)
	#print command
	exeCommand(command)

def _pair_iter(adjacencies):
	for pre, next in zip(adjacencies[:-1], adjacencies[1:]):
		yield pre, next

def read_links(ifile):
    """Parser for scaffolds_links file
    @param ifile: _scaffolds.links file
    @return
    """
	ofile = open(ifile)
	scaffolds_links = ofile.read().strip().split('\n\n')
	chr_scaffolds = {}
	for block in scaffolds_links:
		arr = block.split('\n')
		chr_scaffolds[arr[0]] = arr[1:]
	return chr_scaffolds

def generate_path(chr_scaffolds):
	chr_paths = {}
	for chr, contigs in chr_scaffolds.items():
		chr_name = chr
		chr_paths[chr_name] = [scf.strip().split() for scf in contigs ]
	return chr_paths

def verify_read_span(file, bam):
	chr_scaffolds = read_links(file)
	for chr, contigs in chr_scaffolds.items():
		chr_name = chr
		for pre, next in _pair_iter(contigs):
			try:
				preName, preStart, preLen, preGaps,__ = pre.strip().split()
				nextName, nextStart, nextLen, nextGaps,__ = next.strip().split()
				if int(preGaps) > 50:
					continue
				wall = int(nextStart)
				run_command(chr_name, wall, preName, nextName, bam,
								info=[preStart, preLen, nextLen, preGaps)
			except KeyboardInterrupt:
				sys.exit(0)

def reduceN50(LINK_STREAM, N50_STREAM, BAM_FILE)
	chr_scaffolds = read_links(sys.argv[1])
	ologger = LogInstance(LINK_STREAM+".log")
	fout = open(LINK_STREAM+".red", "w")
	fN50 = open(N50_STREAM+".n50", "w")
	bam = BAM_FILE
	chr_paths = generate_path(chr_scaffolds)
	chr_new_paths = dict()
	for chr, scfs in chr_paths.items():
		ologger.info("Infering %s" %chr)
		new_paths = []
		temp_path = [scfs[0]]
		for pre, next in _pair_iter(scfs):
            ###Add bamtools 2.4.0 to LD_LIBRARY_PATH first
			status, ret = commands.getstatusoutput("./verify %s %s %s %s %s" %(bam, chr, pre[1], next[1], next[2]))
			if ret == "True":
				temp_path.append(next)
				#print "./verify %s %s %s %s %s" %(bam, chr, pre[1], next[1], next[2])
				#print pre[0], next[0]
			else:
				new_paths.append(copy(temp_path))
				temp_path = [next]
				#print "./verify %s %s %s %s %s" %(bam, chr, pre[1], next[1], next[2])
                #print pre[0], next[0]
				#print ret
		new_paths.append(copy(temp_path))
		chr_new_paths[chr] = copy(new_paths)
		ologger.info([len(path) for path in new_paths])
		for path in new_paths:
			fN50.write("%d\n" %(sum([int(scf[2])+int(scf[3]) for scf in path])))
			#print "%d" %(sum([int(scf[2])+int(scf[3]) for scf in path]))

	for chr, paths in chr_new_paths.items():
		for idx, path in enumerate(paths):
			fout.write("%s.%s\n" %(chr, idx))
			fout.write('\n'.join(['\t'.join(scf) for scf in path])+'\n')
