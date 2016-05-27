#!/usr/bin/env python

import sys
import os
import networkx as nx
import subprocess
import commands
from copy import copy

import logging

def LogInstance(log_file):
    """
    Turns on logging, sets debug levels and assigns a log file
    """
    logger = logging.getLogger()
    log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
                                      "%(message)s", "%H:%M:%S")
    console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: "
                                          "%(message)s", "%H:%M:%S")
    console_log = logging.StreamHandler()
    console_log.setFormatter(console_formatter)
    console_log.setLevel(logging.INFO)

    file_handler = logging.FileHandler(log_file, mode="w")
    file_handler.setFormatter(log_formatter)

    logger.setLevel(logging.DEBUG)
    logger.addHandler(console_log)
    logger.addHandler(file_handler)
    return logger

def exeCommand(sCommand):
    ###Get all output data
    #print sCommand
    process = subprocess.Popen(sCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        close_fds=True)
    out, errData = process.communicate()

    ###Get all response data
    """for lineData in outData.splitlines():
        #if(self.RUNNING_DEBUG_FLAG == 1):
        outStringData = str(lineData)
        print("%s" % (outStringData))"""
    #print "outdata: "+outData
    ###If there is error
    if ((errData != None) and (len(errData) > 0)):
        print("Command has error:{0}".format(errData))
	print out
	return out.splitlines()

def run_command(chr, wall, pre, next, bam, info = None):
	start, end = wall-500, wall+500
	command = ' '.join(['sambamba', 'view', '-t', '32', '-F',  'proper_pair',  bam,  "".join([chr, ':', 
						str(start),'-', str(end)]), '|', './check_span.py', str(wall), pre, next, chr, " "]) + " ".join(info)
	#print command
	os.system(command)

def _pair_iter(adjacencies):
	for pre, next in zip(adjacencies[:-1], adjacencies[1:]):
		yield pre, next

def read_links(ifile):
	ofile = open(ifile)
	scaffolds_links = ofile.read().strip().split('\n\n')
	chr_scaffolds = []
	for block in scaffolds_links:
		arr = block.split('\n')
		chr_scaffolds.append([arr[0]]+arr[1:])
	return chr_scaffolds

def generate_path(chr_scaffolds):
	chr_paths = {}
	for chr in chr_scaffolds:
		chr_name = chr[0]
		chr_paths[chr_name] = [scf.strip().split() for scf in chr[1:]]
	return chr_paths

def main(file, bam):
	chr_scaffolds = read_links(file)
	for chr in chr_scaffolds:
		chr_name = chr[0]
		for pre, next in _pair_iter(chr[1:]):
			try:
				prescf = pre.strip().split()
				nextscf = next.strip().split()
				if int(prescf[3]) > 50:
					continue
				wall = int(nextscf[1])
				run_command(chr_name, wall, prescf[0], nextscf[0], bam, 
								info=[prescf[1], prescf[2], nextscf[2], prescf[3]])
			except KeyboardInterrupt:
				sys.exit(0)

if __name__ == "__main__":
	chr_scaffolds = read_links(sys.argv[1])
	ologger = LogInstance(sys.argv[1]+".log")
	fout = open(sys.argv[1]+".red", "w")
	fN50 = open(sys.argv[1]+".n50", "w")
	bam = sys.argv[2]
	chr_paths = generate_path(chr_scaffolds)
	chr_new_paths = dict()
	for chr, scfs in chr_paths.items():
		ologger.info("Infering %s" %chr)
		new_paths = []
		temp_path = [scfs[0]]
		for pre, next in _pair_iter(scfs):
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
