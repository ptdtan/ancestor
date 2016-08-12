#!/usr/bin/env python

import sys
import os

def check_spaning(strwall, pre, next, chr, info1, info2, info3, info4):
	"""
	check if any read pair that span 2 sides of the wall
	"""
	wall = int(strwall)
	strMapinfo = sys.stdin.readline()
	count=0
	while strMapinfo:
		arrMaprecord = strMapinfo.strip().split("\t")
		start, end = int(arrMaprecord[3]), int(arrMaprecord[7])
		if (start-wall)*(end-wall) < 0 :
			count+=1
			#print pre, next, True
		strMapinfo = sys.stdin.readline()
	"""if not count:
		print pre, next, False, count"""
	if count > 100:
		print pre, next, info2, info3
		os.system(" ".join(["./verify_ghaGan1", chr, info1, info2, info3, info4])) #change the verification command to make consistent with your Genome
wall, pre, next, chr, info1, info2, info3, info4 = sys.argv[1:]
check_spaning(wall, pre, next, chr, info1, info2, info3, info4)
