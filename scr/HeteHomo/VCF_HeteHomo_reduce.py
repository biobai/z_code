#!/usr/bin/env python
# encoding: utf-8
"""
VCF_HeteHomo.py
"""
import argparse
import sys
import logging
import os
import csv
import gzip
import re
from collections import defaultdict

def initialize_logger(logfile, logname, isDebug):
  logger = logging.getLogger(logname)
  loglevel = logging.DEBUG if isDebug else logging.INFO
  logger.setLevel(loglevel)

  formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')    
 
  # create console handler and set level to info
  handler = logging.StreamHandler()
  handler.setLevel(loglevel)
  handler.setFormatter(formatter)
  logger.addHandler(handler)
 
  # create error file handler and set level to error
  handler = logging.FileHandler(logfile, "w")
  handler.setLevel(loglevel)
  handler.setFormatter(formatter)
  logger.addHandler(handler)
 
  return(logger)

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="filter SV that produced by smoove",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input VCF file', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output file name", required=NotDEBUG)
parser.add_argument('--debug', action='store_true', help="Output debug information", default=False)

args = parser.parse_args()

if DEBUG:
  args.input = "./IPF_batch2.smoove.square.anno.vcf.gz"
  args.output = "./IPF_batch2.smoove.square.anno.vcf.Filtered.vcf"


logger = initialize_logger(args.output + ".log", 'SV_filter', args.debug)
logger.info(str(args))

sample_name={}
with open(args.output, "w") as fout:
    if args.input.endswith(".gz"):
      fin = gzip.open(args.input, 'rb')
    else:
      fin = open(args.input, "r")
    try:
      while True:
        line = fin.readline()
        if line.startswith("Chr"):
          #fout.write(line.rstrip())
          vcfheaders = line.rstrip().split("\t")
          format_index = vcfheaders.index("FORMAT")
          sample_index = format_index + 1
	  fout.write('\t'.join(line.rstrip().split('\t')[:format_index])+"\tPositive_Samples")
          for vi in range(sample_index, len(vcfheaders)):
	    sample_name[vi]=vcfheaders[vi]
          break
        else:
          fout.write(line)
      fout.write("\t"+'0'+"\t"+'1'+"\t"+'2'+"\t"+"unknown"+"\n")
      totalsnv = 0
      for line in fin:
        table=defaultdict(int)
        name=defaultdict(list)
	fname=defaultdict()
	snv = line.rstrip().split('\t')
	#header=snv[:format_index-1]
        for si in range(sample_index, len(vcfheaders)):
          sampleData = snv[si]
          heteHomo=sampleData.split(":")[0]
          table[heteHomo] += 1
          name[heteHomo].append(sample_name[si])
        for tkey in table:
	    fname[tkey]=','.join(name[tkey])
	middle=' -- '.join('{}:{}'.format(key, value) for key, value in sorted(fname.items()) if key != './.' and key != '0/0' ) 
	#for key, value in fname.items():
	#    if key != './.' and key != '0/0':
	#	middle=' -- '.join('%s:%s' % (key, value))
	fout.write("%s\t%s\t%s\t%s\t%s\t%s\n" %('\t'.join(snv[:format_index]),str(middle),str(table['0/0']),str(table['0/1']),str(table['1/1']),str(table['1/0']+table['./.'])))

        #  else:
        #    fout.write(line)
        #else:
        #  fout.write(line)
          
      logger.info("Done")
    finally:
      fin.close()
