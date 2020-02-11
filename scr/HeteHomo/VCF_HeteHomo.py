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


with open(args.output, "w") as fout:
    if args.input.endswith(".gz"):
      fin = gzip.open(args.input, 'rb')
    else:
      fin = open(args.input, "r")
    try:
      while True:
        line = fin.readline()
        if line.startswith("Chr"):
          fout.write(line.rstrip())
          vcfheaders = line.rstrip().split("\t")
          format_index = vcfheaders.index("FORMAT")
          sample_index = format_index + 1
          
          break
        else:
          fout.write(line)
      fout.write("\t"+'0'+"\t"+'1'+"\t"+'2'+"\t"+"unknown"+"\n")
      totalsnv = 0
      for line in fin:
        table={'0/0': 0, '0/1': 0, '1/1': 0, '1/0': 0,'./.': 0}
        snv = line.rstrip().split('\t')
        for si in range(sample_index, len(vcfheaders)):
          sampleData = snv[si]
          heteHomo=sampleData.split(":")[0]
          if heteHomo in table:
            table[heteHomo] += 1
          else:
            table[heteHomo] = 1
        fout.write("%s\t%s\t%s\t%s\t%s\n" %(line.rstrip(),table['0/0'],table['0/1'],table['1/1'],table['1/0']+table['./.']))

        #  else:
        #    fout.write(line)
        #else:
        #  fout.write(line)
          
      logger.info("Done")
    finally:
      fin.close()