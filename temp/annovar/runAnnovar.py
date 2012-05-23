#!/usr/bin/env python

########################################################################
##
## Author: Ryan
## Version: 0.1a
## This program execute annovar pacakge with common usages
##
#########################################################################

import os, sys, math, re
from optparse import OptionParser


parser = OptionParser(usage="usage: python %prog [options]")
parser.add_option( "-i", action="store", dest="inputVCF", help="VCF input files", metavar="FILE" )
parser.add_option( "-t", action="store", dest="annotationType", help="annotationType : gene, region, filter (default Gene-based)", type="string", default="gene" )
parser.add_option( "-d", action="store", dest="dbDir", help="Directory that the database are written to", type="string" )
parser.add_option( "-b", action="store", dest="dbType", help="Datebase type", type="string" )
parser.add_option( "-r", action="store", dest="reference", help="Reference i.e., hg19, hg18, (default hg19)", type="string", default="hg19" )
parser.add_option( "-g", action="store", dest="gffdb", help="dff3db file, it is required when dbType is gff3", type="string")

(options,args)=parser.parse_args(sys.argv)
annovarDir="/ifs/data/shares/deepsequencing/solid/src/annovar"



cmd="perl %s/annotate_variation.pl -downdb %s -buildver %s %s" %(annovarDir, options.dbType, options.reference, options.dbDir)
os.system(cmd)
print "Download the requested databases"

"""gene-based annotation of variants in the varlist file"""
if options.annotationType=="gene":
  cmd="perl %s/annotate_variation.pl -%sanno --buildver %s %s %s" %(annovarDir, options.annotationType, options.reference, options.inputVCF, options.dbDir)
  os.system(cmd)
  print "Annotate Gene-based VCF file"

if options.dbType=="":
  parser.error("Please choose dbType -b")
else:
  """region-based annotate variants"""
  if options.annotationType=="region":
    if options.dbType=="gff3":
      if options.gffdb=="":
        parser.error("Please choose dff3db -g")
      else:
        cmd="perl %s/annotate_variation.pl -%sanno --buildver %s -dbtype %s -gff3dbfile %s %s %s" %(annovarDir, options.annotationType, options.reference, options.dbType, options.gffdb, options.inputVCF, options.dbDir)
        os.system(cmd)
        print "Annotate Region-based VCF file"
    else:
      cmd="perl %s/annotate_variation.pl -%sanno --buildver %s -dbtype %s %s %s" %(annovarDir, options.annotationType, options.reference, options.dbType, options.inputVCF, options.dbDir)
      os.system(cmd)
      print "Annotate Region-based VCF file"

  """filter rare or unreported variants (in 1000G/dbSNP) or predicted deleterious variants"""
  if options.annotationType=="filter":
    cmd="perl %s/annotate_variation.pl -%s --buildver %s -dbtype %s %s %s" %(annovarDir, options.annotationType, options.reference, options.dbType, options.inputVCF, options.dbDir)
    os.system(cmd)
    print "Annotate Filter-based VCF file"

