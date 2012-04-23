########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2012 J.Craig Venter Institute.
########################################################################################



#!/usr/bin/python
#################################################
## 	A new program
#################################################
import sys
from optparse import OptionParser
from collections import defaultdict
_author="Sebastian Szpakowski"
_date="2011/01/01"
_version="Version 1"

#################################################
##		Classes
##

#################################################
##		Functions
##
	################################################
	### Read in a file and return a list of lines
	###
def	loadLines(x):
	try:
		fp = open(x, "r")
		cont=fp.readlines()
		fp.close()
		#print "%s line(s) loaded."  % (len(cont))
	except:
		cont=""
		#print "%s cannot be opened, does it exist? " % ( x )	
	return cont
#################################################
##		Arguments
##
parser = OptionParser()

#parser.add_option("-f", "--file", dest="filename",
#                  help="write report to FILE", metavar="FILE")
#parser.add_option("-q", "--quiet",
#                  action="store_false", dest="verbose", default=True,
#                  help="don't print status messages to stdout")


(options, args) = parser.parse_args()

#################################################
##		Begin
##

mappingfile = loadLines(sys.argv[1])
column = sys.argv[2]
headers = defaultdict(str)
mapping = defaultdict()

counter=0
for header in mappingfile[0].strip().split(",") :
	headers[header] = counter
	counter+=1
	
offset = headers[column]
if offset == "":
	print "requested column header not found (%s)" % (column)
	print headers
	sys.exit(1)
groups = loadLines(sys.argv[3])

for maps in mappingfile[1:]:
	maps = maps.strip().split(",")
	mapping[maps[0]] = maps[offset].replace(" ", "_")


otpt = open("%s.%s.group" % (sys.argv[1], column), "w")
for line in groups:
	line = line.strip().split("\t")
	x = mapping[line[1]]
	if x == "":
		x = "other"
	o = "%s\t%s\n" % (line[0], x)
	otpt.write(o)
otpt.close()	
	


#################################################
##		Finish
#################################################
