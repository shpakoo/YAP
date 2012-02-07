########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2012 Sebastian Szpakowski, J.Craig Venter Institute.
########################################################################################



#!/usr/bin/python
#################################################
## 	A new program
#################################################
import sys, tempfile, shlex, glob, os
from subprocess import *
from collections import defaultdict

_author="Sebastian Szpakowski"
_date="2011/09/14"
_version="Version 1"

#################################################
##		Classes
##


class	MothurCommand():
	def	__init__(self):
		self.name = ""
		self.input_files = dict()
		self.output_files = dict()
		self.args = dict()			
		self.annotation = dict()		
		self.text = ""
		self.linkbacks = defaultdict(list)
				
	def	addArg(self, arg, val):
		self.text = "%s\n%s\t%s" % (self.text, arg[0:50].ljust(50), val[0:50].ljust(50)) 
		if arg == "commandName":
			self.name = val
		else:
			self.determineType(arg, val)
		
	def	determineType(self, arg, val):
		vals = val.split("|")
		
		if arg in ["citation", "description", "commandCategory", "help"]:
			self.annotation[arg]=(val)
			
		elif arg == "outputTypes":
			for o in val.split("-"):
				self.output_files[o]= ""
				
		elif arg == "inputTypes":	
			for i in val.split("-"):
				self.input_files[i]= ()
					
		elif len(vals)==3:
			vals[2] = self.booleanify(vals[2])
			vals[0] = vals[0].strip().split("-")
			self.args[arg] = vals
			
		elif arg in self.input_files.keys():
			req, one, oneplus, linked = vals
			self.input_files[arg] = (req, one, oneplus, linked)	
			
			for K in one, oneplus, linked:
				if K !="none":
					for L in K.strip().split("-"):
						self.linkbacks[L].append(arg) 
							
		else:	 
			self.args[arg] = vals
			#print self.name, ":", arg, "=", val
			
			
	def	booleanify(self, x):
			if x in ["T", "t", "True", "true", "TRUE"]:
				return True
			elif x in ["F", "f", "False", "false", "FALSE"]:
				return False	
			else:
				raise ValueError
	
	def	getInputs(self):
		return self.input_files.keys()
		
	def	getArguments(self):
		return self.args.keys()
		
	def	isAnArgument(self, x):
		return (x in self.args.keys())
	
	def	getOutputs(self):
		return self.output_files.keys()
		
	def	isRequired(self, filetype):
		return (self.booleanify(self.input_files[filetype][0]))
		
	def	getAlternativeInput(self, filetype):
		x = set(self.linkbacks[self.input_files[filetype][1]])
		x = x.difference(set([filetype]))
		return (x)
		
	def	getOtherInput(self, filetype):
		x = set(self.linkbacks[self.input_files[filetype][2]])
		x = x.difference(set([filetype]))
		return (x)
		
	def	getLinkedInput(self, filetype):
		x = set(self.linkbacks[self.input_files[filetype][3]])
		x = x.difference(set([filetype]))
		return (x)
		
		
	def	getDefault(self, arg):
		return (self.args[arg])	
				
	def	__str__(self):
		otpt = "%s:" % (self.name)
		otpt = "%s\n%s\t%s" % (otpt, "input".rjust(10)  ,"".join([x.ljust(20) for x in self.input_files.keys()]))
		otpt = "%s\n%s\t%s" % (otpt, "req".rjust(10)  ,"".join(["".join("%s" % self.isRequired(x))[:18].ljust(20) for x in self.input_files.keys()]))
		otpt = "%s\n%s\t%s" % (otpt, "alt".rjust(10)  ,"".join(["-".join(self.getAlternativeInput(x))[:18].ljust(20) for x in self.input_files.keys()]))
		otpt = "%s\n%s\t%s" % (otpt, "multi".rjust(10)  ,"".join(["-".join(self.getOtherInput(x))[:18].ljust(20) for x in self.input_files.keys()]))
		otpt = "%s\n%s\t%s" % (otpt, "link".rjust(10)  ,"".join(["-".join(self.getLinkedInput(x))[:18].ljust(20) for x in self.input_files.keys()]))
		
		otpt = "%s\n%s\t%s" % (otpt, "output".rjust(10)  ,"".join([x.ljust(20) for x in self.output_files.keys()]))
		otpt = "%s\n%s\t%s" % (otpt, "---".rjust(10)  ,"".join([self.output_files[x].ljust(20) for x in self.output_files.keys()]))
		otpt = "%s\n%s" % (otpt, self.text)	
		return otpt
		
		
class	MothurCommandInfo():
	def	__init__(self,path=""):
		tmpdir = tempfile.mkdtemp(prefix="output")
		p    = Popen(shlex.split("%smothur #get.commandinfo(output=mothurinfo)"% (path)), stdout=PIPE, stderr=PIPE, cwd=tmpdir)
		log  = p.stdout.readlines()
		info = loadLines("%s/mothurinfo" % (tmpdir))	
		
		self.version, self.date = self.parseLog(log)
		
		self.commands = dict()	
		self.parseInfo(info)
		
		self.cleanUpDir(tmpdir)	
				
				
				
	def parseLog(self, log):
		version = ""
		date = ""
		for line in log:
			if line.startswith("mothur v"):
				version = line[9:].strip()
			if line.startswith("Last updated"):
				date = line[14:].strip()
		
		print ("Found Mothur: %s [%s]" % (version, date))
		return (version, date)	
		
	def parseInfo(self, info):
		current = MothurCommand()
		counter = 0 
		for line in info:
			line = line.strip().split("=")
			
			if len(line)>1:
				if line[0]=="commandName":
					self.addCommand(current)
					current = MothurCommand()					
				current.addArg(line[0], "=".join(line[1:]))					
			else:
				counter = int(line[0])
				
			self.addCommand(current)	
			
		print 	"Identified %s commands [reported: %s]" % (len(self.commands.keys()), counter)
			
	def	addCommand(self, K):
		if K.name!="":
			self.commands[K.name] = K	
																	
	def cleanUpDir(self, path, counter=1):	
		try:	
			for file in glob.glob("%s/*" % (path)):
				os.remove(file)
			os.rmdir(path)
		except OSError:
			print "Tmp cleanup delayed... x%s [%s]" % (counter, path)
			time.sleep(1)
			self.cleanUpDir(path, counter=counter+1)
			
	def	getCommandInfo(self, x):
		return self.commands[x]
		

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
	
def test():
	M = MothurCommandInfo()
	x = M.getCommandInfo("dist.seqs")
	print x
	
	x = M.getCommandInfo("unique.seqs")
	print x
	
	x = M.getCommandInfo("cluster")
	print x
	
	x = M.getCommandInfo("trim.seqs")
	print x
	
	x = M.getCommandInfo("align.seqs")
	print x
	
	x = M.getCommandInfo("remove.seqs")
	print x

#################################################
##		Arguments
##

#################################################
##		Begin
##

#test()
#################################################
##		Finish
#################################################
