########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2013 Sebastian Szpakowski
########################################################################################

#################################################
## A library of experimental "steps" or program wrappers to construct pipelines
## Pipeline steps orchestration, grid management and output handling.
#################################################
import sys, tempfile, shlex, glob, os, stat, hashlib, time, datetime, re, curses
from threading import *
from subprocess import *
from MothurCommandInfoWrapper import *
from StepsLibrary import *
from collections import defaultdict
from collections import deque
from random import *
from Queue import *


_author="Sebastian Szpakowski"
_date="2012/09/20"
_version="Version X"

#################################################
##		Classes
##



class GroupSplit(DefaultStep):
	def __init__(self, INS, ARGS, PREV):
		DefaultStep.__init__(self)
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("GroupSplit")
		self.start()
		
	def	performStep(self):
		pattern=[]
		fasta = self.find("fasta")[0]
		groups = self.find("group")[0]
		mapping = defaultdict(str)
		
		fasta = "%s/%s" % (self.stepdir, fasta)
		groups = "%s/%s" % (self.stepdir, groups)
		
		for read, group in GeneralPurposeParser(groups, sep="\t"):
			mapping[read] = group
				
		samples = defaultdict(list)
		for cur in set(mapping.values()):	
			
			fn= "%s/%s.fasta" % (self.stepdir, cur)
			f=open(fn, 'w')
			counter=0
			
			for head, seq in FastaParser(fasta):
				if mapping[head] == cur:
					f.write(">%s\n%s\n" % (head, seq))	
					counter+=1
					
			self.message( "%s\t%s" % (cur, counter))
			f.close()
				
class	TCOFFEE(DefaultStep):
	def __init__(self, INS, ARGS, PREV):
		DefaultStep.__init__(self)
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("TCOFFEE")
		self.start()
						
	def	performStep(self):
		pattern=[]
		files = self.find("fasta")
		tasks = list()
		argstring = ""
		for arg, val in self.arguments.items():
			if arg == "pattern":
				pattern = val.strip().split(",")
			else:	
				argstring = "%s %s %s " % (argstring, arg, val)
			
		template = "t_coffee " 
			
		for f in files:
			if len(pattern)>0:
				for k in pattern:
					if f.find(k)>-1:
						command = "%s %s %s " % (template, f, argstring)
						self.message(command)
						task = GridTask(template="pick", name=self.stepname, command=command, cpu=4,  dependson=list(), cwd = self.stepdir, debug=False)
						task.wait()
			else:
				command = "%s %s %s " % (template, f, argstring)
				self.message(command)
				task = GridTask(template="pick", name=self.stepname, command=command, cpu=4,  dependson=list(), cwd = self.stepdir, debug=False)
				task.wait()	

class	CLUSTALW2(DefaultStep):
	def __init__(self, INS, ARGS, PREV):
		DefaultStep.__init__(self)
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("CLUSTALW2")
		self.start()
						
	def	performStep(self):
		pattern = []
		files = self.find("fasta")
		tasks = list()
		argstring = ""
		for arg, val in self.arguments.items():
			if arg == "pattern":
				pattern = val.strip().split(",")
			else:	
				argstring = "%s %s %s " % (argstring, arg, val)
			
		template = "%sclustalw2 %s " % (binpath, argstring)
			
		for f in files:
			if len(pattern)>0:
				for k in pattern:
					if f.find(k)>-1:
						command = "%s -INFILE=%s" % (template, f)
						self.message(command)
						task = GridTask(template="pick", name=self.stepname, command=command, cpu=1,  dependson=list(), cwd = self.stepdir, debug=False)
						tasks.append(task)
			else:
				command = "%s -INFILE=%s" % (template, f)
				self.message(command)
				task = GridTask(template="pick", name=self.stepname, command=command, cpu=1,  dependson=list(), cwd = self.stepdir, debug=False)
				tasks.append(task)
				
		for task in tasks:
			task.wait()						
							
class	Bowtie1 (DefaultStep):
	def __init__(self, INS, ARGS, PREV):		
		DefaultStep.__init__(self)
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("BOWTIE1aligner")
		self.counter=0
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		
		m1 = self.find("mate1")
		m2 = self.find("mate2")
		m3 = self.find("fasta")
		
		cpus = 1
		
		argstring = ""
		for arg, val in self.arguments.items():
			argstring = "%s %s %s " % (argstring, arg, val) 
			if arg =="-p":
				cpus = val
		
		
		tasks = list()
		self.message("processing %s file(s)..." % (len(m1)+len(m3)))
		for f in m1:
			name = ".".join(f.strip().split(".")[:-1])

			if "%s.mate2" % (name) in m2:
				k = "bowtie %s -1 %s.mate1 -2 %s.mate2 > %s.bowtie1alignment" % (argstring, name, name, name)
		 		#self.message(k)
		 		task = GridTask(template="pick", name=self.stepname, command=k, cpu=cpus, dependson=list(), cwd = self.stepdir, debug=False)	
				tasks.append(task)
			else:
				self.message("skipping: %s" % (name))
				
		for f in m3:
			name = ".".join(f.strip().split(".")[:-1])
			k = "bowtie %s  %s > %s.bowtie1alignment" % (argstring, f, name)
			#self.message(k)
			task = GridTask(template="pick", name=self.stepname, command=k, cpu=cpus, dependson=list(), cwd = self.stepdir, debug=False)	
			tasks.append(task)
								
		for task in tasks:
			task.wait()	

### remove spaces from header (i.e. keep first oken only)
### and make all 50 base fragments (overlapping by 25)	
class	TilingFasta(DefaultStep):
	def __init__(self, INS, PREV):		
		DefaultStep.__init__(self)
		self.setInputs(INS)
		#self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("TilingFasta")
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		f = self.find("fasta")
		for file in f:
			input = FastaParser("%s/%s" % (self.stepdir, file))
			output = open("%s/%s.tile.fasta"% (self.stepdir, file), "w")
			
			for head, seq in input:
				head = head.split()[0]
				counter=1
				while len(seq)>50:
					tmphead = "%s:%s-%s" % (head, counter, counter+100)
					tmpseq = seq[:50]
					
					seq = seq[25:]
					counter+=25

					
					output.write(">%s\n%s\n" % (tmphead, tmpseq))
			
			output.close()

class	Bowtie2 (DefaultStep):
	def __init__(self, INS, ARGS, PREV):		
		DefaultStep.__init__(self)
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("BOWTIE2aligner")
		self.counter=0
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		
		m1 = self.find("mate1")
		m2 = self.find("mate2")
		m3 = self.find("fasta")
		
		cpus = 1
		
		argstring = ""
		for arg, val in self.arguments.items():
			argstring = "%s %s %s " % (argstring, arg, val) 
			if arg =="-p":
				cpus = val
		
		
		tasks = list()
		jobs = len(m1)
		counter=0
		for f in m1:
			counter+=1
			name = ".".join(f.strip().split(".")[:-1])

			if "%s.mate2" % (name) in m2:
				k = "~/bin/bowtie2 %s -q -1 %s.mate1 -2 %s.mate2 -S %s.sam" % (argstring, name, name, name)
		 		if counter==1:
		 			self.message(k)
		 		elif counter==2:
		 			self.message("processing %s file(s)..." % (jobs))
		 			
		 		task = GridTask(template="pick", name=self.stepname, command=k, cpu=cpus, dependson=list(), cwd = self.stepdir, debug=False)	
				tasks.append(task)
			else:
				self.message("skipping: %s" % (name))
		
		jobs = len(m3)	
		counter=0		
		for f in m3:
			counter+=1
			name = ".".join(f.strip().split(".")[:-1])
			k = "~/bin/bowtie2 %s -f -U %s -S %s.sam" % (argstring, f, name)
			if counter==1:
		 		self.message(k)
	 		elif counter==2:
	 			self.message("processing %s file(s)..." % (jobs))
			task = GridTask(template="pick", name=self.stepname, command=k, cpu=cpus, dependson=list(), cwd = self.stepdir, debug=False)	
			tasks.append(task)
								
		for task in tasks:
			task.wait()	
			
class	AwkCommand(DefaultStep):
	def __init__(self, INS, ARGS, PREV,):
		DefaultStep.__init__(self)
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("AWK")
		#self.nodeCPUs=nodeCPUs
		self.start()
		
	def	performStep(self):
		tasks = list()
		
		oldtype = self.getInputValue("type")
		newtype = self.getInputValue("newtype")
		files = self.find(oldtype)
		awk = self.getInputValue("awk")
		
		postprocessing = ""
		if self.getInputValue("sort") != None:
			postprocessing = "%s | sort" % (postprocessing)
			
		if self.getInputValue("uniq") != None:
			postprocessing = "%s | uniq" % (postprocessing)	
			
		if self.getInputValue("postprocess") != None:
			postprocessing = "%s | %s " % (postprocessing, self.getInputValue("postprocess"))	
		
		counter=0
		for f in files:
			counter+=1
			newname = f[0:-len(oldtype)]
			newname = "%s%s" % (newname, newtype)
			k = "awk '%s' %s %s > %s" % (awk, f, postprocessing, newname)
			if counter==1:
		 		self.message(k)
	 		elif counter==2:
	 			self.message("processing %s file(s)..." % (len(files)))	
			task = GridTask(template="pick", name=self.stepname, command=k, cpu=1,  cwd = self.stepdir)
			tasks.append(task)
				
		for task in tasks:
			task.wait()	
			
class	ContaminantRemoval(DefaultStep):
	def __init__(self, INS, ARGS, PREV):		
		DefaultStep.__init__(self)
		
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("decontaminate")

	 	self.start()	
		 	
	def	performStep(self):
		tasks = list()

		m1 = self.find("fasta")
		m1.extend(self.find("mate1"))
		m1.extend(self.find("mate2"))
		
		filter = self.find("bowtie1alignment")	
		filter.extend(self.find("filter"))
		
		### index a filename using the filename sans the step id (0) and extension (-1)
		#filters = {".".join(key.strip().split(".")[1:-1]) : key for key in filter}
		
		filters = dict()
		for key in filter:
			filters[".".join(key.strip().split(".")[1:-1])] = key
		
		argstring = ""
		for arg, val in self.arguments.items():
			argstring = "%s %s %s " % (argstring, arg, val) 
		
		tasks = list()
		missing = 0
#		if len(filter)>1:
#			self.message("too many (%s) filters..." % (len(filter)))
#			self.failed = True
			
		for file in m1:
			
			### find appropriate filter	
			name = ".".join(file.strip().split(".")[1:-1])
			self.message("%s -> %s" % (name, filters[name]))
			
			if name in filters.keys():
				k = "%spython %sMateFilter.py %s -i %s -f %s " % (binpath, scriptspath, argstring, file, filters[name])
				if len(m1)==1:
					self.message(k)
			
				task = GridTask(template="pick", name="%s" % (self.stepname), command=k, cpu=1,  cwd = self.stepdir, debug=False)
				tasks.append(task)
			else:
				missing +=1	
			
			if missing>0:
				self.message("%s missing filters observed..." % (missing))
				self.failed = True
			
		for task in tasks:
			task.wait()	

class	SingletonsFishOut(DefaultStep):
	def __init__(self, INS, ARGS, PREV):		
		DefaultStep.__init__(self)
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("fishing")
	 	self.start()
		
	def	performStep(self):
		singletons = self.find("singletons")
		assembled = self.find("assembled")
		
		fasta = self.find("fasta")
		tasks = list()
		
		if len (singletons) !=0 and len(assembled) !=0:
			singletons = singletons[0]	
			assembled = assembled[0]
		
			for f in fasta:	
				
								
					if f.find("contigs")==-1:				
						k = "%spython %sMateFilter.py -i %s -k %s -t fasta -s singletons " % (binpath, scriptspath, f, singletons )
						self.message(k)
						task = GridTask(template="pick", name="%s" % (self.stepname), command=k, cpu=1,  cwd = self.stepdir, debug=False)
						task.wait()
						
						k = "%spython %sMateFilter.py -i %s -k %s -t fasta -s assembled " % (binpath, scriptspath, f, assembled )
						self.message(k)
						task = GridTask(template="pick", name="%s" % (self.stepname), command=k, cpu=1,  cwd = self.stepdir, debug=False)	
						task.wait()

						
				

class	fastx_quality_stats(DefaultStep):
	def __init__(self, INS, ARGS, PREV):		
		DefaultStep.__init__(self)
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("fastx_qstats")
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		
		m1 = self.find("mate1")
		m1.extend(self.find("mate2"))
				
		argstring = ""
		for arg, val in self.arguments.items():
			argstring = "%s %s %s " % (argstring, arg, val) 
			
		self.message(m1)
		for file in m1:
				
			k = "fastx_quality_stats %s -i %s -o %s.fastx_stats" % (argstring, file, file)
			self.message(k)
			task = GridTask(template="pick", name="%s" % (self.stepname), command=k, cpu=1,  cwd = self.stepdir)
			tasks.append(task)
			
		for task in tasks:
			task.wait()	
				
class	fastq_quality_filter(DefaultStep):
	def __init__(self, INS, ARGS, PREV):		
		DefaultStep.__init__(self)
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("fastq_qfilter")
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		
		m1 = self.find("mate1")
		m1.extend(self.find("mate2"))
				
		argstring = ""
		for arg, val in self.arguments.items():
			argstring = "%s %s %s " % (argstring, arg, val) 

		tasks = list()
		self.message("processing %s files..." % len(m1))
		for file in m1:
			suffix = file.split(".")[-1]
			prefix = ".".join(file.split(".")[:-1])
			k = "fastq_quality_filter %s -i %s -o %s.q.%s" % (argstring, file, prefix, suffix)
			if len(m1)<10:
				self.message(k)
			task = GridTask(template="pick", name="%s" % (self.stepname), command=k, cpu=1,  cwd = self.stepdir)
			tasks.append(task)
			
		for task in tasks:
			task.wait()	
			
class	fastq2fasta(DefaultStep):
	def __init__(self, INS, ARGS, PREV):		
		DefaultStep.__init__(self)
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("fastq2fasta")
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		
		m1 = self.find("fastq")
		if len(m1)==0:
			m1.extend(self.find("mate1"))
			m1.extend(self.find("mate2"))
				
		argstring = ""
		for arg, val in self.arguments.items():
			argstring = "%s %s %s " % (argstring, arg, val) 

		tasks = list()
		self.message("processing %s files..." % len(m1))
		for file in m1:
			suffix = file.split(".")[-1]
			prefix = ".".join(file.split(".")[:-1])
			if suffix=="fastq":
				k = "%sfastq_to_fasta %s -i %s -o %s.fasta" % (binpath, argstring, file, prefix)
			else:
				k = "%sfastq_to_fasta %s -i %s -o %s.fasta.%s" % (binpath, argstring, file, prefix, suffix)
			if len(m1)<10:
				self.message(k)
			task = GridTask(template="pick", name="%s" % (self.stepname), command=k, cpu=1,  cwd = self.stepdir, debug=True)
			tasks.append(task)
			
		for task in tasks:
			task.wait()	
		
class	mateInterweave(DefaultStep):
	def __init__(self, INS, ARGS, PREV):		
		DefaultStep.__init__(self)
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("mateInterweave")
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		m1 = self.find("mate1")
		m2 = self.find("mate2")		
		tasks = list()
		
		self.message("processing %s/%s files..." % (len(m1), len(m2)))
		
		for f in m1:
			f = ".".join(f.strip().split(".")[:-1])
			if "%s.mate2" %( f) in m2:
		
				argstring = ""
				for arg, val in self.arguments.items():
					argstring = "%s %s %s " % (argstring, arg, val) 
		
				argstring = "%s -f %s.mate1,%s.mate2" % (argstring, f, f)
				# if len(cluster)==2:
		# 			argstring = "%s -c %s " % (argstring, ",".join(cluster))
		# 		
				
				k = "%spython %sinterweaveMates.py %s" % (binpath, scriptspath, argstring)
				if len(m1)<10:
					self.message(k)
				task = GridTask(template="pick", name="%s" % (self.stepname), command=k, cpu=1,  cwd = self.stepdir)
				tasks.append(task)
			
		for task in tasks:	
			task.wait()
			
class	MateMerge(DefaultStep):
	def __init__(self, INS, ARGS, PREV, prefix="files"):
	
		DefaultStep.__init__(self)
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("MATE_cat")
		#self.nodeCPUs=nodeCPUs
		self.prefix = prefix
	 	self.start()
		
	def	performStep(self):
		m1 = self.find("mate1")
		m2 = self.find("mate2")	
		
		k = "cat *.fasta.mate1 *.fasta.mate2 > %s.mates.fasta" % (self.prefix)
		self.message(k)	
		task = GridTask(template="pick", name="cat", command=k, cpu=1,  cwd = self.stepdir)
		task.wait()		
						
class	CLC_Assemble(DefaultStep):
	def __init__(self, INS, ARGS, PREV):		
		DefaultStep.__init__(self)
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("CLC_Assemble")
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		x= self.find("properties")
		m1 = self.find("fasta")
		cpus=1
		template="pick"
		if len(x)!=1 or len(m1)!=1:
			self.failed=True
		else:
			m1 = m1[0]			
			prefix = ".".join(m1.split(".")[:-1])
			argstring=""
			for arg, val in self.arguments.items():
				argstring = "%s %s %s " % (argstring, arg, val) 
				if arg=="--cpus":
					cpus=val
					
			done = False
			
			while not done:
			
				k = "/usr/local/packages/clc-ngs-cell/clc_novo_assemble -o %s.contigs.fasta %s -q %s" % (prefix, argstring, m1)
				self.message(k)
				
				task = GridTask(template="pick", name="%s" % (self.stepname), command=k, cpu=cpus,  cwd = self.stepdir, debug=True)
				task.wait()
				
				for file in glob.glob("%s/*.e*" % (self.stepdir)):
					#self.message(file)
					contents = "\n".join(loadLines("%s" % (file)))
					if contents.find("No more available licenses")>-1:	
						self.message("No more available licenses, retrying in a bit...")
						command = "rm %s" % (file)
						p = Popen(shlex.split(command), close_fds=True)
						p.wait()
						time.sleep(60)
					else:
						done = True
								
class	CLC_Assemble_Ref(DefaultStep):
	def __init__(self, INS, ARGS, PREV):		
		DefaultStep.__init__(self)
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("CLC_Assemble_Ref")
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		x= self.find("properties")
		cpus=1
		template="pick"
		
		if len(x)!=1:
			self.failed=True
		else:
			fastas = self.find("fasta")
			self.message(fastas)
			contigs = list()
			reads = list()
			if len(fastas)==2:
				
				for f in fastas:
					if f.find("contig")>-1:
						contigs.append(f)
					else:
						reads.append(f)
				
				if len(contigs)==1 and len(reads)==1:
					contigs = contigs[0]
					reads =reads[0]	
					prefix = ".".join(reads.split(".")[:-1])
					argstring=""
					for arg, val in self.arguments.items():
						argstring = "%s %s %s " % (argstring, arg, val) 
						if arg=="--cpus":
							cpus=val
							
					### clean reads		
					done = False
					while not done:
						k = "/usr/local/packages/clc-ngs-cell/clc_ref_assemble_long  %s -o %s.clean.cas -q %s -d %s" % (argstring, prefix, reads, contigs)
						self.message(k)
						
						task = GridTask(template="pick", name="%s" % (self.stepname), command=k, cpu=cpus,  cwd = self.stepdir, debug=True)
						task.wait()
						
						for file in glob.glob("%s/*.e*" % (self.stepdir)):
							#self.message(file)
							contents = "\n".join(loadLines("%s" % (file)))
							if contents.find("No more available licenses")>-1:	
								self.message("No more available licenses, retrying in a bit...")
								command = "rm %s" % (file)
								p = Popen(shlex.split(command), close_fds=True)
								p.wait()
								time.sleep(60)
							else:
								done = True
					
				else:
					self.failed=True			
			else:
				self.failed=True

class	CLC_Assemble_Info(DefaultStep):
	def __init__(self, ARGS, PREV):
		DefaultStep.__init__(self)
		#self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("CLC_Assemble_Info")
		#self.nodeCPUs=nodeCPUs
	 	self.start()
		
	def	performStep(self):
		for f in self.find("cas"):	
			k = "/usr/local/packages/clc-ngs-cell/assembly_info %s > %s.clcassemblystats" % (f,f)
			self.message(k)
			task = GridTask(template="pick", name=self.stepname, command=k, cwd = self.stepdir)
			task.wait()

class	ContigCoverageUpdate(DefaultStep):
	def __init__(self, ARGS, PREV):
		DefaultStep.__init__(self)
		#self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("ContigCovUp")
		#self.nodeCPUs=nodeCPUs
	 	self.start()
		
	def	performStep(self):
		table = self.find("clcassemblystats")
		fasta = self.find("fasta")
		tasks = list()
		for f in fasta:
			#self.message(f)
			if f.find("contigs")>-1:
				for t in table: 
					if ".".join(f.split(".")[1:-3]) in t:
						k = "%spython %sContigCoverageUpdate.py %s  %s" % (binpath, scriptspath, t, f)
						self.message(k)
						task = GridTask(template="pick", name=self.stepname, command=k, cwd = self.stepdir, debug=False)
						tasks.append(task)
		for t in tasks:
			t.wait()
		#self.failed=True	

class	ORFCoverage(DefaultStep):
	def __init__(self, INS, ARGS, PREV):
		DefaultStep.__init__(self)
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("ORFCoverage")
		#self.nodeCPUs=nodeCPUs
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		table = self.find("clcassemblytable")[0]
		fastas = self.find("fasta")
		fasta= ""
		for f in fastas:
			if f.find("orf")>-1:
				fasta = f
				id = fasta.strip().split(".")[1]

		k = "%spython %sORFCoverage.py -o %s -a %s -e %s.orfs.weight" % (binpath, scriptspath, fasta, table, id)
		self.message(k)
		task = GridTask(template="pick", name=self.stepname, command=k, cwd = self.stepdir)
		tasks.append(task)
		
		for t in tasks:
			t.wait()
			
class	ORFCoverageNorm(DefaultStep):
	def __init__(self, INS, ARGS, PREV):
		DefaultStep.__init__(self)
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("ORFCoverageNorm")
		#self.nodeCPUs=nodeCPUs
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		weight = self.find("weight")
		
		k = "%sR CMD BATCH %sORFweights.r" % (binpath, scriptspath)
		self.message(k)
		task = GridTask(template="pick", name=self.stepname, command=k, cwd = self.stepdir)
		tasks.append(task)
		for t in tasks:
			t.wait()			

class	PROKModify(DefaultStep):
	def __init__(self, INS, ARGS, PREV):
		DefaultStep.__init__(self)
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("ProkModify")
		#self.nodeCPUs=nodeCPUs
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		w = self.find("normweight")[0]
		p = self.find("txt")[0]
		id = p.strip().split(".")[1]
		
		k = "%spython %sORFCoverageModifyProk.py -p %s -w %s -o %s.jcvinorm" % (binpath, scriptspath, p, w, id)
		self.message(k)
		
		task = GridTask(template="pick", name=self.stepname, command=k, cwd = self.stepdir)
		tasks.append(task)
	
		for t in tasks:
			t.wait()
			
		
class	CLC_Assemble_Table(DefaultStep):
	def __init__(self, ARGS, PREV):
		DefaultStep.__init__(self)
		#self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("CLC_Assemble_Table")
		#self.nodeCPUs=nodeCPUs
	 	self.start()
		
	def	performStep(self):
		x = self.find ("fasta")
		tasks = list()
		
		for f in self.find("cas"):
			k = "/usr/local/packages/clc-ngs-cell/assembly_table -n -p -s %s > %s.clcassemblytable" % (f,f)
			self.message(k)
			task = GridTask(template="pick", name=self.stepname, command=k, cwd = self.stepdir)
			task.wait()	
			
			k = "awk ' $5 == -1 && $4 == -1 {print $2}' %s.clcassemblytable > %s.singletons" % (f,f)
			self.message(k)	
			task = GridTask(template="pick", name=self.stepname, command=k, cwd = self.stepdir)
			tasks.append(task)
			
			k = "awk ' $5 > -1 && $4 > -1 {print $2}' %s.clcassemblytable > %s.assembled" % (f,f)
			self.message(k)
			task = GridTask(template="pick", name=self.stepname, command=k, cwd = self.stepdir)
			tasks.append(task)	
			
		for task in tasks:
			task.wait()
							
class	FastaSummaryRPlots(DefaultStep):
	def __init__(self, ARGS, PREV):
		DefaultStep.__init__(self)
		#self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("FastaSummary")
		#self.nodeCPUs=nodeCPUs
	 	self.start()
		
	def	performStep(self):
		self.find("fasta")
		k = "%sR CMD BATCH %sStatFastaFiles.R" % (binpath, scriptspath)
		self.message(k)
		task = GridTask(template="pick", name=self.stepname, command=k, cwd = self.stepdir)
		task.wait()
		
class	ClearcutTree(DefaultStep):
	def __init__(self, ARGS, PREV):
		DefaultStep.__init__(self)
		#self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("ClearcutTree")
		#self.nodeCPUs=nodeCPUs
	 	self.start()
		
	def	performStep(self):
		x = self.find("fasta")
		for fasta in x:
			k = "%sclearcut --alignment --DNA --in=%s --out=%s.tre" % (binpath, fasta, fasta)
			self.message(k)
			task = GridTask(template="himem.q", name=self.stepname, command=k, cwd = self.stepdir)
			task.wait()		

class	SQA(DefaultStep):
	def __init__(self, ARGS, PREV):
		DefaultStep.__init__(self)
		#self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("SolexaQA")
		#self.nodeCPUs=nodeCPUs
	 	self.start()
		
	def	performStep(self):
		files= list()
		for type in ("mate1", "mate2", "fastq"):
			tmp = self.find(type)
			if tmp!=None:
				files.extend(tmp)
		
		argstring = ""
		for arg, val in self.arguments.items():
			argstring = "%s %s %s " % (argstring, arg, val) 
		
		tasks=list()
		for f in files:
			k = "perl %sSolexaQA.pl %s %s" % (sqapath, f, argstring)
			self.message(k)
			task = GridTask(template="pick", name=self.stepname, command=k, cwd = self.stepdir, debug=False)
			tasks.append(task)
			
		for task in tasks:
			task.wait()	
			
class	SQAtrim(DefaultStep):
	def __init__(self, ARGS, PREV):
		DefaultStep.__init__(self)
		#self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("DynamicTrim")
		#self.nodeCPUs=nodeCPUs
	 	self.start()
		
	def	performStep(self):
		files= list()
		for type in ("mate1", "mate2", "fastq"):
			tmp = self.find(type)
			if tmp!=None:
				files.extend(tmp)
		
		argstring = ""
		for arg, val in self.arguments.items():
			argstring = "%s %s %s " % (argstring, arg, val) 
		
		jobs = len(files)
		counter=0
		tasks=list()
		for f in files:
			counter+=1
			tokens = f.strip().split(".")
			newf = "%s.trimmed.%s" % ( ".".join(tokens[:-1]), tokens[-1])	
				
			k = "perl %sDynamicTrim.pl %s %s; mv %s.trimmed %s" % (sqapath, f, argstring, f, newf)
			if counter==1:
				self.message(k)
			elif counter==2:
				self.message("processing %s files..." % (jobs))	
			task = GridTask(template="pick", name=self.stepname, command=k, cwd = self.stepdir, debug=False)
			tasks.append(task)
			
		for task in tasks:
			task.wait()	

class	SQAlenfil(DefaultStep):
	def __init__(self, ARGS, PREV):
		DefaultStep.__init__(self)
		#self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("LengthSort")
		#self.nodeCPUs=nodeCPUs
	 	self.start()
		
	def	performStep(self):
		M1 = list()
		M2 = list()
		pairedindex = defaultdict(list)
		singlefiles = list()
		
		for tp in ("mate1", "mate2"):
			tmp = self.find(tp)
			if tmp!=None:
				if tp.endswith("1"):
					M1.extend(tmp)
				elif tp.endswith("2"):
					M2.extend(tmp)	
				
		for f1 in M1:
			core1 = f1.strip().split(".")[1:-1]
			for f2 in M2:
				core2 = f2.strip().split(".")[1:-1]
				if core1==core2:
					pairedindex[f1].append(f2)
				
		
		
		for type in ("fastq"):
			tmp = self.find(type)
			if tmp!=None:
				singlefiles.extend(tmp)		
				
		
		argstring = ""
		for arg, val in self.arguments.items():
			argstring = "%s %s %s " % (argstring, arg, val) 
				
			
		tasks=list()
		for f1, f2s in pairedindex.items():
			if len(f2s)!=1:
				singlefiles.append(f1)
				singlefiles.extend(f2s)
			else:	
				newf1 = "%s.len.mate1" % (".".join(f1.strip().split(".")[:-1]))
				newf2 = "%s.len.mate2" % (".".join(f2s[0].strip().split(".")[:-1]))

				k = "perl %sLengthSort.pl %s %s %s; mv %s.paired1 %s; mv %s.paired2 %s" % (sqapath, f1, f2s[0], argstring, f1, newf1, f1, newf2)
				#k = "perl %sLengthSort.pl %s %s %s" % (sqapath, f1, f2s[0], argstring)
				
				self.message(k)
				task = GridTask(template="pick", name=self.stepname, command=k, cwd = self.stepdir, debug=False)
				tasks.append(task)
			
		for f in singlefiles:
			newf = "%s.len.fastq" % (".".join(f.strip().split(".")[:-1]))
			k = "perl %sLengthSort.pl %s %s; mv %s.single %s" % (sqapath, f, argstring, f, newd)
			self.message(k)
			task = GridTask(template="pick", name=self.stepname, command=k, cwd = self.stepdir, debug=False)
			tasks.append(task)
			
		for task in tasks:
			task.wait()							


class	GuessFastQEncoding(DefaultStep):
	def __init__(self, ARGS, PREV):
		DefaultStep.__init__(self)
		#self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("FastQEnc")
		#self.nodeCPUs=nodeCPUs
	 	self.start()
		
	def	performStep(self):
	
		x = self.find("mate1")
		k = "%spython %sFastQEncoding.py %s > %s.offset" % (binpath, scriptspath, x[0], x[0])
		self.message(k)
		task = GridTask(template="pick", name=self.stepname, command=k, cwd = self.stepdir, debug=False)
		task.wait()
		
		otpt = ""
		for line in loadLines("%s/%s.offset" % (self.stepdir, x[0])):
			otpt = "%s%s" % (otpt, line.strip())
		
		self.message("%s -> %s" % (x[0], otpt)) 
		self.setOutputValue("-Q", otpt)
			
class	MascotReportLifter(DefaultStep):
	def __init__(self, ARGS, PREV):
		DefaultStep.__init__(self)
		#self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("MascotLifter")
		#self.nodeCPUs=nodeCPUs
	 	self.start()
		
	def	performStep(self):
		one_at_a_time.acquire()
		
		file = self.getInputValue("file")
		id = self.getInputValue("id")
		script =  "/usr/local/projects/T1D-CNMC/sszpakow/MASCOT_FANTOMJS/MascotAutomaton.js"
		self.message("caching and reporting on %s in file %s" % (id, file))
		
		k = "/home/sszpakow/bin/phantomjs %s %s %s" % (script, file, id)
		self.message(k)
		
		#task = GridTask(template="default", name=self.stepname, command=k, cwd = self.stepdir, debug=True)
		#task.wait()
		
		p = Popen(shlex.split(k), stdout = PIPE, stderr = PIPE, close_fds=True, cwd=self.stepdir)
		out, err = p.communicate()
		p.wait()
		
		f = open("%s/%s_%s.out.log" % (self.stepdir, file.strip().split("/")[-1], id ), "w")
		f.write(out)
		f.write("\n")
		f.close()
		
		f = open("%s/%s_%s.err.log" % (self.stepdir, file.strip().split("/")[-1], id ), "w")
		f.write(err)
		f.write("\n")
		f.close()
		
		if err.find("'waitFor()' timeout")>-1 or out.find("'waitFor()' timeout"):
			self.message("Time out detected! %s - %s" % (id, file) )
			#self.failed=True
		
		one_at_a_time.release()

class	SED_replace(DefaultStep):
	def __init__(self, INS, ARGS, PREV,):
		DefaultStep.__init__(self)
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("SED_replace")
		#self.nodeCPUs=nodeCPUs
		self.start()
		
	def	performStep(self):
		tasks = list()
		
		for t in self.getInputValue("types").strip().split(","):
			files = self.find(t)
			old = self.getInputValue("old")
			new = self.getInputValue("new")
			
			for f in files:
				newname = f[0:-len(t)]
				newname = "%str.%s" % (newname, t)
				k = "sed 's/%s/%s/g' %s > %s" % (old, new, f, newname)
				self.message(k)	
				task = GridTask(template="pick", name=self.stepname, command=k, cpu=1,  cwd = self.stepdir)
				tasks.append(task)
			
		for task in tasks:
			task.wait()			


class	FileMiniImport(DefaultStep):	
	def	__init__(self, INS, ARGS):
		DefaultStep.__init__(self)
		self.setInputs(INS)
		self.setArguments(ARGS)
		#self.setPrevious(PREV)
		self.setName("FILE_mini_input")
	 	self.start()
			
	def	performStep(self):
		
		lines = self.getInputValue("lines")
		if lines == None:
			lines = 100000
			
		for type in self.inputs.keys():
			files = self.inputs[type]
			for file in files:
				pool_open_files.acquire()
				file = file.split("~")
				if len(file)>1:
					file, newname = file
					tmp = file.strip().split("/")[-1]
					k = "head -n %s %s" % (lines, file)
					outname = "%s.%s" % (newname, type)
					
				else:
					file = file[0]
					tmp = file.strip().split("/")[-1]
					k ="head -n %s %s "	% (lines, file, tmp)
					outname = "imported.%s" % (tmp)	
								
				p = Popen(shlex.split(k), stdout=PIPE, stderr=PIPE, cwd=self.stepdir, close_fds=True)
				self.message(k)
				out,err = p.communicate()
				p.wait()
				
				o = open("%s/%s" % (self.stepdir, outname), "w")
				o.write(out)
				o.close()
				
				pool_open_files.release()
					
				

class	FileSplit(DefaultStep):
	def __init__(self, ARGS, PREV):		
		DefaultStep.__init__(self)
		#ARGS = 	{"types": TYPES}
		#self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("FILE_split")
		#self.nodeCPUs=nodeCPUs
		self.chunk = 0
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		self.chunk = self.getInputValue("chunk")
		for t in self.getInputValue("types").strip().split(","):
			files = self.find(t)
			for f in files:
				k = "split -a 5 -l %s %s %s.split. " % (self.chunk, f, f)
				self.message(k)	
				task = GridTask(template="pick", name=self.stepname, command=k, cpu=1,  cwd = self.stepdir)
				tasks.append(task)
		for task in tasks:
			task.wait()
			time.sleep(1)	
			
		### rename the files
		for file in glob.glob("%s/*.split.*" % (self.stepdir)):
			newfile = file.strip().split(".")
			suffix = newfile[-3:]
			suffix.reverse()
 			newfile = "%s.%s" % ( ".".join(newfile[:-3]), ".".join(suffix))
 			#self.message("%s -> %s" % (file, newfile) )
			command = "mv %s %s" % (file, newfile)
			p = Popen(shlex.split(command), close_fds=True)
			p.wait()
			
class	FastaSplit(DefaultStep):
	def __init__(self, ARGS, PREV):		
		DefaultStep.__init__(self)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("FASTA_split")
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		chunk = self.getInputValue("chunk")
		for t in self.getInputValue("types").strip().split(","):
			files = self.find(t)
			for f in files:
				k = "{0}python ~/scripts/python/FastaSplitter.py  -f {1} -c {2}".format(binpath, f, chunk)
				self.message(k)	
				task = GridTask(template="pick", name=self.stepname,  command=k, cpu=1,  cwd = self.stepdir)
				tasks.append(task)
		for task in tasks:
			task.wait()
			
class	FastaSort(DefaultStep):
	def __init__(self, ARGS, PREV):		
		DefaultStep.__init__(self)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("FASTA_sort")
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		for t in self.getInputValue("types").strip().split(","):
			files = self.find(t)
			for f in files:
				k = "{0}python ~/scripts/python/FastaSort.py {1}".format(binpath, f)
				self.message(k)	
				task = GridTask(template="pick", name=self.stepname,  command=k, cpu=1,  cwd = self.stepdir)
				tasks.append(task)
		for task in tasks:
			task.wait()			

class	BLAST(DefaultStep):
	def __init__(self, ARGS, PREV):		
		DefaultStep.__init__(self)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("BLAST")
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		cpus = 0
		arguments = ""
		
		mode = "blastn"
		
		for arg, val in self.arguments.items():
			
			if arg == "mode":
				mode = val;
			
			else:
			
				if arg == "-num_threads":
					cpus = val
					
				arguments = "{0} {1} {2} ".format(arguments, arg, val) 	
			
		if cpus ==0:
			cpus = 4
			arguments = "{0} -num_threads 4 ".format(arguments)	
		
		for f in self.find("fasta"):
			k = "/usr/local/packages/ncbi-blast+/bin/{0} {1} -query {2} -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\"  -out {2}.blast6".format(mode, arguments, f)
			self.message(k)	
			task = GridTask(template="pick", name=self.stepname,  command=k, cpu=cpus,  cwd = self.stepdir, debug=True)
			tasks.append(task)
		
		for task in tasks:
			task.wait()	
			
		#self.failed = True

class	FileTypeTrim(DefaultStep):
	def __init__(self, ARGS, PREV):		
		DefaultStep.__init__(self)
		#self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("FILE_typetrim")
		#self.nodeCPUs=nodeCPUs
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		for input in self.arguments.keys():
			files = self.find(input)
			for file in files:
				outname = "%s" % (file[:-len(input)])
				outname = outname.strip(".")
				k = "cp %s %s" % (file, outname)
				self.message(k)	
				task = GridTask(template="pick", name=self.stepname, command=k, cpu=1,  cwd = self.stepdir)
				tasks.append(task)
		for task in tasks:
			task.wait()
			time.sleep(1)

class	Flash (DefaultStep):
	def __init__(self, INS, ARGS, PREV):		
		DefaultStep.__init__(self)
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("Flash")
		self.counter=0
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		
		m1 = self.find("mate1")
		m2 = self.find("mate2")
		argstring = ""
		for arg, val in self.arguments.items():
			argstring = "%s %s %s " % (argstring, arg, val) 
	
		tasks = list()
		jobs = len(m1)
		counter=0
		for f in m1:
			counter+=1
			name = ".".join(f.strip().split(".")[:-1])

			if "%s.mate2" % (name) in m2:
				k = "%sflash %s.mate1 %s.mate2 %s -o %s; mv %s.notCombined_1.fastq %s.notC.mate1; mv %s.notCombined_2.fastq %s.notC.mate2" % (binpath, name, name, argstring, name, name, name, name, name)
		 		if counter==1:
		 			self.message(k)
		 		elif counter==2:
		 			self.message("processing %s file(s)..." % (jobs))
		 			
		 		task = GridTask(template="pick", name=self.stepname, command=k, cpu=1, dependson=list(), cwd = self.stepdir, debug=False)	
				tasks.append(task)
			else:
				self.message("skipping: %s" % (name))
								
		for task in tasks:
			task.wait()
	
class	PrimerClipper(DefaultStep):
	def __init__(self, INS, ARGS, PREV):		
		DefaultStep.__init__(self)
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("PrimerClipper")
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		
		m1 = self.find("fasta")
				
		argstring = ""
		for arg, val in self.arguments.items():
			argstring = "%s %s %s " % (argstring, arg, val) 

		tasks = list()
		self.message("processing %s files..." % len(m1))
		for file in m1:
			suffix = file.split(".")[-1]
			prefix = ".".join(file.split(".")[:-1])
			k = "%spython %sPrimerClipper.py %s -i %s" % (binpath, scriptspath, argstring, file)
			if len(m1)<10:
				self.message(k)
			task = GridTask(template="pick", name="%s" % (self.stepname), command=k, cpu=1,  cwd = self.stepdir, debug=True)
			tasks.append(task)
			
		for task in tasks:
			task.wait()	
		if len(m1)==0:
			self.message("No files for clipping...")	
			
class	FastaHeadHash(DefaultStep):
	def __init__(self, INS, ARGS, PREV):		
		DefaultStep.__init__(self)
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("FastaHeadHash")
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		
		m1 = self.find("fasta")
				
		argstring = ""
		for arg, val in self.arguments.items():
			argstring = "%s %s %s " % (argstring, arg, val) 

		tasks = list()
		self.message("processing %s files..." % len(m1))
		for file in m1:
			suffix = file.split(".")[-1]
			prefix = ".".join(file.split(".")[:-1])
			k = "%spython %sFastaUniversalRenamer.py %s --fasta %s" % (binpath, scriptspath, argstring, file)
			if len(m1)<10:
				self.message(k)
			task = GridTask(template="pick", name="%s" % (self.stepname), command=k, cpu=1,  cwd = self.stepdir)
			tasks.append(task)
			
		for task in tasks:
			task.wait()				


class	OtuTable(DefaultStep):
	def __init__(self, INS, ARGS, PREV):		
		DefaultStep.__init__(self)
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("OtuTable")
	 	self.start()
		
	def	performStep(self):
		
		list = self.find("list")[0]
		group = self.find("group")[0]	
		k = "%spython %sOTUtableMaker.py -l %s -g %s " % (binpath, scriptspath, list, group )
		self.message(k)
		task = GridTask(template="pick", name="%s" % (self.stepname), command=k, cpu=1,  cwd = self.stepdir, debug=False)
		task.wait()
						
#################################################
##		FUNCTIONS
#################################################

def getQ(file):
	k = "%spython %sFastQEncoding.py %s" % (binpath, scriptspath, file)
	p = Popen(shlex.split(k), stdout=PIPE, stderr=PIPE, close_fds=True)
	out,err = p.communicate()
	p.wait()
	return "%s" % (out.strip())

#################################################
##		ARGS
#################################################

one_at_a_time = BoundedSemaphore(value=2, verbose=False)
sqapath = "/usr/local/devel/ANNOTATION/sszpakow/YAP/bin/solexaQA-current/"


#################################################
##		Finish
#################################################
