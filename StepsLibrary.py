########################################################################################
## This file is a part of YAP package of scripts. https://github.com/shpakoo/YAP
## Distributed under the MIT license: http://www.opensource.org/licenses/mit-license.php
## Copyright (c) 2011-2012 J.Craig Venter Institute.
########################################################################################

#################################################
## A library of "steps" or program wrappers to construct pipelines
## Pipeline steps orchestration, grid management and output handling.
#################################################
import sys, tempfile, shlex, glob, os, stat, hashlib, time, datetime, re, curses
from threading import *
from subprocess import *
from MothurCommandInfoWrapper import *
from collections import defaultdict
from collections import deque
from random import *
from Queue import *

_author="Sebastian Szpakowski"
_date="2012/09/20"
_version="Version 2"


#################################################
##		Classes
##

class	BufferedOutputHandler(Thread):
	def __init__(self, usecurses=False):
		Thread.__init__(self)
		self.shutdown=False
		self.cache = deque()
		self.registered=0
	
		self.ids = list()
		self.wrap = 140
		
		self.starttime = time.time()
		
		#### init log
		try:
			self.otptfile = open("logfile.txt", 'a')
			self.toPrint("-----", "GLOBAL", "Appending to a logfile.txt...")
		except:
			self.otptfile = open("logfile.txt", 'w')
			self.toPrint("-----", "GLOBAL", "Creating a new logfile.txt...")
		
		command = " ".join(sys.argv)
		self.otptfile.write("command: %s\n" % command)
		
		#### init output (curses)
		
		
		self.usecurses = usecurses
		
		if (self.usecurses):
			
			self.stdscr=curses.initscr()
			curses.savetty()
			curses.noecho()
			curses.cbreak()
			curses.curs_set(0) 
 
			self.textbuffer= list()
			self.stdscr.clear()
			self.stdscr.refresh()
			self.cursestrackbuffer = 100
			self.scrollpad = curses.newpad(self.cursestrackbuffer*2, self.wrap*2)
			self.spacerpad = curses.newpad(1,1000)
			self.updatepad = curses.newpad(10,1000)
			self.rows, self.cols = self.stdscr.getmaxyx()
			
		else:
			self.stdscr=None
			
		self.start()
		
	def	run(self):
		self.toPrint("-----", "GLOBAL", "Setting up the pipeline...")
		self.flush()
		
		time.sleep(5)
		
		while activeCount()>3 or self.registered>0 or len(self.cache) > 0:
			self.flush()
			time.sleep(1)
		
		self.flush()
		endtime = time.time()
		text =  "+%s [fin]" % (str(datetime.timedelta(seconds=round(endtime-self.starttime,0))).rjust(17)) 
		self.toPrint("-----", "GLOBAL", text)	
		
		command = "python %straverser.py" % (scriptspath)
		p = Popen(shlex.split(command), stdout = PIPE, stderr = PIPE, close_fds=True)
		dot, err = p.communicate()
		p.wait()
		
		x = open("workflow.dot", "w")
		x.write(dot)
		x.write("\n")
		x.close()
		
		for format in ["svg", "png", "pdf", "jpg"]:
			command = "dot -T%s -o workflow.%s" % (format, format) 
			p = Popen(shlex.split(command), stdin = PIPE, stdout = PIPE, stderr = PIPE, close_fds=True)
			out, err = p.communicate(dot)
			p.wait()
		
		self.toPrint("-----", "GLOBAL", "Check out workflow.{svg,png,jpg} for an overview of what happened.")
		self.flush()
		self.otptfile.close()
		self.closeDisplay()
		
	def	register(self, id):
		self.registered+=1
		self.ids.append(id)
	
	def	deregister(self):
		self.registered-=1
	
	def	collapseIDs(self, text ):
		for id in self.ids:
			if len(id)>5:
				text = re.sub(id, "[{0}~]".format(id[:5]), text)
		return (text)	
	
	def flush(self):
			
		while len(self.cache) > 0:
			id, name, line = self.cache.popleft()
			tag = "[{0}] {1:<20} > ".format( id[:5], name) 
			line = "{0!s}".format(line)
			line = self.collapseIDs(line)
			
			
			otpt = "{0}{1}".format(tag, line[:self.wrap])
			self.otptfile.write("{0}{1}\n".format(tag, line))
			
			line = line[self.wrap:]
			self.outputScroll(otpt)
			
			
			while len(line)>=self.wrap:
				otpt = "{0}\t{1}".format(tag, line[:self.wrap])
				line = line[self.wrap:]
				self.outputScroll(otpt)	
				
			if len(line)>0:	
				otpt = "{0:<30}\t\t{1}".format("", line)
				line = line
				self.outputScroll(otpt)	
		
		self.redrawScreen()						
	
	def	redrawScreen(self):
		try:
			y,x = self.stdscr.getmaxyx()
			### enough screen to print:
			if y>20 and x>20:
				
				if len(self.textbuffer) < (y-10):
					self.scrollpad.refresh(0, 0, 0, 0, y-10, x-5)
				else:
					self.scrollpad.refresh(self.cursestrackbuffer-y+10 , 0, 0, 0, y-10, x-5)	
				
				self.updatepad.refresh(0, 0, y-8, 10 , y-3, x-5)
				
			### when screen too small
			else:	
				self.scrollpad.refresh(0,0,0,0,0,0)
				self.updatepad.refresh(0,0,0,0,0,0)
		except:
			self.closeDisplay()
			self.usecurses=False		
#														
	def	toPrint(self, id, name, line):
		self.cache.append((id, name, line))
	
	def outputScroll(self, k):	
		if self.usecurses:
			self.textbuffer.append("%s\n" %(k))
			self.scrollpad.clear()
			for k in self.textbuffer[-self.cursestrackbuffer:]:
				self.scrollpad.addstr(k)		
		else:
			print k	
		
	def outputUpdate(self,k):
		if self.usecurses:
			self.updatepad.clear()
			for k in k.strip().split("\n"):
				self.updatepad.addstr("%s\n" % k)
				
			
	def closeDisplay(self):	
		if self.usecurses:
			
			self.stdscr.clear()
			self.stdscr.refresh()
			 
			curses.curs_set(1)
			curses.nocbreak()
			curses.echo()
			curses.resetty()
			curses.endwin()
			

class   TaskQueueStatus(Thread):
	def __init__(self, update=1, maxnodes=10):
		Thread.__init__(self)
		self.active=True
		
		self.maxnodes = maxnodes 
		self.available = self.maxnodes
		
		self.update = update	
			
		#### queue of grid jobs to run
		self.scheduled = Queue() 
		#### to keep track of things popped off the queue
		self.processing = dict()
		
		#### inventory of what ran
		#### tuple (jid, status) indexed by command
		#### status: new/running/done/remove
		#### new		 upon registering
		#### running	 when submitted to the grid
		#### done		 when completed
		
		self.registered = dict()
				
		#### inventory of completed jobs		
		self.bestqueue = "default.q"
		self.pollqueues()
		
		self.running=0
		self.stats=dict()
		
		self.previous =""
		
		self.start()
			
	def run(self):
		BOH.outputUpdate("Setting up the grid...")
		print "Setting up grid..."
		time.sleep(5)
		while activeCount()>3 or self.running>0 or self.scheduled.qsize()>0:

			self.pollfinished()
			self.pollqueues()
			self.pollrunning()
			self.dispatch()
			self.cleanup()
			
			BOH.outputUpdate("%s" % (self))
			#print self
			
			time.sleep(self.update)
		
		BOH.outputUpdate("%s\nGrid Offline." % (self))
		
		print self	
		print "Queue status shutting down."
		
	def cleanup(self):
		toremove = set()
		for key, tup in self.registered.items():
			id, status = tup
			if status == "remove":
				toremove.add(key)
		for key in toremove:
			del self.registered[key]
	
	def flagRemoval(self, task):
		id, status = self.registered[task.getUniqueID()]
		if status =="done":
			 self.registered[task.getUniqueID()] = [id, "remove"]
		else:
			print "cannot flag yet:", id, status
			
	def pollfinished(self):
					
#		donejobs = set() 
#		
#		### only 100 recent jobs shown, which could be a problem ;-)   
#		p = Popen(shlex.split("qstat -s z"), stdout=PIPE, stderr=PIPE, close_fds=True)
#		p.wait()
#		out,err = p.communicate()
#
#		lines = out.split("\n")
#		tmp = set()
#		if len(lines)>2:
#			for line in lines[2:]:
#				line = line.strip().split()
#				if len(line)>0:
#					donejobs.add(line[0])		
#		
		#if len(donejobs)>0:
		for key, tup in self.registered.items():
			id, status = tup
			#if (status == "running") and (id in donejobs):
			if (status == "running") and (self.isJobDone(id)):
				tmp = self.registered[key][1]= "done"
				self.processing[key].setCompleted()
				self.available += 1
				del self.processing[key]
				
	def isJobDone(self, jid):
		p = Popen(shlex.split("qstat -j %s" % jid), stdout=PIPE, stderr=PIPE, close_fds=True)
		p.wait()
		out,err = p.communicate()		
		return err.find("jobs do not exist")>-1
		   
	def pollqueues(self):
		command="qstat -g c" 
		p = Popen(shlex.split(command), stdout=PIPE, stderr=PIPE, close_fds=True )
		p.wait()
		out,err = p.communicate()
		
		if err.find("neither submit nor admin host")==-1:
			queues = dict()
			out = out.strip().split("\n")
			for q in out[2:]:
				queue, cqload, used, res, avail, total, acds, cdsu = q.split()
				avail = float(avail)
				total = float(total)
				if total>0:
					queues[queue] = avail
			
			#for k in ("default.q", "medium.q", "fast.q", "himem.q"):
			for k in ("himem.q", "medium.q", "default.q"):
			#for k in ("fast.q", "medium.q", "default.q"):
				if queues[k] > queues[self.bestqueue]:
					self.bestqueue= k
			  
### sanity check, this should match the counters				
	def pollrunning(self):
		tmp=defaultdict(int)
		for jid, value in self.registered.values():
			tmp[value]+=1
		self.stats = tmp	
		self.running = self.stats["running"]
	
	def dispatch(self): 
	
		while self.nodesAvailable():
			if not self.scheduled.empty():
				
				tmp = self.scheduled.get()
				self.processing[tmp.getUniqueID()]=tmp
				#print "submitting", tmp.getUniqueID()
				
				
				jid = tmp.submit()
				#print jid
				
				if jid==-1:
					print "???", tmp
				
		
				
				self.registered[tmp.getUniqueID()] = [tmp.getGridId(), "running"]	
				self.available-=1
			else:
				break
			
			
		
	def pickQ(self):
		return self.bestqueue
	
	def register(self, task):
		self.scheduled.put(task)
		self.registered[task.getUniqueID()]=[-1, "new"]
								
	def shutdown(self):
		self.active=False
		print "Queue status shutting down..."
	
	def nodesAvailable(self):
		return (self.available > 0)
	
	def __str__(self):
		otpt ="Currently running/waiting: %s/%s\n" % (self.running, self.scheduled.qsize())
		otpt ="%savailable/total: %s/%s" % (otpt, self.available, self.maxnodes)
		# for key, tup in self.registered.items():
#			 id, status = tup
#			 if id != -1:
#				 otpt = "%s\n\t%s\t%s\t%s" % (otpt, id, status, key[0:10])
		for key, val in self.stats.items():
			otpt = "%s\n\t%s\t%s" % (otpt, key, val) 

		otpt = "%s\n\nbest queue: %s" % (otpt, self.bestqueue)
		return (otpt)	


	#################################################
	### a thread that will track of a qsub job
	### templates adapted to JCVIs grid
	###	
class GridTask():
	def __init__(self, template="default.q", command = "", name="default", cpu="1", dependson=list(), cwd=".", debug=False):
		
		self.gridjobid=-1
		self.completed=False
		self.queue=template
		self.inputcommand = command
		
		self.cwd=cwd
		self.project = __projectid__
		self.email = __email__
	  
		### remove *e##, *pe## *o## *po##  
		self.retainstreams=" -o /dev/null -e /dev/null "
		
		### debug flag
		self.debugflag = debug
			 
		### the only queue that has more than 4 CPUs...
		if int(cpu)>4:
			self.queue = "himem.q"
				
		if len(dependson)>0:
			holdfor = "-hold_jid "
			for k in dependson:
				holdfor = "%s%s," % (holdfor, k.getJobid())
			holdfor=holdfor.strip(",")		
		else:
			holdfor = ""	
		
		### keep po pe o e streams for debugging purposes
		if self.debugflag:
			self.retainstreams=""
		
			
		### to avoid long command problems, create a script with the command, and invoke that instead of the command directyly.	
		px = "tmp.%s.%s.%s.%s." % (randrange(1,100),randrange(1,100),randrange(1,100),randrange(1,100))
		sx = ".%s.%s.%s.%s.sh" % (randrange(1,100),randrange(1,100),randrange(1,100),randrange(1,100))
		

		##### to avoid too many opened files OSError
		pool_open_files.acquire()
		### bounded semaphore should limit throttle the files opnening for tasks created around the same time	
		scriptfile, scriptfilepath = tempfile.mkstemp(suffix=sx, prefix=px, dir=self.cwd, text=True)
		os.close(scriptfile)
		self.scriptfilepath = scriptfilepath	
		
		os.chmod(self.scriptfilepath,  0777 )
		input= "%s\n" % (self.inputcommand) 
		scriptfile = open(self.scriptfilepath, "w")
		scriptfile.write(input)
		scriptfile.close()
		pool_open_files.release()
		####
		
		
			
		self.templates=dict()
		self.templates["himem.q"]	  = 'qsub %s -P %s -N jh.%s -cwd -pe threaded %s -l "himem" -M %s -m a  %s "%s" ' % (self.retainstreams, self.project, name, cpu, self.email, holdfor, self.scriptfilepath)	
		self.templates["default.q"]	 = 'qsub %s -P %s -N jd.%s -cwd -pe threaded %s -M %s -m a %s "%s" ' % (self.retainstreams, self.project, name, cpu,  self.email, holdfor,  self.scriptfilepath)
		self.templates["fast.q"]	   = 'qsub %s -P %s -N jf.%s -cwd -pe threaded %s -l "fast" -M %s -m a  %s "%s" ' % (self.retainstreams, self.project, name,cpu,  self.email, holdfor, self.scriptfilepath)	
		self.templates["medium.q"]	 = 'qsub %s -P %s -N jm.%s -cwd -pe threaded %s -l "medium" -M %s -m a  %s "%s" ' % (self.retainstreams, self.project, name, cpu,  self.email, holdfor, self.scriptfilepath)
		self.templates["himemCHEAT"] = 'qsub %s -P %s -N jH.%s -cwd -pe threaded %s -l "himem" -M %s -m a  %s "%s" ' % (self.retainstreams, self.project, name, 1, self.email, holdfor, self.scriptfilepath)	
		self.templates["mpi"]		 = 'qsub %s -P %s -N jP.%s -cwd -pe orte %s -M %s -m a  %s mpirun -np %s "%s" ' % (self.retainstreams, self.project, name, cpu, cpu, self.email, holdfor, self.scriptfilepath )
		self.command = ""
		QS.register(self);
						
	def submit(self):

		if not self.queue in self.templates.keys():
			self.queue = QS.pickQ()
		self.command = self.templates[self.queue]
		
		p = Popen(shlex.split(self.command), stdout=PIPE, stderr=PIPE, cwd=self.cwd, close_fds=True)
		
		p.wait()
		out, err = p.communicate()
		
		
			
		err = err.strip()
		out = out.strip()
		

		
		if err!="":
			print err
		
		if out.endswith("has been submitted"):
			self.gridjobid = out.split(" ")[2]
		else:
			print ">>>", out
			print "#FAIL"
				
		return (self.getGridId())

	
	def getGridId(self):
		return self.gridjobid	
	
	def getUniqueID(self):
		return "%s_%s_%s" % (id(self), self.cwd, self.inputcommand)

	def setCompleted(self):
		self.completed=True
		try:
			if not self.debugflag:
				os.remove(self.scriptfilepath)
		except OSError, error:
			print( "%s already gone" % self.scriptfilepath)
		QS.flagRemoval(self)

	def isCompleted(self):
		return (self.completed)
	
	def wait(self):
		while not self.isCompleted():
			time.sleep(0.1)
					

	#################################################
	### Iterator over input fasta file.
	### Only reading when requested
	### Useful for very large FASTA files
	### with many sequences
class	FastaParser:
	def	__init__ (self, x):
		self.filename = x
		self.fp = open(x, "r")	
		self.currline = "" 
		self.currentFastaName = ""
		self.currentFastaSequence = ""
		self.lastitem=False	
			
	def	__iter__(self):
		return(self)	
				
		##### 
	def	next(self):
		for self.currline in self.fp:
			if self.currline.startswith(">"):
				self.currline = self.currline[1:]
				if self.currentFastaName == "":
					self.currentFastaName = self.currline
				else:
					otpt = (self.currentFastaName.strip(), self.currentFastaSequence.strip())
					self.currentFastaName = self.currline
					self.currentFastaSequence = ""	
					self.previoustell = self.fp.tell()
					return (otpt)
				
			else:
				self.addSequence(self.currline)	
		
		if not self.lastitem:
			self.lastitem=True			
			return (self.currentFastaName.strip(), self.currentFastaSequence.strip())
		else:
			raise StopIteration	
				   				
	def	addSequence(self, x):
	   		self.currentFastaSequence = "%s%s" % (self.currentFastaSequence, x.strip())			
	   					
	def	__str__():
		return ("reading file: %s" %self.filename)	



	#################################################
	### The mother of all Steps:
	###
class	DefaultStep(Thread):
	def __init__(self):
		
		#### thread init
		Thread.__init__(self)
		self.random = uniform(0, 10000)
		self.name = ("%s[%s]" % (self.name, self.random))

		#### hash of the current step-path (hash digest of previous steps + current inputs + arguments?)
		self.workpathid = ""
		#### path where the step stores its files
		self.stepdir = ""
					
		#### what needs to be completed for this step to proceed
		#### a list of steps
		self.previous = list()
		
		#### mapping type - path for files
		self.inputs = defaultdict(set)
		
		#### mapping type - name for files
		self.outputs = defaultdict(set)
		
		#### mapping arg  val for program's arguments
		self.arguments= dict()
		
		#### ID of the step...
		self.stepname = ""	
		
		#### flag for completion
		self.completed = False
		self.completedpreviously=False
		self.failed = False
		
		#### keep track of time elapsed
		self.starttime = 0
		self.endtime = 0
		
		#### special flag, some steps might not want to delete the inputs (argcheck)
		self.removeinputs = True
		
		####
		
	def	setInputs(self, x):
		for k,v in x.items():
			for elem in v:
				self.inputs[k].add(elem)
	def setArguments(self, x):
		for k,v in x.items():
			if v=="":
				v=" "
			self.arguments[k] = v
	def setPrevious(self, x):
		if not type(x) is list:
			self.previous.append(x)
		else:
			for elem in x:
				self.previous.append(elem)
	def setName(self, x):
		self.stepname=x
	
	def run(self):	
		self.init()

		if self.failed:
			#self.message("Error detected... ")
			BOH.deregister()
			self.completed=True
			
		elif not self.isDone():
			try:
				self.performStep()
				self.finalize()
			except Exception, inst:	
				self.message("...")
				self.message( type(inst))
				self.message( inst)
				BOH.deregister()
				self.completed=True
				self.failed=True
		else:
			self.message("Completed (previously).")	
			BOH.deregister()
			self.completed=True
			self.completedpreviously=True
		
	def	performStep():
		self.message("in a step...")
				
	def	init(self):
	
		redo=False
		### wait for previous steps to finish
		for k in self.previous:
			while not k.isDone():
				#self.message( "waiting" )
				time.sleep(1)
			if k.hasFailed():
				self.failed=True
			redo=redo or (not k.isDonePreviously())	
		
		#self.message("needs a redo %s" % (redo))
		if not self.failed:
			### time stamp
			self.starttime = time.time()
		
			#### hash of the current step-path (hash digest of previous steps + current inputs + arguments?)
			self.workpathid = self.makeWorkPathId()
			####
			
			### output handler
			BOH.register(self.workpathid)
			###
			
			#self.message("Initializing %s %s" % (self.workpathid, self.name))
			
			#### create directories if necessary
			self.stepdir =""
			self.prepareDir(redo=redo)	
			
	def	makeWorkPathId(self):
		tmp = list()
		tmp.append(self.stepname)
		if self.previous!=None:
			for k in self.previous:
				while k.getWorkPathId()==-1:
					time.wait(1)
				tmp.extend([k.getWorkPathId()])
			
		for k,v in self.inputs.items():
			tmp.extend(["%s=%s" % (k, ",".join(v) ) ] )
		
		for k,v in self.arguments.items():
			tmp.extend(["%s=%s" % (k, v) ] )
				
		tmp.sort()	
		tmp = "\n".join(tmp)
		
		workpathid = hashlib.sha224(tmp).hexdigest()[0:5]
		return (workpathid)
				
	def	getWorkPathId(self):	
		return (self.workpathid)
				
	def	prepareDir(self, redo=False):
		### make step's directory
		self.stepdir = "Step_%s_%s" % (self.stepname, self.workpathid)
		
		
		flush_old = False
		try:	
			os.mkdir(self.stepdir)				
		except OSError, error:
			self.message( "Step directory already exists...")
			flush_old=True
		
		
		if redo:
			if flush_old:
				self.message("Updating...")
				k = "rm -r *"
				task = GridTask(template="pick", name="redo_clean", command=k, cpu=1,  cwd = self.stepdir)
				task.wait()
			else:
				###supposedly no old step-data to flush
				pass	
			
		else:	
			### has analysis been done already?
			try:
				self.parseManifest()
				self.completed=True		
				self.completedpreviously=True
				self.message("Using data generated previously...")
	
			except	IOError, inst:
				#self.message("Will make new manifest...")
				pass
			except Exception, inst:	
				self.message("****ERROR***")	
				self.message(type(inst))
				self.message(inst.args)
				self.message(inst)
				self.message("************")					
																	
	def	finalize(self):
		
		if not self.failed:			
			self.categorizeAndTagOutputs()
			self.makeManifest()
			
			self.endtime = time.time()
			self.message( "+%s\t[Done]" % (str(datetime.timedelta(seconds=round(self.endtime-self.starttime,0))).rjust(17)) ) 
		else:
			self.endtime = time.time()
			self.message( "+%s\t[Fail]" % (str(datetime.timedelta(seconds=round(self.endtime-self.starttime,0))).rjust(17)) ) 
			
		self.completed=True	
		BOH.deregister()
				
	def	makeManifest(self):
		m = open("%s/%s.manifest" % (self.stepdir, self.workpathid), "w")
				
		for type, files in self.inputs.items():
			if len(files)>0:
				m.write("input\t%s\t%s\n" % (type, ",".join(files)) )
				
		for arg, val in self.arguments.items():
			m.write("argument\t%s\t%s\n" % (arg, val ) )		
				
		for type, files in self.outputs.items():
			if len(files)>0:
				m.write("output\t%s\t%s\n" % (type, ",".join(files)) )	
		m.close()
	
	def	determineType(self, filename):
		filename = filename.strip().split(".")
		extension = filename[-1]
		preextension = filename[-2]
		if preextension == "scrap":
			return "scrap"
			
		elif preextension == "align" and extension == "report":
			return "alignreport"
	
		elif extension == "dist" and preextension == "phylip":
			return "phylip"			
		elif extension == "dist":
			return "column"
				
		elif preextension == "tax" and extension =="summary":
			return "taxsummary"
		
		elif preextension == "cdhit" and extension =="clstr":
			return "cdhitclstr"
		elif preextension == "bak" and extension =="clstr":
			return "cdhitbak"
		elif extension == "cdhit":
			return "fasta"
		
		elif extension in ["align", "fna", "fa", "seq", "aln"]:
			return "fasta"
			
		elif extension == "qual":
			return "qfile"
		
		elif extension == "tax":
			return "taxonomy"
		
		elif extension == "names":
			return "name"
		elif extension == "groups":
			return "group"
		elif extension == "files":
			return "file"
		
		elif extension in ["tre", "tree", "dnd"]:
			return "tre"
		
		### sge job files
		elif re.search("po\d{3}", extension) != None:
			return "po"
		elif re.search("pe\d{3}", extension) != None:
			return "pe"		
		elif re.search("o\d{3}", extension) != None:
			return "o"
		elif re.search("e\d{3}", extension) != None:
			return "e"	
		else:
			return extension
		
	def	categorizeAndTagOutputs(self):
		inputs = [x.split("/")[-1] for x in unlist( self.inputs.values()) ]
		for file in glob.glob("%s/*" % self.stepdir):
			file = file.split("/")[-1]
			if file in inputs:
				if self.removeinputs:
					command = "unlink %s" % (file)
					p = Popen(shlex.split(command), stdout=PIPE, stderr=PIPE, cwd=self.stepdir, close_fds=True)	
					out,err = p.communicate()
					p.wait()

				else:
					self.message("kept %s" % file)
					#pass
			elif not file.endswith("manifest"):
				#self.message( "output: %s" % (file))
				
				### make sure that every output file except for the manifest starts with the workpathID
				file = file.split(".")
				if len(file[0]) == len(self.workpathid):
					newfilename = "%s.%s"  % (self.workpathid, ".".join(file[1:]))
				else:
					newfilename = "%s.%s"  % (self.workpathid, ".".join(file[0:])) 
				
				
				if ".".join(file) != newfilename:
					k="mv %s %s" % (".".join(file), newfilename)
					p = Popen(shlex.split(k), stdout=PIPE, stderr=PIPE, cwd=self.stepdir, close_fds=True)	
					out,err = p.communicate()
					p.wait()

				self.outputs[self.determineType(newfilename)].add(newfilename)
						
	def	find(self, arg, ln=True, original=False):
		files=list()	
		if not original:		
			if len(self.inputs[arg])==0:
				tmp = {arg: self.getOutputs(arg)}
				self.setInputs(tmp)	
		else:
			tmp = {arg: self.getOriginal(arg)}		
			self.setInputs(tmp)	

		files = self.inputs[arg]
		
		toreturn=list()	
		
		for file in files:
			if self.isVar(file):
				toreturn.append(file[5:])
			else:	
				tmp = file.strip().split("/")[-1] 
				if (ln):
					command = "cp -s %s %s" % (file, tmp )
				else:
					command = "cp %s %s" % (file, tmp )
				p = Popen(shlex.split(command), stdout=PIPE, stderr=PIPE, cwd=self.stepdir, close_fds=True )	
				out,err = p.communicate()
				p.wait()
				toreturn.append(tmp)		
		#unique	
		toreturn = set(toreturn)
		return list(toreturn)
				
	def	isVar(self,x):
		return x.startswith("[var]")
	
	def	getOutputs(self, arg):
		if self.outputs.has_key(arg):
			otpt = list()
			for x in unlist(self.outputs[arg]):
				if self.isVar(x):
					otpt.append(x)
				else:
					otpt.append("../%s/%s" % (self.stepdir, x))
			return otpt
			
		elif self.previous!=None:
			otpt = list()
			for k in self.previous:
				otpt.extend(k.getOutputs(arg))
			return otpt			
		else:
			return list()
		
	def	getOriginal(self, arg):
		if self.previous == None:
			return self.getOutputs(arg)
		else:	
			current = self.getOutputs(arg)
			otpt = list()			
			for k in self.previous:
				otpt.extend(k.getOriginal(arg))
			if len(otpt)>0:
				return otpt
			else:
				return current		
			
	def	parseManifest(self):
		fp = open("%s/%s.manifest" % (self.stepdir, self.workpathid), "r")
		lines=fp.readlines()
		fp.close()
		for line in lines:
			line = line.strip("\n").split("\t")
			if line[0] == "output":
				type = line[1]
				files = line[2].split(",")
				for file in files:
					self.outputs[type].add(file)
			elif line[0] == "input":
				type = line[1] 
				files = line[2].split(",")
				for file in files:
					self.inputs[type].add(file)	
			elif line[0] == "argument":
				if len(line)==2:
					self.arguments[line[1]] = " "
				else:
					self.arguments[line[1]]=line[2]
				
	def	message(self, text):
			if type(text) == list:
				for line in text: 
					BOH.toPrint(self.workpathid, self.stepname, line)
			else:	
				BOH.toPrint(self.workpathid, self.stepname, text)		
				
	def	isDone(self):
		return self.completed 
	
	def	isDonePreviously(self):
		return self.completedpreviously 

	def	hasFailed(self):
		return self.failed
		
	def	getInputValue(self, arg):
		if self.arguments.has_key(arg):
			return self.arguments[arg]
		else:
			return None
			
	def	setOutputValue(self, arg, val):
		self.outputs[arg] = ["[var]%s" % (val)] 
				
	def	__str__(self):
		otpt = "%s\t%s" % (self.stepname, self.name)
			
		for val in self.previous:
			otpt += "%s\n%s" % (otpt, val.__str__())	
			
		#otpt = "\n".join(set(otpt.strip().split("\n")))		
		return otpt	


	
class	FileImport(DefaultStep):	
	def	__init__(self, INS):
		DefaultStep.__init__(self)
		self.setInputs(INS)
		#self.setArguments(ARGS)
		#self.setPrevious(PREV)
		self.setName("FILE_input")
	 	self.start()
			
	def	performStep(self):
		for type in self.inputs.keys():
			files = self.inputs[type]
			for file in files:
				file = file.split("~")
				if len(file)>1:
					file, newname = file
					tmp = file.strip().split("/")[-1]
					k = "cp %s %s.%s" % (file, newname, type)
					
				else:
					file = file[0]
					tmp = file.strip().split("/")[-1]
					k ="cp %s imported.%s"	% (file, tmp)
								
				p = Popen(shlex.split(k), stdout=PIPE, stderr=PIPE, cwd=self.stepdir, close_fds=True)
				self.message(k)
				out,err = p.communicate()
				p.wait()
#				
class	ArgumentCheck(DefaultStep):
	def	__init__(self, SHOW, PREV):
		ARGS = {"show":SHOW}
		DefaultStep.__init__(self)
		#self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("ArgCheck")
		#self.nodeCPUs=nodeCPUs
		self.removeinputs=False
	 	self.start()
		
	def	performStep(self):
		x = self.getInputValue("show")
		if x!=None:
			for type in x.split(","):
				for file in self.find(type):
					self.message("%s: %s" % (type,file))	

class	OutputStep(DefaultStep):
	def	__init__(self, NAME, SHOW, PREV):
		ARGS = {"show":SHOW}
		DefaultStep.__init__(self)
		#self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("OUTPUT_%s" % (NAME))
		#self.nodeCPUs=nodeCPUs
		self.removeinputs=False
	 	self.start()
		
	def	performStep(self):
		x = self.getInputValue("show")
		if x!=None:
			for type in x.split(","):
				for file in self.find(type.strip(), ln = False):
					self.message("%s: %s" % (type,file))
		
class	SFFInfoStep(DefaultStep):
	def	__init__(self, INS, ARGS, PREV):
		DefaultStep.__init__(self)
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("sffinfo")
		self.start()
		
	def	performStep(self):
		steps = list()		
		for sff in self.find("sff"):
			
			k = "/usr/local/bin/sffinfo -s %s > %s.fasta" % (sff, sff)
			self.message(k)
			task = GridTask(template="pick", name=self.stepname, command=k, cpu=1, dependson=list(), cwd = self.stepdir)
			steps.append(task)
			
			k = "/usr/local/bin/sffinfo -q %s > %s.qual" % (sff, sff)
			self.message(k)
			task = GridTask(template="pick", name=self.stepname, command=k, cpu=1, dependson=list(), cwd = self.stepdir)
			steps.append(task)
			
			k = "/usr/local/bin/sffinfo -f %s > %s.flow" % (sff, sff)
			self.message(k)
			task = GridTask(template="pick", name=self.stepname, command=k, cpu=1, dependson=list(), cwd = self.stepdir)
			steps.append(task)
			
		for s in steps:
			s.wait()
	
class	MothurStep(DefaultStep):
	def __init__(self, NM, nodeCPUs, INS, ARGS, PREV):
	 	DefaultStep.__init__(self)
	 	self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName(NM)
		self.nodeCPUs=nodeCPUs
		self.start()
	 			 
	def	makeCall(self):
		FORCE = self.getInputValue("force")
		x = MOTHUR.getCommandInfo(self.stepname)
		#self.message(self.inputs)
		if FORCE != None :
			TYPES = FORCE.strip().split(",")
		else:
			TYPES = x.getInputs()
			
		mothurargs = list()
		
		for TYPE in TYPES:
			#### on occasion, mothur overwrites the original file - namely names file
			#### FALSE creates a copy
			#### TRUE creates a link
			if TYPE=="name":
				tmp = self.find(TYPE, False)
			else:
				tmp = self.find(TYPE, True)
				
			if len(tmp)>0:
				mothurargs.append ("%s=%s" % (TYPE, "-".join(tmp)))
			else:
				if x.isRequired(TYPE):
					self.message("Required argument '%s' not found!" % (TYPE))	
					raise Exception 
				else:
					self.message("Optional argument '%s' not found, skipping" % (TYPE))	
		
		for arg, val in self.arguments.items():
		
			if x.isAnArgument(arg):
				mothurargs.append("%s=%s" % (arg, val))
			elif arg=="find":
				for a in val.strip().split(","):
					self.message(a)
					valstoinsert = self.find(a)
					self.message(valstoinsert)
					if len(valstoinsert)>0:
						mothurargs.append("%s=%s" % (a, "-".join(valstoinsert)))
					else:
						self.message("skipping '%s' - not found" % (a))
			else:
				self.message("skipping '%s', as it is not an argument for %s" % (arg, self.stepname))
	
		### method is parallelizable,   
		if x.isAnArgument("processors") and "processors" not in self.arguments.keys():
			mothurargs.append("%s=%s" % ("processors", self.nodeCPUs ))
			self.message("Will run on %s processors" % (self.nodeCPUs))
		
		himemflag=False
		### steps requiring lots of memory
		if self.stepname in ("clearcut", "align.seq"):
			himemflag=True
			self.message("Needs lots of memory")
		
		command = "%s(%s)" % (self.stepname, ", ".join(mothurargs))
		
		return (command, x.isAnArgument("processors"), himemflag)
	
	def	performStep(self):
		call, parallel, himem = self.makeCall()	
		k = "%smothur \"#%s\"" % (mothurpath, call)
		self.message(k)
		if (parallel and self.nodeCPUs>1):
			task = GridTask(template=defaulttemplate, name=self.stepname, command=k, cpu=self.nodeCPUs, dependson=list(), cwd = self.stepdir)
		elif (himem):
			task = GridTask(template="himem.q", name=self.stepname, command=k, cpu=self.nodeCPUs, dependson=list(), cwd = self.stepdir)
		else:
			task = GridTask(template="pick", name=self.stepname, command=k, cpu=1, dependson=list(), cwd = self.stepdir)
		
		task.wait()
		self.parseLogfile()
		
	def	parseLogfile(self):
		for f in glob.glob("%s/*.logfile" % (self.stepdir)):
			line = ""
			for line in loadLines(f):
				### UCHIME throws an error when it does not find chimeras, even though it completes.
				if line.find ("ERROR")>-1 and line.find("uchime")==-1:
					self.failed=True
			
			### last line
			if line.find("quit()")==-1:
				self.failed=True

class	MothurSHHH(DefaultStep):
	def __init__(self,  PREV):		
		DefaultStep.__init__(self)
		#self.setInputs(INS)
		#self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("MPyro")
		#self.nodeCPUs=nodeCPUs
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		
		TOC = self.find("file")
		flows = self.find("flow")
		
		TOC = loadLines("%s/%s" % (self.stepdir, TOC[0]))
		TOC = [ ".".join(x.strip().split(".")[1:]) for x in TOC]
			
		for f in flows:
			tmp = ".".join(f.split(".")[1:])
			
			if tmp in TOC:
				
				### split tmp into 10,000 lines chunks
				k = "split -l 7000 -a 3 %s %s.split." % (f, f)
				task = GridTask(template="pick", name="MPyroSplit", command=k, cpu=1,  cwd = self.stepdir)
				tasks.append(task)				
			else:
				self.message("skipping %s" % (f))  
		
		self.message("splitting %s file(s)" % len(tasks))
		
		for task in tasks:
			task.wait()
			time.sleep(1)	
		
		################################################
		tasks = list()
		
		for chunk in glob.glob("%s/*.split.*" % (self.stepdir)):
			chunk = chunk.split("/")[-1]
			#self.message(chunk)	
			call = "shhh.flows(flow=%s)" % (chunk)
			k = "%smothur \\\"#%s\\\"" % (mothurpath, call)
			task = GridTask(template="pick", name="Mpyro", command=k, cpu=1,  cwd = self.stepdir)
			tasks.append(task)
				
		self.message("processing %s file(s)" % len(tasks))
		
		for task in tasks:
			task.wait()
			time.sleep(1)		
					
class	LUCYcheck(DefaultStep):
	def __init__(self, nodeCPUs, PREV):
		DefaultStep.__init__(self)
		self.nodeCPUs=nodeCPUs
		#self.setInputs(INS)
		#self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("LUCY_check")
		self.nodeCPUs=nodeCPUs
		if self.nodeCPUs>32:
			self.nodeCPUs=30
	 	self.start()
				
	def	performStep(self):
		f = self.find("fasta")[0]
		q = self.find("qfile")[0]
		k ="lucy -error 0.002 0.002 -bracket 20 0.002 -debug -xtra %s -output %s.fastalucy %s.qfilelucy %s %s" % (self.nodeCPUs, f,q, f,q)
		self.message(k)
		if self.nodeCPUs>2:
			task = GridTask(template=defaulttemplate, name=self.stepname, command=k, cpu=self.nodeCPUs,  dependson=list(), cwd = self.stepdir)
		else:
			task = GridTask(template="pick", name=self.stepname, command=k, cpu=self.nodeCPUs,  dependson=list(), cwd = self.stepdir)
		task.wait()	
		
class	LUCYtrim(DefaultStep):	
	def __init__(self, PREV):
		DefaultStep.__init__(self)
		#self.setInputs(INS)
		#self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("LUCY_trim")
		#self.nodeCPUs=nodeCPUs
	 	self.start()
		
				
	def	performStep(self):
		f = self.find("fastalucy")[0]
		q = self.find("qfilelucy")[0]
		
		k = "python %s/fastAtrimmer.py -l %s %s %s " % (scriptspath, f.split(".")[0], f, q)
		self.message(k)
		task = GridTask(template="pick", name=self.stepname, command=k, cpu=1,  cwd = self.stepdir)
		task.wait()

class	MatchGroupsToFasta(DefaultStep):
	def __init__(self, INS, PREV):		
		DefaultStep.__init__(self)
		self.setInputs(INS)
		#self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("MatchGroups")
		#self.nodeCPUs=nodeCPUs
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		f = self.find("fasta")
		f = f[0]
		g = self.find("group")
		g = g[0]
		
		n = self.find("name")
		if len(n)>0:
			n = "-n %s" % (n[0])
		else:
			n = ""
		
		k = "python %s/MatchGroupsToFasta.py %s -f %s -g %s -o %s.matched.group" % (scriptspath, n, f, g, ".".join(g.split(".")[:-1]))
		self.message(k)	
		task = GridTask(template="pick", name=self.stepname, command=k, cpu=1,  cwd = self.stepdir)
		task.wait()	

class	MatchGroupsToList(DefaultStep):
	def __init__(self, INS, PREV):		
		DefaultStep.__init__(self)
		self.setInputs(INS)
		#self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("MatchGroups")
		#self.nodeCPUs=nodeCPUs
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		f = self.find("list")
		f = f[0]
		g = self.find("group")
		g = g[0]

		
		k = "python %s/MatchGroupsToFasta.py -l %s -g %s -o %s.matched.group" % (scriptspath, f, g, ".".join(g.split(".")[:-1]))
		self.message(k)	
		task = GridTask(template="pick", name=self.stepname, command=k, cpu=1,  cwd = self.stepdir)
		task.wait()			
			
class	FileMerger(DefaultStep):
	def __init__(self, TYPES, PREV, prefix="files"):
		ARGS = 	{"types": TYPES}		
		DefaultStep.__init__(self)
		#self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("FILE_cat")
		#self.nodeCPUs=nodeCPUs
		self.prefix = prefix
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		for t in self.getInputValue("types").strip().split(","):
			files = self.find(t)
			if len(files)>0 and len(files)<25:
				k = "cat %s > %s_x%s.merged.%s" % (" ".join(files), self.prefix, len(files), t)
				self.message(k)	
				task = GridTask(template="pick", name="cat", command=k, cpu=1,  cwd = self.stepdir)
				tasks.append(task)
			elif len(files)>=25:
				k = "cat *.%s* > %s_x%s.merged.%s" % (t, self.prefix, len(files), t)
				self.message(k)	
				task = GridTask(template="pick", name="cat", command=k, cpu=1,  cwd = self.stepdir)
				tasks.append(task)
			#else:
			#	self.failed=True
			
		for task in tasks:
			task.wait()
			time.sleep(1)

class	FileSort(DefaultStep):
	def __init__(self, TYPES, PREV):
		ARGS = 	{"types": TYPES}		
		DefaultStep.__init__(self)
		#self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("FILE_sort")
		#self.nodeCPUs=nodeCPUs
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		for t in self.getInputValue("types").strip().split(","):
			files = self.find(t)
			if len(files)>0:
				k = "sort -n %s > files_x%s.sorted.%s" % (" ".join(files), len(files), t)
				self.message(k)	
				task = GridTask(template="pick", name="sort", command=k, cpu=1,  cwd = self.stepdir)
				tasks.append(task)
		for task in tasks:
			task.wait()
			time.sleep(1)
			
class	FileSplit(DefaultStep):
	def __init__(self,  TYPES, PREV, chunk=40000):		
		DefaultStep.__init__(self)
		ARGS = 	{"types": TYPES}
		#self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("FILE_split")
		#self.nodeCPUs=nodeCPUs
		self.chunk = chunk
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		
		for t in self.getInputValue("types").strip().split(","):
			files = self.find(t)
			for f in files:
				k = "split -a 5 -l %s %s %s.split. " % (self.chunk, f, f)
				self.message(k)	
				task = GridTask(template="himem.q", name="filesplit", command=k, cpu=1,  cwd = self.stepdir)
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
			
class	FileType(DefaultStep):
	def __init__(self, ARGS, PREV):		
		DefaultStep.__init__(self)
		#self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("FILE_type")
		#self.nodeCPUs=nodeCPUs
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		for input, output in self.arguments.items():
			files = self.find(input)
			for file in files:
				outname = "%s.%s" % (file, output)
				k = "cp %s %s" % (file, outname)
				self.message(k)	
				task = GridTask(template="pick", name="sort", command=k, cpu=1,  cwd = self.stepdir)
				tasks.append(task)
		for task in tasks:
			task.wait()
			time.sleep(1)			
			
class	CleanFasta(DefaultStep):
	def __init__(self, INS, PREV):		
		DefaultStep.__init__(self)
		self.setInputs(INS)
		#self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("CleanFasta")
		#self.nodeCPUs=nodeCPUs
	 	self.start()
		
	def	performStep(self):
		tasks = list()
		f = self.find("fasta")
		f = f[0]
		k = "python %s/CleanFasta.py -i %s -o %s.dash_stripped.fasta" % (scriptspath,f, ".".join(f.split(".")[:-1]))
		self.message(k)	
		task = GridTask(template="pick", name=self.stepname, command=k, cpu=1,  cwd = self.stepdir)
		task.wait()			
		
class	MakeNamesFile(DefaultStep):
	def __init__(self, PREV):
		DefaultStep.__init__(self)
		#self.setInputs(INS)
		#self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("FILE_names")
		#self.nodeCPUs=nodeCPUs
	 	self.start()
		
	def	performStep(self):
		f = self.find("fasta")[0]
		self.message("Creating 'names' file for sequences in {0}".format( f))

		newname = f.strip().split(".")[:-1]
		newname = "%s.names" % (".".join(newname)) 
		otpt = open("%s/%s" % (self.stepdir,newname ), 'w')
		for head, seq in FastaParser("%s/%s" % (self.stepdir, f)):
			head = head.strip().split()[0]
			otpt.write("%s\t%s\n" % (head, head))
		otpt.close()	
			
class	MakeGroupsFile(DefaultStep):
	def __init__(self, PREV, id):
		ARGS = 	{"groupid": id}
		DefaultStep.__init__(self)
		#self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("FILE_groups")
		#self.nodeCPUs=nodeCPUs
	 	self.start()
		
	def	performStep(self):
		f = self.find("fasta")[0]
		id = self.getInputValue("groupid")
		self.message("Creating 'groups' file; '{0}' for sequences in {1}".format(id, f))
		newname = f.strip().split(".")[:-1]
		newname = "%s.groups" % (".".join(newname)) 
		otpt = open("%s/%s" % (self.stepdir, newname ), 'w')
		for head, seq in FastaParser("%s/%s" % (self.stepdir, f)):
			head = head.strip().split()[0]
			otpt.write("%s\t%s\n" % (head, id))
		otpt.close()			

class	MakeQualFile(DefaultStep):
	def __init__(self, PREV, q):
		ARGS = 	{"qual": q}
		DefaultStep.__init__(self)
		#self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("FILE_qfile")
		#self.nodeCPUs=nodeCPUs
	 	self.start()
		
		
	def	performStep(self):
		f = self.find("fasta")[0]
		q = self.getInputValue("qual")
		self.message("Creating 'qual' file; '{0}' for sequences in {1}".format(q, f))
		newname = f.strip().split(".")[:-1]
		newname = "%s.qual" % (".".join(newname)) 
		otpt = open("%s/%s" % (self.stepdir, newname ), 'w')
		for head, seq in FastaParser("%s/%s" % (self.stepdir, f)):
			otpt.write(">%s\n" % (head))
			for k in seq:
				otpt.write("%s " % (q))
			otpt.write("\n")	
		otpt.close()			
		
class 	AlignmentSummary(DefaultStep):
	def __init__(self, ARGS, PREV):
		DefaultStep.__init__(self)
		#self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("AlignmentSummary")
		#self.nodeCPUs=nodeCPUs
	 	self.start()
		
	def	performStep(self):
		self.project = __projectid__
		self.mailaddress = __email__
			
		f = self.find("fasta")[0] 
		
		ref =  self.getInputValue("ref")
		if ref == None:
			ref="e_coli2"
			
		self.message("summarizing an alignment in %s" % (f) )
		k = "python %s/summarizeAlignment.py -P %s -M %s -t 500 -p %s -i %s -o %s.alsum" % (scriptspath, self.project, self.mailaddress, ref, f,f)
		self.message(k)
		task = GridTask(template="pick", name=self.stepname, command=k, cwd = self.stepdir)
		task.wait() 	

class	AlignmentPlot(DefaultStep):
	def __init__(self, ARGS, PREV):
		DefaultStep.__init__(self)
		#self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("AlignmentPlot")
		#self.nodeCPUs=nodeCPUs
	 	self.start()
		
	def	performStep(self):
		f = self.find("alsum")[0]
		
		ref =  self.getInputValue("ref")
		if ref == None:
			ref="e_coli"
			
		trimstart = self.getInputValue("trimstart")
		if trimstart==None:
			trimstart=0
			
		trimend = self.getInputValue("trimend")
		if trimend ==None:
			trimend=0
		
		tmp = open("%s/alsum.r" % (self.stepdir), "w")
		tmp.write("source(\"%s/alignmentSummary.R\")\n" % (scriptspath))	
		tmp.write("batch2(\"%s\", ref=\"%s\", trimstart=%s, trimend=%s )\n" % (f, ref, trimstart, trimend))
		tmp.close()
		k = "R CMD BATCH alsum.r"
		task = GridTask(template="pick", name=self.stepname, command=k, cwd = self.stepdir)
		task.wait()
		
class	GroupRetriever(DefaultStep):
	def __init__(self, PREV):
		DefaultStep.__init__(self)
		#self.setInputs(INS)
		#self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("Retrieve_Groups")
	 	self.start()
	 	
	def performStep(self):
		group = self.find("group")[0]
		groups = set()
		for line in loadLines("%s/%s" % (self.stepdir, group)):
			groups.add(line.strip().split("\t")[1])
		groups  = "-".join(groups)
		self.message(groups)
		self.setOutputValue("groups", groups)
	 	
class	CDHIT_454(DefaultStep):
	def __init__(self, nodeCPUs, ARGS, PREV):
		DefaultStep.__init__(self)
		if ARGS.has_key("T"):
			self.nodeCPUs = ARGS["T"]
		else:
			self.nodeCPUs=nodeCPUs
			ARGS["T"]=self.nodeCPUs		

		#self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("CDHIT_454")
		self.start()
						
	def	performStep(self):

		args = ""	
				
		for arg, val in self.arguments.items():
			args = "%s -%s %s" % (args, arg, val) 

		fs = self.find("fasta")
		fs.extend(self.find("mate1"))
		fs.extend(self.find("mate2"))
		tasks=list()
		for f in fs:
			k ="%scd-hit-454 -i %s -o %s.cdhit %s" % (cdhitpath, f, f, args)
			self.message(k)
			if self.nodeCPUs>2:
				task = GridTask(template=defaulttemplate, name=self.stepname, command=k, cpu=self.nodeCPUs,  dependson=list(), cwd = self.stepdir, debug=True)
			else:
				task = GridTask(template="himem.q", name=self.stepname, command=k, cpu=self.nodeCPUs,  dependson=list(), cwd = self.stepdir, debug=True)
			tasks.append(task)
		
		for task in tasks:
			task.wait()	

class	CDHIT_Mothurize(DefaultStep):
	def __init__(self, ARGS, PREV):
		DefaultStep.__init__(self)
		#self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("CDHIT_Mothurise")
		self.start()
		
	def	performStep(self):
	
		### mode can be "name" or None (to produce a name file)
		### mode can be "#.## (to produce a labeled ra/sa/list combo)
		m = self.getInputValue("mode")
		if m == None:
			m = "name"
		
		modeswitch = "-o %s" % (m)
		
		### is there an optional names file?
		n = self.find("name")	
		if len(n)>0:
			n = n[0]
			nameswitch = "-n %s" % (n)
		else:
			nameswitch = "" 
	
		### is there a required cluster file
		clst = self.find("cdhitclstr")
		
		if len(clst)>0:
			k = "python %sCDHIT_mothurize_clstr.py -c %s  %s %s" % (scriptspath, clst[0], nameswitch, modeswitch)	
			self.message(k)
			task = GridTask(template="pick", name=self.stepname, command=k, dependson=list(), cwd = self.stepdir)
			task.wait()	
			
		else:
			self.failed=True
			
class	CDHIT_EST(DefaultStep):
	def __init__(self, nodeCPUs, ARGS, PREV):
		DefaultStep.__init__(self)
		
		if ARGS.has_key("T"):
			self.nodeCPUs = ARGS["T"]
		else:
			self.nodeCPUs=nodeCPUs
			ARGS["T"]=self.nodeCPUs		
		
		#self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("CDHIT_EST")
		self.start()
				
	def	performStep(self):
		f = self.find("fasta")[0]
		args = ""
		dist = 1
		for arg, val in self.arguments.items():
			args = "%s -%s %s" % (args, arg, val) 
			if arg == "c":
				dist = dist - (float(val)) 
		
		k ="%scd-hit-est -i  %s -o %s._%s_.cdhit %s" % (cdhitpath, f, f, dist, args)
		
		self.message(k)
		if self.nodeCPUs>2:
			task = GridTask(template=defaulttemplate, name=self.stepname, command=k, cpu=self.nodeCPUs,  dependson=list(), cwd = self.stepdir, debug=True)
		else:
			task = GridTask(template="himem.q", name=self.stepname, command=k, cpu=self.nodeCPUs,  dependson=list(), cwd = self.stepdir, debug=True)
		task.wait()	

class	R_defaultplots(DefaultStep):
	def __init__(self, INS, ARGS, PREV):
		DefaultStep.__init__(self)	
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("R_plots")
		self.start()
	def	performStep(self):
		f = self.find("taxsummary")
		anno = self.find("annotation")[0]
		
		tasks = list()
	
		script = open("%s/script.r" % (self.stepdir), "w")
		script.write("""source("%sConsTaxonomyPlots.R")\n""" % (scriptspath))
		
		for file in f:
			dist = ".%s"% (self.getInputValue("dist"))
			
			if file.find(dist)>-1 and file.find("seq")>-1 :
				script.write("""makeDefaultBatchOfPlots("%s", "%s", fileprefix="SEQnum")\n""" % (anno, file))
	
				
			elif file.find(dist)>-1 and file.find("otu")>-1 :
				script.write("""makeDefaultBatchOfPlots("%s", "%s", fileprefix="OTUnum")\n""" % (anno, file))
			
		script.close()
		
		k =	"R CMD BATCH script.r"	
		
		self.message(k)
		task = GridTask(template="pick", name=self.stepname, command=k, dependson=list(), cwd = self.stepdir)
		task.wait()
				
class	R_OTUplots(DefaultStep):
	def __init__(self, INS, ARGS, PREV):
		DefaultStep.__init__(self)	
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("R_plots_otu")
		self.start()
	def	performStep(self):
	
		####OTUS
		f = self.find("fasta")
		tasks = list()

		#script = open("%s/script.r" % (self.stepdir), "w")
		#script.write("""source("%sOtuReadPlots.r")\n""" % (scriptspath))
		
		for file in f:
			if file.find("annotated.fasta")>0:
				k = """grep ">" %s | awk '{FS = "|"; OFS="\t"} {print $4, $5}' > %s.otustats""" % (file, file)
				task = GridTask(template="pick", name=self.stepname, command=k, dependson=list(), cwd = self.stepdir, debug=False)
				tasks.append(task)
				#script.write("""makeBatch("%s.otustats")\n""" % (file))
				
		####COVERAGE
		f = self.find("clcassemblystats")
		
		#for file in f:
				#script.write("""makeBatchCoverage("%s")\n""" % (file))		
		
		#script.close()
		
		### make sure all conversions are complete
		for task in tasks:
			task.wait()
		
		k =	"R CMD BATCH %sOtuReadPlots.r" % (scriptspath)
		self.message(k)
		
		task = GridTask(template="pick", name=self.stepname, command=k, dependson=list(), cwd = self.stepdir)
		task.wait()				
				
class	R_rarefactions(DefaultStep):
	def __init__(self, INS, ARGS, PREV):
		DefaultStep.__init__(self)	
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("R_rarefactions")
		self.start()
	def	performStep(self):
		for k in "r_nseqs,rarefaction,r_simpson,r_invsimpson,r_chao,r_shannon,r_shannoneven,r_coverage".strip().split(","):
			f = self.find(k)
		k =	"R CMD BATCH %srarefactions.R" % (scriptspath)	
		self.message(k)
		task = GridTask(template="pick", name=self.stepname, command=k, dependson=list(), cwd = self.stepdir)
		task.wait()					

class	AlignmentTrim(DefaultStep):
	def __init__(self, INS, ARGS, PREV):
		DefaultStep.__init__(self)
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("AlTrim")
		self.start()
						
	def	performStep(self):
		f = self.find("fasta")[0]
		args = ""					
		for arg, val in self.arguments.items():
			args = "%s -%s %s" % (args, arg, val) 
		
		k ="python %salignmentTrimmer.py %s -I %s" % (scriptspath, args, f)
		self.message(k)
		task = GridTask(template="pick", name=self.stepname, command=k, cpu=1,  dependson=list(), cwd = self.stepdir)
		task.wait()	

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
					
class	AnnotateClusters(DefaultStep):
	def __init__(self, INS, ARGS, PREV):
		DefaultStep.__init__(self)
		self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("Annotate")
		self.start()
						
	def	performStep(self):
		l = self.find("list")
		t = self.find("taxonomy")
		f = self.find("fasta")
		g = self.find("group")
		
		self.message(l)
		self.message(t)
		self.message(f)
		self.message(g)
			
		if len(l)==0 or len(t)==0 or len(f)==0 or len(g) == 0:
			self.failed=True
		else:
		
			tasks=list()	
			for fasta in f:
				dist = fasta.split("_")[-2]
				
				for tax in t:
					if tax.find(dist)>-1 and tax.find("otu")==-1:
						k = "python %sRetrieve.py %s %s %s %s %s" % (scriptspath, dist, l[0], tax, g[0], fasta)
						self.message(k)
						task = GridTask(template="pick", name=self.stepname, command=k, cpu=1,  dependson=list(), cwd = self.stepdir)
						tasks.append(task)
			
			for task in tasks:
				task.wait()	
		#self.failed=True

######		
######	these are mostly for contaminant checking/metagenomcis		
###########################################################							
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
		
		### index a filename using the filename sans the step id (0) and extension (-1)
		filters = {".".join(key.strip().split(".")[1:-1]) : key for key in filter}
		
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
				k = "python %sMateFilter.py %s -i %s -f %s " % (scriptspath, argstring, file, filters[name])
				if len(m1)==1:
					self.message(k)
			
				task = GridTask(template="pick", name="%s" % (self.stepname), command=k, cpu=1,  cwd = self.stepdir)
				tasks.append(task)
			else:
				missing +=1	
			
			if missing>0:
				self.message("%s missing filters observed..." % (missing))
				self.failed = True
			
		for task in tasks:
			task.wait()	

################################################
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
		
		m1 = self.find("mate1")
		m1.extend(self.find("mate2"))
		m1.extend(self.find("fastq"))
				
		argstring = ""
		for arg, val in self.arguments.items():
			argstring = "%s %s %s " % (argstring, arg, val) 

		tasks = list()
		self.message("processing %s files..." % len(m1))
		for file in m1:
			suffix = file.split(".")[-1]
			prefix = ".".join(file.split(".")[:-1])
			if suffix=="fastq":
				k = "fastq_to_fasta %s -i %s -o %s.fasta" % (argstring, file, prefix)
			else:
				k = "fastq_to_fasta %s -i %s -o %s.fasta.%s" % (argstring, file, prefix, suffix)
			if len(m1)<10:
				self.message(k)
			task = GridTask(template="pick", name="%s" % (self.stepname), command=k, cpu=1,  cwd = self.stepdir)
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
				
				k = "python %sinterweaveMates.py %s" % (scriptspath, argstring)
				if len(m1)<10:
					self.message(k)
				task = GridTask(template="pick", name="%s" % (self.stepname), command=k, cpu=1,  cwd = self.stepdir)
				tasks.append(task)
			
		for task in tasks:	
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
		cpus=1
		template="pick"
		if len(x)!=1:
			self.failed=True
		else:
			m1 = self.find("fasta")[0]			
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
			contigs = list()
			reads = list()
			if len(fastas)<=2:
				for f in fastas:
					if f.find("contig")>-1:
						contigs.append(f)
					else:
						reads.append(f)
				
						
				contigs = contigs[0]
				reads = reads[0]		
				m1 = self.find("mate1", original=True)[0]	
				m2 = self.find("mate2", original=True)[0]			
							
				prefix = ".".join(contigs.split(".")[:-1])
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
				
				### original / all reads		
#				done = False
#				while not done:			
#					k = "/usr/local/packages/clc-ngs-cell/clc_ref_assemble_long  %s -o %s.original.cas -q -i %s %s -d %s" % (argstring, prefix, m1, m2, contigs)
						

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
		k = "R CMD BATCH %sStatFastaFiles.R" % (scriptspath)
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
			
class	CDHIT_Perls(DefaultStep):
	def __init__(self, ARGS, PREV):
		DefaultStep.__init__(self)
		#self.setInputs(INS)
		self.setArguments(ARGS)
		self.setPrevious(PREV)
		self.setName("CDHITperls")
		#self.nodeCPUs=nodeCPUs
	 	self.start()
		
	def	performStep(self):
		x = self.find("cdhitclstr")
		tasks=list()
		for cluster in x:
			k = "%sclstr2tree.pl %s > %s.tre" % (cdhitpath, cluster, cluster)
			self.message(k)
			task = GridTask(template="pick", name=self.stepname, command=k, cwd = self.stepdir, debug=False)
			tasks.append(task)
			
			k = "%sclstr_size_histogram.pl %s > %s.hist.tab.txt " % (cdhitpath, cluster, cluster)
			self.message(k)
			task = GridTask(template="pick", name=self.stepname, command=k, cwd = self.stepdir, debug=False)
			tasks.append(task)
			
			k = "%sclstr_size_stat.pl %s  > %s.size.tab.txt" % (cdhitpath, cluster, cluster)
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
		k = "python %sFastQEncoding.py %s > %s.offset" % (scriptspath, x[0], x[0])
		self.message(k)
		task = GridTask(template="pick", name=self.stepname, command=k, cwd = self.stepdir, debug=False)
		task.wait()
		
		otpt = ""
		for line in loadLines("%s/%s.offset" % (self.stepdir, x[0])):
			otpt = "%s%s" % (otpt, line.strip())
		
		self.message("%s -> %s" % (x[0], otpt)) 
		self.setOutputValue("-Q", otpt)
		
	
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

def	unlist(struct):
	for x in struct:
		if type(x) is tuple or type(x) is list or type(x) is set :
			for y in unlist(x):
				yield y
		else:
			yield x
	
def	init(id, e):
	global __projectid__
	global __email__
	
	__projectid__ = id
	__email__ = e

def	getNext():
	global __counter__
	__counter__+=1
	
	return __counter__

def	revComp(string):
	global transtab
	string=string.upper()
	#reverse
	string = string [::-1]
	return string.translate(transtab)

def getQ(file):
	k = "python %sFastQEncoding.py %s" % (scriptspath, file)
	p = Popen(shlex.split(k), stdout=PIPE, stderr=PIPE, close_fds=True)
	out,err = p.communicate()
	p.wait()
	return "%s" % (out.strip())
	
#################################################
##		Arguments
##

#################################################
##		Begin
##

from string import maketrans
inttab=  "ACGTN"
outtab = "TGCAN"
transtab = maketrans(inttab, outtab)

pool_open_files = BoundedSemaphore(value=10, verbose=False)

mothurpath  = "/usr/local/devel/ANNOTATION/sszpakow/YAP/bin/mothur-current/"
cdhitpath 	= "/usr/local/devel/ANNOTATION/sszpakow/YAP/bin/cdhit-current/"
scriptspath = "/usr/local/devel/ANNOTATION/sszpakow/YAP/scripts/"
binpath = "/usr/local/devel/ANNOTATION/sszpakow/YAP/bin/"

__counter__ = 0

BOH = BufferedOutputHandler()
MOTHUR = MothurCommandInfo(path=mothurpath)
QS = TaskQueueStatus(update = 0.1, maxnodes=250)

defaulttemplate = "himem.q"




#################################################
##		Finish
#################################################
