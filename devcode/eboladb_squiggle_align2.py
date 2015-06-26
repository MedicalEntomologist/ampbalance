#!/usr/bin/env python
import sklearn.preprocessing
import random
import math
import csv
import numpy as np
import array as ar
import configargparse
import subprocess
import shutil
import glob
import h5py
from itertools import islice
from collections import OrderedDict

global oper



parser = configargparse.ArgParser(description='ampbalance: A program designed to balance amplicons from a specific reference sequence post sequencing on ONT minIONs but prebasecalling. Developed by Matt Loose @mattloose or matt.loose@nottingham.ac.uk for help!')
parser.add('-fasta', '--reference_fasta_file', type=str, dest='fasta', required=True, default=None, help="The fasta format file for the reference sequence for your organism.")
parser.add('-ids', nargs = '*', dest='ids',required=True, help = 'A list of start and stop positions for each amplicon from the reference genome - should be space separated with fasta_name:start-stop.\n e.g.\n EM_079517:27-1938 EM_078517:1927-3828 EM_078517:3823-5718 EM_078517:5759-7633 EM_078517:7601-10007 EM_078517:9550-10921 EM_078517:10944-12354 EM_078517:12354-14252 EM_078517:14253-15680 EM_078517:15691-17087 EM_078517:16632-18553\n ')
parser.add('-w', '--watch-dir', type=str, required=True, default=None, help="The path to the folder containing the downloads directory with fast5 reads to analyse - e.g. C:\data\minion\downloads (for windows).", dest='watchdir')
parser.add('-d', '--depth',type=int, required=True, default=None, help = 'The desired coverage depth for each amplicon. Note this is unlikely to be achieved for each amplicon and should probably be an overestimate of the minimum coverage required.', dest='depth')
args = parser.parse_args()



######################################################
# Connect to the database 							 #
######################################################

### This is only used to update the model file.

dbhost = 'localhost'
dbusername = ''
dbpass = ''
dbport = 3306

limitbases = "50,200"

######################################################
	
######################################################
def get_amplicons():
	for sequence in args.ids:
		print sequence
		start = int(float(sequence.split(':', 1 )[1].split('-',1)[0]))
		stop = int(float(sequence.split(':', 1 )[1].split('-',1)[1]))
		print start
		print stop
		REVERSE_stop = seqlengths['EM_079517']-start
		REVERSE_start = seqlengths['EM_079517']-stop
		print REVERSE_stop
		print REVERSE_start

######################################################
def get_seq_len(ref_fasta):
	seqlens=dict()
	for record in SeqIO.parse(ref_fasta, 'fasta'):
		seq=record.seq
		seqlens[record.id]=len(seq)
	return seqlens

#######################################################################
def runProcess(exe):
	p=subprocess.Popen(exe, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	while(True):
		retcode= p.poll()
		line=p.stdout.readline()
		yield line
		if(retcode is not None):
			break
	
#######################################################################
def squiggle_search2(squiggle,kmerhash2,channel_id,read_id,seqlen):
	result=[]
	for id in kmerhash2:
		for ref in kmerhash2[id]:
			#print len(kmerhash2[id][ref]['F'])
			queryfile=str(channel_id)+"_"+str(read_id)+"_query.bin"
			#We are going to normalise this sequence with the sklearn preprocessing algorithm to see what happens.		
			queryarray = sklearn.preprocessing.scale(np.array(squiggle),axis=0,with_mean=True,with_std=True,copy=True)
			
			dist, cost, path = mlpy.dtw_subsequence(queryarray,kmerhash2[id][ref]['Fprime'])
			result.append((dist,id,"F",path[1][0],ref))
			dist, cost, path = mlpy.dtw_subsequence(queryarray,kmerhash2[id][ref]['Rprime'])
			result.append((dist,id,"R",path[1][0],ref))
			

	return sorted(result,key=lambda result: result[0])[0][1],sorted(result,key=lambda result: result[0])[0][0],sorted(result,key=lambda result: result[0])[0][2],sorted(result,key=lambda result: result[0])[0][3],sorted(result,key=lambda result: result[0])[0][4]

######################################################################

def checkhairpin(squiggle,hairpin):
	queryarray = sklearn.preprocessing.scale(np.array(hairpin),axis=0,with_mean=True,with_std=True,copy=True)
	subjectarray = sklearn.preprocessing.scale(np.array(squiggle),axis=0,with_mean=True,with_std=True,copy=True)
	dist,cost,path=mlpy.dtw_subsequence(queryarray,subjectarray)
	return (dist,cost,path)

#######################################################
# Retrieve a model from the database rather than the  #
# expected data 									  #

def retrieve_model():
	model_kmers = dict()
	db = MySQLdb.connect(host=dbhost, user=dbusername, passwd=dbpass, port=dbport)
	sql = "SELECT * FROM minion_LomanLabz_013731_11rx_v2_3135.model_data where model like '%template%'"
	cursor.execute(sql)			
	kmerresults = cursor.fetchall()
	for line in kmerresults:
		kmer = line[2]
		mean = line[4]
		#print kmer,mean
		model_kmers[kmer]=mean
	return model_kmers

							
	kmer_means2=dict()
	kmer_means3=dict()
	
	for record in SeqIO.parse(ref_fasta, 'fasta'):
		counter = 1
		kmer_means3[record.id]=dict()
		for amplicon in args.ids:
			kmer_means3[record.id][counter]=dict()
			kmer_means3[record.id][counter]["F"]=list()
			kmer_means3[record.id][counter]["R"]=list()
			print amplicon
			seq = record.seq
			#print seq
			start = int(float(amplicon.split(':', 1 )[1].split('-',1)[0]))
			stop = int(float(amplicon.split(':', 1 )[1].split('-',1)[1]))
			newseq=seq[start:stop]
			print "Length of newseq:",len(newseq)
			shortregion=round(len(newseq)/2)
			print shortregion
			for x in range(int(shortregion),len(newseq)-kmer_len):
			print counter,"F",newseq
			newseq2 = revcomp = newseq.reverse_complement()
			#for x in range(len(newseq2)+1-kmer_len):
			for x in range(int(shortregion),len(newseq2)-kmer_len):
			print counter,"R",newseq2
			kmer_means3[record.id][counter]["Fprime"]=sklearn.preprocessing.scale(kmer_means3[record.id][counter]["F"], axis=0, with_mean=True, with_std=True, copy=True)
			kmer_means3[record.id][counter]["Rprime"]=sklearn.preprocessing.scale(kmer_means3[record.id][counter]["R"], axis=0, with_mean=True, with_std=True, copy=True)
			counter += 1
	
	for record in SeqIO.parse(ref_fasta, 'fasta'):
		counter = 1
		kmer_means2[record.id]=dict()
		for amplicon in args.ids:
			kmer_means2[record.id][counter]=dict()
			kmer_means2[record.id][counter]["F"]=list()
			kmer_means2[record.id][counter]["R"]=list()
			print amplicon
			seq = record.seq
			#print seq
			start = int(float(amplicon.split(':', 1 )[1].split('-',1)[0]))
			stop = int(float(amplicon.split(':', 1 )[1].split('-',1)[1]))
			newseq=seq[start:stop]
			print "Length of newseq:",len(newseq)
			shortregion=round(len(newseq)/2)
			print shortregion
			for x in range(int(shortregion)):
			print counter,"F",newseq
			newseq2 = revcomp = newseq.reverse_complement()
			#for x in range(len(newseq2)+1-kmer_len):
			for x in range(int(shortregion)):
			print counter,"R",newseq2
			kmer_means2[record.id][counter]["Fprime"]=sklearn.preprocessing.scale(kmer_means2[record.id][counter]["F"], axis=0, with_mean=True, with_std=True, copy=True)
			kmer_means2[record.id][counter]["Rprime"]=sklearn.preprocessing.scale(kmer_means2[record.id][counter]["R"], axis=0, with_mean=True, with_std=True, copy=True)
			counter += 1
			
			
		kmer_means[record.id]["Fprime"]=list()
		kmer_means[record.id]["Rprime"]=list()
		kmer_means[record.id]["Rprime"]=sklearn.preprocessing.scale(kmer_means[record.id]["R"], axis=0, with_mean=True, with_std=True, copy=True)


fasta_file = args.fasta
#fasta_file = "EM_079517.fasta"

if not os.path.exists('model.csv'):
	model_kmer_means = retrieve_model()
	# Saving the objects:
	#with open('model.pickle', 'wb') as f:
	#    pickle.dump([model_kmer_means], f)
	writer = csv.writer(open('model.csv', 'wb'))
	for key, value in model_kmer_means.items():
		writer.writerow([key, value])
else:
	# Getting back the objects:
	reader = csv.reader(open('model.csv','rb'))
	read = [((x),float(y)) for (x,y) in reader]
	model_kmer_means = dict(read)
	#with open('model.pickle', 'rb') as f:
	#	model_kmer_means = pickle.load(f)
	print "Loaded Model Correctly"
#print type(model_kmer_means)
#print model_kmer_means

seqlengths = get_seq_len(fasta_file)
#print seqlengths['EM_079517']
get_amplicons()

hairpinseq = "hairpin.fasta"

hairpinhash,hairpinhash2,hairpinhash3 = process_ref_fasta(hairpinseq,model_kmer_means)

#print hairpinhash

hairpin = hairpinhash['hairpin']['Fprime']
#print hairpin


for id in kmerhash2:
		for ref in kmerhash2[id]:
			#print id,ref
			#print (kmerhash2[id][ref]['F'])
			testarray = kmerhash2[id][ref]['Fprime']
			filename = id+"_"+str(ref)+"_F_subject.bin"
			filename = re.sub('\|','_',filename)			
			with open(filename, "wb") as f:
				f.write(ar.array("f", testarray))
			filename = id+"_"+str(ref)+"_F_subject.txt"
			filename = re.sub('\|','_',filename)
			np.savetxt(filename, testarray, delimiter=',')
			testarray = kmerhash2[id][ref]['Rprime']
			filename = id+"_"+str(ref)+"_R_subject.bin"
			filename = re.sub('\|','_',filename)			
			with open(filename, "wb") as f:
				f.write(ar.array("f", testarray))
			filename = id+"_"+str(ref)+"_R_subject.txt"
			filename = re.sub('\|','_',filename)
			np.savetxt(filename, testarray, delimiter=',')

goodcount=0
badcount=0
truegood=0
truebad=0
mismatch=0
notbasecalled=0	
reallygood=0
reallybad=0
reallytruegood=0
reallymismatch=0
reallynotbasecalled=0

###A dictionary to store our results in

readprediction=dict()
	
print "Now we are going to try and open the raw reads and do the same as we have done above..."
for filename in glob.glob(os.path.join(args.watchdir, '*.fast5')):
	print filename
	hdf = h5py.File(filename, 'r')
	for read in hdf['Analyses']['EventDetection_000']['Reads']:
		events = hdf['Analyses']['EventDetection_000']['Reads'][read]['Events'][()]
		event_collection=list()
		for event in events:
			event_collection.append(event[2])
#		print len(event_collection)
		matchseq,matchdist,orientation,matchstart,amplicon = squiggle_search2(event_collection[50:550],kmerhash2,'channel_id','read_id',seqlengths)
		print matchseq,matchdist,orientation,matchstart,amplicon		
		matchseq2,matchdist2,orientation2,matchstart2,amplicon2 = squiggle_search2(event_collection[-550:-50],kmerhash3,'channel_id','read_id',seqlengths)
		print matchseq2,matchdist2,orientation2,matchstart2,amplicon2
#		print len(kmerhash2[matchseq][amplicon][orientation])
		###An ideal read should be at least twice as many events as the length of the sequence so...
		testlen = len(kmerhash2[matchseq][amplicon][orientation])
#		print "Amplicon half length", testlen
#		print "Amplicon True length", (2 * testlen)
#		print "Query length", len(event_collection)
		###Test a custom array		
		if (amplicon2 == amplicon) and (orientation != orientation2):
			print "**************************MATCH***********************"
			if (amplicon not in readprediction):
			if ("0" not in readprediction[amplicon]):
				readprediction[amplicon]["0"]=dict()
			if (filename not in readprediction[amplicon]["0"]):
				readprediction[amplicon]["0"][filename]=dict()
			readprediction[amplicon]["0"][filename]["name"]=filename
			readprediction[amplicon]["0"][filename]["matchdistance"]=matchdist
			break
		hairpin2 = np.array([72.29113922,73.52673668,73.55946867,67.93627677,70.78341073,69.77046482,59.71001375,66.88501831,107.0488161,113.3369954,117.0188287,123.0165146,124.2285211,124.9807286,125.1277794,123.500317,124.180417,121.7045039,119.8741114,108.5086081,74.11850918,71.26983532,72.25233459,53.74366443,53.8004576,54.82023287,54.88335392,55.01291259,55.17210148,52.3131189,55.09850506,53.61134266,53.8489204,52.46566847,53.51238622,67.85040665,76.43446973,94.73609366,92.39555688,95.32671653,105.49724,112.8739329,112.0287627,118.7370562,117.5573792,119.1706588,122.5249539,121.5482435,122.3349303,123.0971417,120.6362789,120.5928762,121.932797,120.9672013,119.6665197,118.0263136,117.8794998,117.0123523,117.7958905,104.752265,81.78663782,76.73602154,84.10819387,76.19999462,70.38227417])
		hairpin3 = sklearn.preprocessing.scale(hairpin2, axis=0, with_mean=True, with_std=True, copy=True)
		dist,cost,path=checkhairpin(event_collection,hairpin3)
		print dist,path[1][0],path[1][-1],(path[1][-1]-path[1][1])
		filetocheck = os.path.split(filename)
#		try:
#			print "Template Event Guess:", (path[1][0]-1),"True:",templatelength[filetocheck[1]]
#			print "Complement Event Guess:", (len(event_collection)-path[1][-1]+1),"True:",complementlength[filetocheck[1]]
#		except Exception, err:
#			print "Problem determining read lengths"
		
		
		if ( 0.5 <  ((path[1][0]-1)/(len(event_collection)-path[1][-1]+1)) < 2) or (len(event_collection)-path[1][-1]+1)>(path[1][0]-1) or ((8 * testlen) > len(event_collection) > (2 * testlen)):
			###If this criterion is satisfied we are going to put this read in tier 1 - i.e the best quality
			reallygood += 1			
			if (8 * testlen) > len(event_collection) > (2 * testlen):
				if (amplicon not in readprediction):
				if ("1" not in readprediction[amplicon]):
					readprediction[amplicon]["1"]=dict()					
				if (filename not in readprediction[amplicon]["1"]):
					readprediction[amplicon]["1"][filename]=dict()					
				readprediction[amplicon]["1"][filename]["name"]=filename
				readprediction[amplicon]["1"][filename]["matchdistance"]=matchdist
			else:
				if (amplicon not in readprediction):
				if ("2" not in readprediction[amplicon]):
					readprediction[amplicon]["2"]=dict()					
				if (filename not in readprediction[amplicon]["2"]):
					readprediction[amplicon]["2"][filename]=dict()					
				readprediction[amplicon]["2"][filename]["name"]=filename
				readprediction[amplicon]["2"][filename]["matchdistance"]=matchdist
		else:
			if (8 * testlen) > len(event_collection) > (2 * testlen):
				if (amplicon not in readprediction):
					
				if ("3" not in readprediction[amplicon]):
					readprediction[amplicon]["3"]=dict()
					
				if (filename not in readprediction[amplicon]["3"]):
					readprediction[amplicon]["3"][filename]=dict()
					
				readprediction[amplicon]["3"][filename]["name"]=filename
				readprediction[amplicon]["3"][filename]["matchdistance"]=matchdist
				#### A read here is in tier 2 - we don't know if it is any good or not...
				print "Looking promising"
			else:
				if (amplicon not in readprediction):
					
				if ("4" not in readprediction[amplicon]):
					readprediction[amplicon]["4"]=dict()
					
				if (filename not in readprediction[amplicon]["4"]):
					readprediction[amplicon]["4"][filename]=dict()
					
				readprediction[amplicon]["4"][filename]["name"]=filename
				readprediction[amplicon]["4"][filename]["matchdistance"]=matchdist
				#### A read here we think is really bad.
				print "**********BAD READ************"
				badcount+=1
				###Even though the read looks bad lets check it anyway!	
				dist,cost,path=checkhairpin(event_collection,hairpin3)
				print dist,path[1][0],path[1][-1],(path[1][-1]-path[1][1])
				filetocheck = os.path.split(filename)
				
				
				
#			print "THIS DOESN'T SATISFY ME"
			reallybad += 1
		
		
		### We want to know if this read really is good or bad as defined by nanopore - we already have this in the database so we can check from there!
		filetocheck = os.path.split(filename)
		
	hdf.close()
print "Amplicon Read Counts"	
for amplicon in readprediction:
	numberofreads = 0
	try:
		if len(readprediction[amplicon]["0"].keys()) > 0:
			numberofreads += len(readprediction[amplicon]["0"].keys())
	except Exception, err:
			print "",
	try:
		if len(readprediction[amplicon]["1"].keys()) > 0:
			numberofreads += len(readprediction[amplicon]["1"].keys())
	except Exception, err:
			print "",
	try:
		if len(readprediction[amplicon]["2"].keys()) > 0:
			numberofreads += len(readprediction[amplicon]["2"].keys())
	except Exception, err:
			print "",
	try:
		if len(readprediction[amplicon]["3"].keys()) > 0:
			numberofreads += len(readprediction[amplicon]["3"].keys())
	except Exception, err:
			print "",
	try:
		if len(readprediction[amplicon]["4"].keys()) > 0:
			numberofreads += len(readprediction[amplicon]["4"].keys())
	except Exception, err:
			print "",
	print ""
	print "Amplicon Number:",amplicon,"Reads:",numberofreads

print ""
print "Copying Amplicon Data"
for amplicon in readprediction:
	print "Amplicon Number",amplicon
	counter = 0
	try:
		if (len(readprediction[amplicon]["0"].keys())>0):
			print len(readprediction[amplicon]["0"].keys())
			if (counter < args.depth):			
				ordered0 = OrderedDict(sorted(readprediction[amplicon]["0"].iteritems(), key=lambda x: x[1]['matchdistance']))
				#print ordered1
				#ordered0 = readprediction[amplicon]["0"].iteritems()
				for read in ordered0:
					print read, ordered0[read]["matchdistance"]
					if not os.path.exists('prefilteredraw'):
						os.makedirs('prefilteredraw')
					destdir = os.path.join('prefilteredraw')    	
					if not os.path.exists(destdir):
						os.makedirs(destdir)
					try:
						filetocheck = os.path.split(read)
						#subdir=readlookup[filetocheck[1]]
						#sourcefile = os.path.join(filetocheck[0],'downloads',subdir,filetocheck[1])
						sourcefile = read
						destfile = os.path.join(destdir,filetocheck[1])
			
						#sourcefile = filename
						#filename = ids[5]+'.fast5'
						#destfile = os.path.join(destdir,filename)
#						print "sourcefile is:",sourcefile
#						print "destfile is:",destfile
						try:
							#os.symlink(sourcefile, destfile)
							shutil.copy(sourcefile,destfile)
						except Exception, err:
							print "File Copy Failed",err
					except Exception, err:
						print "Weird bug I don't GROK"
					counter += 1
					if counter >= args.depth:
						break
				#for read in sorted(readprediction[amplicon]["0"].iteritems(), key=lambda x: x[1]["matchdistance"]):
			#	print read, readprediction[amplicon]["0"][read]["matchdistance"]
	except Exception, err:
		print "No reads of class 0"
	try:
		if (len(readprediction[amplicon]["1"].keys())>0):
			print len(readprediction[amplicon]["1"].keys())
			if (counter < args.depth):			
				ordered1 = OrderedDict(sorted(readprediction[amplicon]["1"].iteritems(), key=lambda x: x[1]['matchdistance']))
				#print ordered1
				#ordered1 = readprediction[amplicon]["1"].iteritems()
				for read in ordered1:
					print read, ordered1[read]["matchdistance"]
					if not os.path.exists('prefilteredraw'):
						os.makedirs('prefilteredraw')
					destdir = os.path.join('prefilteredraw')    	
					if not os.path.exists(destdir):
						os.makedirs(destdir)
					try:
						filetocheck = os.path.split(read)
						#subdir=readlookup[filetocheck[1]]
						#sourcefile = os.path.join(filetocheck[0],'downloads',subdir,filetocheck[1])
						sourcefile = read
						destfile = os.path.join(destdir,filetocheck[1])
			
						#sourcefile = filename
						#filename = ids[5]+'.fast5'
						#destfile = os.path.join(destdir,filename)
#						print "sourcefile is:",sourcefile
#						print "destfile is:",destfile
						try:
							#os.symlink(sourcefile, destfile)
							shutil.copy(sourcefile,destfile)
						except Exception, err:
							print "File Copy Failed",err
					except Exception, err:
						print "Weird bug I don't GROK"
					counter += 1
					if counter >= args.depth:
						break
				#for read in sorted(readprediction[amplicon]["1"].iteritems(), key=lambda x: x[1]["matchdistance"]):
			#	print read, readprediction[amplicon]["1"][read]["matchdistance"]
	except Exception, err:
		print "No reads of class 1"
	try:
		if (len(readprediction[amplicon]["2"].keys())>0):
			print len(readprediction[amplicon]["2"].keys())
			if (counter < args.depth):
				ordered2 = OrderedDict(sorted(readprediction[amplicon]["2"].iteritems(), key=lambda x: x[2]['matchdistance']))
				#ordered2 = readprediction[amplicon]["2"].iteritems()
				#print ordered1
				for read in ordered2:
					print read, ordered2[read]["matchdistance"]
					if not os.path.exists('prefilteredraw'):
						os.makedirs('prefilteredraw')
					destdir = os.path.join('prefilteredraw')    	
					if not os.path.exists(destdir):
						os.makedirs(destdir)
					try:
						filetocheck = os.path.split(read)
						#subdir=readlookup[filetocheck[1]]
						#sourcefile = os.path.join(filetocheck[0],'downloads',subdir,filetocheck[1])
						sourcefile = read
						destfile = os.path.join(destdir,filetocheck[1])
						#sourcefile = filename
						#filename = ids[5]+'.fast5'
						#destfile = os.path.join(destdir,filename)
#						print "sourcefile is:",sourcefile
#						print "destfile is:",destfile
						try:
							#os.symlink(sourcefile, destfile)
							shutil.copy(sourcefile,destfile)
						except Exception, err:
							print "File Copy Failed",err
					except Exception, err:
						print "Weird bug I don't GROK"
					
					counter += 1
					if counter >= args.depth:
						break
	except Exception, err:
		print "No reads of class 2"
	try:
		if (len(readprediction[amplicon]["3"].keys())>0):
			print len(readprediction[amplicon]["3"].keys())
			if (counter < args.depth):	
				ordered3 = OrderedDict(sorted(readprediction[amplicon]["3"].iteritems(), key=lambda x: x[3]['matchdistance']))
				#ordered3 = readprediction[amplicon]["3"].iteritems()
				#print ordered1
				for read in ordered3:
					print read, ordered3[read]["matchdistance"]
					if not os.path.exists('prefilteredraw'):
						os.makedirs('prefilteredraw')
					destdir = os.path.join('prefilteredraw')    	
					if not os.path.exists(destdir):
						os.makedirs(destdir)
					try:
						filetocheck = os.path.split(read)
						#subdir=readlookup[filetocheck[1]]
						#sourcefile = os.path.join(filetocheck[0],'downloads',subdir,filetocheck[1])
						sourcefile = read
						destfile = os.path.join(destdir,filetocheck[1])
						#sourcefile = filename
						#filename = ids[5]+'.fast5'
						#destfile = os.path.join(destdir,filename)
#						print "sourcefile is:",sourcefile
#						print "destfile is:",destfile
						try:
							#os.symlink(sourcefile, destfile)
							shutil.copy(sourcefile,destfile)							
						except Exception, err:
							print "File Copy Failed",err
					except Exception, err:
						print "Weird bug I don't GROK"

					counter += 1
					if counter >= args.depth:
						break
	except Exception, err:
		print "No reads of class 3"
	try:
		if (len(readprediction[amplicon]["4"].keys())>0):
			print len(readprediction[amplicon]["4"].keys())
			if (counter < args.depth):
				ordered4 = OrderedDict(sorted(readprediction[amplicon]["4"].iteritems(), key=lambda x: x[4]['matchdistance']))
				#ordered4 = readprediction[amplicon]["4"].iteritems()
				#print ordered1
				for read in ordered4:
					print read, ordered4[read]["matchdistance"]
					if not os.path.exists('prefilteredraw'):
						os.makedirs('prefilteredraw')
					destdir = os.path.join('prefilteredraw')    	
					if not os.path.exists(destdir):
						os.makedirs(destdir)
					try:
						filetocheck = os.path.split(read)
						#subdir=readlookup[filetocheck[1]]
						#sourcefile = os.path.join(filetocheck[0],'downloads',subdir,filetocheck[1])
						sourcefile = read
						destfile = os.path.join(destdir,filetocheck[1])
						#sourcefile = filename
						#filename = ids[5]+'.fast5'
						#destfile = os.path.join(destdir,filename)
#						print "sourcefile is:",sourcefile
#						print "destfile is:",destfile
						try:
							#os.symlink(sourcefile, destfile)
							shutil.copy(sourcefile,destfile)
						except Exception, err:
							print "File Copy Failed",err
					except Exception, err:
						print "Weird bug I don't GROK"

					counter += 1
					if counter >= args.depth:
						break
	except Exception, err:
		print "No reads of class 4"