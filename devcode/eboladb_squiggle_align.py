#!/usr/bin/env pythonimport sys, os, reimport time#import datetime#import logging#from watchdog.observers.polling import PollingObserver as Observer#from watchdog.events import FileSystemEventHandlerimport threading, thread#import h5pyfrom Bio import SeqIOfrom StringIO import StringIOimport MySQLdb#import subprocessimport string#import configargparse#from warnings import filterwarnings#import socket#import hashlib
import mlpy
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



parser = configargparse.ArgParser(description='eboladb_squiggle_align: A program providing the ability to determine which region of the ebola genome and individual read is derived from.')
parser.add('-fasta', '--reference_fasta_file', type=str, dest='fasta', required=True, default=None, help="The fasta format file for the reference sequence for your organism.")
parser.add('-ids', nargs = '*', dest='ids',required=True, help = 'A list of start and stop positions for each amplicon from the reference genome - should be space separated with fasta_name:start-stop.\n e.g.\n EM_079517:27-1938 EM_078517:1927-3828 EM_078517:3823-5718 EM_078517:5759-7633 EM_078517:7601-10007 EM_078517:9550-10921 EM_078517:10944-12354 EM_078517:12354-14252 EM_078517:14253-15680 EM_078517:15691-17087 EM_078517:16632-18553\n ')
parser.add('-w', '--watch-dir', type=str, required=True, default=None, help="The path to the folder containing the downloads directory with fast5 reads to analyse - e.g. C:\data\minion\downloads (for windows).", dest='watchdir')
parser.add('-d', '--depth',type=int, required=True, default=None, help = 'The desired coverage depth for each amplicon. Note this is unlikely to be achieved for each amplicon and should be an overestimate of the minimum coverage required.', dest='depth')
args = parser.parse_args()



######################################################
# Connect to the database 							 #
######################################################
dbhost = 'localhost'
dbusername = 'minion'
dbpass = 'nan0p0re'
dbport = 3306

limitbases = "50,200"

######################################################def process_model_file(model_file):	model_kmers = dict()	with open(model_file, 'rb') as csv_file:		reader = csv.reader(csv_file, delimiter="\t")    		d = list(reader)		#print d		for r in range(1, len(d)):			#print r, d[r]			kmer = d[r][0]			mean = d[r][2]			#print r, kmer, mean			model_kmers[kmer]=mean	return 	model_kmers	
	
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
def squiggle_search(squiggle,kmerhash2,channel_id,read_id,seqlen):
	result=[]
	for id in kmerhash2:
		for ref in kmerhash2[id]:
			#print len(kmerhash2[id][ref]['F'])
			queryfile=str(channel_id)+"_"+str(read_id)+"_query.bin"
			#We are going to normalise this sequence with the sklearn preprocessing algorithm to see what happens.		
			queryarray = sklearn.preprocessing.scale(np.array(squiggle),axis=0,with_mean=True,with_std=True,copy=True)
			#queryarray = np.array(squiggle)
			with open(queryfile, "wb") as f:
				f.write(ar.array("f", queryarray))
			subjectfile = id+"_"+str(ref)+"_"+"F"+"_subject.bin"
			subjectfile = re.sub('\|','_',subjectfile)
			#seqlen2 = str(seqlen[id])
			commands = queryfile+' '+subjectfile+' 200 '+str(len(kmerhash2[id][ref]['F']))+' 0.05'
			#print commands
			#current = str(multiprocessing.current_process())
			#currentnum=int(re.search(r'\d+', current).group())
			gpucode=str()
			#if (currentnum % 2 == 0):
				#print "Even"
			gpucode='./GPU-DTW '
			#else:
			#	#print "Odd"
			#	gpucode='./GPU-DTW '
			runcommand = gpucode+commands
			location = ()
			distance = ()
			for line in runProcess(runcommand.split()):
				if "Location" in line:			
					location = int(line.split(': ',1)[1].rstrip('\n'))
				if "Distance" in line:
					distance = float(line.split(': ',1)[1].rstrip('\n'))
			result.append((distance,id,"F",location,ref))
			subjectfile = id+"_"+str(ref)+"_"+"R"+"_subject.bin"
			subjectfile = re.sub('\|','_',subjectfile)
			#seqlen2 = str(seqlen[id])
			commands = queryfile+' '+subjectfile+' 200 '+str(len(kmerhash2[id][ref]['R']))+' 0.05'
			#print commands
			#current = str(multiprocessing.current_process())
			#currentnum=int(re.search(r'\d+', current).group())
			gpucode=str()
			#if (currentnum % 2 == 0):
				#print "Even"
			gpucode='./GPU-DTW '
			#else:
			#	#print "Odd"
			#	gpucode='./GPU-DTW '
			runcommand = gpucode+commands
			location = ()
			distance = ()
			for line in runProcess(runcommand.split()):
				if "Location" in line:			
					location = int(line.split(': ',1)[1].rstrip('\n'))
				if "Distance" in line:
					distance = float(line.split(': ',1)[1].rstrip('\n'))
			result.append((distance,id,"R",location,ref))
			os.remove(queryfile)

	return sorted(result,key=lambda result: result[0])[0][1],sorted(result,key=lambda result: result[0])[0][0],sorted(result,key=lambda result: result[0])[0][2],sorted(result,key=lambda result: result[0])[0][3],sorted(result,key=lambda result: result[0])[0][4]

#######################################################################
def squiggle_search3(squiggle,kmerhash2,channel_id,read_id,seqlen):
	result=[]
	for id in kmerhash2:
		for ref in kmerhash2[id]:
			#print len(kmerhash2[id][ref]['F'])
			queryfile=str(channel_id)+"_"+str(read_id)+"_query.txt"
			#We are going to normalise this sequence with the sklearn preprocessing algorithm to see what happens.		
			queryarray = sklearn.preprocessing.scale(np.array(squiggle),axis=0,with_mean=True,with_std=True,copy=True)
			#queryarray = np.array(squiggle)
			np.savetxt(queryfile, queryarray, delimiter=',')
			subjectfile = id+"_"+str(ref)+"_"+"F"+"_subject.txt"
			subjectfile = re.sub('\|','_',subjectfile)
			#seqlen2 = str(seqlen[id])
			commands = subjectfile+' '+queryfile+' 200	 0.1'
			#print commands
			#current = str(multiprocessing.current_process())
			#currentnum=int(re.search(r'\d+', current).group())
			gpucode=str()
			#if (currentnum % 2 == 0):
				#print "Even"
			gpucode='./UCR_DTW '
			#else:
			#	#print "Odd"
			#	gpucode='./GPU-DTW '
			runcommand = gpucode+commands
			#print runcommand
			location = ()
			distance = ()
			for line in runProcess(runcommand.split()):
				#print line
				if "Location" in line:			
					location = int(line.split(': ',1)[1].rstrip('\n'))
					#print location
				if "Distance" in line:
					distance = float(line.split(': ',1)[1].rstrip('\n'))
					#print distance
			result.append((distance,id,"F",location,ref))
			subjectfile = id+"_"+str(ref)+"_"+"R"+"_subject.txt"
			subjectfile = re.sub('\|','_',subjectfile)
			#seqlen2 = str(seqlen[id])
			commands = subjectfile+' '+queryfile+' 200'+' 0.1'
			#print commands
			#current = str(multiprocessing.current_process())
			#currentnum=int(re.search(r'\d+', current).group())
			gpucode=str()
			#if (currentnum % 2 == 0):
				#print "Even"
			gpucode='./UCR_DTW '
			#else:
			#	#print "Odd"
			#	gpucode='./GPU-DTW '
			runcommand = gpucode+commands
			#print runcommand
			location = ()
			distance = ()
			for line in runProcess(runcommand.split()):
				if "Location" in line:			
					location = int(line.split(': ',1)[1].rstrip('\n'))
				if "Distance" in line:
					distance = float(line.split(': ',1)[1].rstrip('\n'))
			result.append((distance,id,"R",location,ref))
			os.remove(queryfile)

	return sorted(result,key=lambda result: result[0])[0][1],sorted(result,key=lambda result: result[0])[0][0],sorted(result,key=lambda result: result[0])[0][2],sorted(result,key=lambda result: result[0])[0][3],sorted(result,key=lambda result: result[0])[0][4]

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
# expected data 									  ########################################################

def retrieve_model():
	model_kmers = dict()
	db = MySQLdb.connect(host=dbhost, user=dbusername, passwd=dbpass, port=dbport)	cursor = db.cursor() 
	sql = "SELECT * FROM minion_LomanLabz_013731_11rx_v2_3135.model_data where model like '%template%'"
	cursor.execute(sql)			
	kmerresults = cursor.fetchall()
	for line in kmerresults:
		kmer = line[2]
		mean = line[4]
		#print kmer,mean
		model_kmers[kmer]=mean
	return model_kmers

							######################################################def process_ref_fasta(ref_fasta,model_kmer_means):	print "processing the reference fasta."	#ref_kmers=dict()	kmer_len=5	kmer_means=dict()
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
			for x in range(int(shortregion),len(newseq)-kmer_len):				kmer = str(newseq[x:x+kmer_len])				kmer_means3[record.id][counter]["F"].append(float(model_kmer_means[kmer]))
			print counter,"F",newseq
			newseq2 = revcomp = newseq.reverse_complement()
			#for x in range(len(newseq2)+1-kmer_len):
			for x in range(int(shortregion),len(newseq2)-kmer_len):				kmer = str(newseq2[x:x+kmer_len])				kmer_means3[record.id][counter]["R"].append(float(model_kmer_means[kmer]))
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
			for x in range(int(shortregion)):				kmer = str(newseq[x:x+kmer_len])				kmer_means2[record.id][counter]["F"].append(float(model_kmer_means[kmer]))
			print counter,"F",newseq
			newseq2 = revcomp = newseq.reverse_complement()
			#for x in range(len(newseq2)+1-kmer_len):
			for x in range(int(shortregion)):				kmer = str(newseq2[x:x+kmer_len])				kmer_means2[record.id][counter]["R"].append(float(model_kmer_means[kmer]))
			print counter,"R",newseq2
			kmer_means2[record.id][counter]["Fprime"]=sklearn.preprocessing.scale(kmer_means2[record.id][counter]["F"], axis=0, with_mean=True, with_std=True, copy=True)
			kmer_means2[record.id][counter]["Rprime"]=sklearn.preprocessing.scale(kmer_means2[record.id][counter]["R"], axis=0, with_mean=True, with_std=True, copy=True)
			counter += 1
			
					for record in SeqIO.parse(ref_fasta, 'fasta'):		kmer_means[record.id]=dict()		kmer_means[record.id]["F"]=list()		kmer_means[record.id]["R"]=list()
		kmer_means[record.id]["Fprime"]=list()
		kmer_means[record.id]["Rprime"]=list()		print "ID", record.id		print "length", len(record.seq)			print "FORWARD STRAND"		seq = record.seq		for x in range(len(seq)+1-kmer_len):			kmer = str(seq[x:x+kmer_len])			kmer_means[record.id]["F"].append(float(model_kmer_means[kmer]))			#if model_kmer_means[kmer]:				#print x, kmer, model_kmer_means[kmer]		print "REVERSE STRAND"		seq = revcomp = record.seq.reverse_complement()		for x in range(len(seq)+1-kmer_len):			kmer = str(seq[x:x+kmer_len])			kmer_means[record.id]["R"].append(float(model_kmer_means[kmer]))		kmer_means[record.id]["Fprime"]=sklearn.preprocessing.scale(kmer_means[record.id]["F"], axis=0, with_mean=True, with_std=True, copy=True)
		kmer_means[record.id]["Rprime"]=sklearn.preprocessing.scale(kmer_means[record.id]["R"], axis=0, with_mean=True, with_std=True, copy=True)	return kmer_means,kmer_means2,kmer_means3######################################################


fasta_file = args.fasta
#fasta_file = "EM_079517.fasta"
model_kmer_means = retrieve_model()kmerhash,kmerhash2,kmerhash3 = process_ref_fasta(fasta_file,model_kmer_means)
seqlengths = get_seq_len(fasta_file)
print seqlengths['EM_079517']
get_amplicons()

hairpinseq = "hairpin.fasta"

hairpinhash,hairpinhash2,hairpinhash3 = process_ref_fasta(hairpinseq,model_kmer_means)

#print hairpinhash

hairpin = hairpinhash['hairpin']['Fprime']
print hairpin


for id in kmerhash2:
		for ref in kmerhash2[id]:
			print id,ref
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


#ecit




db = MySQLdb.connect(host=dbhost, user=dbusername, passwd=dbpass, port=dbport)cursor = db.cursor() 
#sql = "use minion_PLSP57501_2014_10_10_DSmin1_run2_LambdaSK002_5041"
#sql = "use minion_PLSP57501_20140909_JA_defA_4434"
sql = "use minion_LomanLabz_013731_11rx_v2_3135"
print sql
cursor.execute(sql)

numbers = range(0,10)

### To store the ONT good/fail characteristic
readlookup=dict()

### Lookup for readfile
readfilelookup=dict()

### Lookup for template length
templatelength=dict()

### Lookup for complement length
complementlength=dict()

for number in numbers:
	#print number

	#sql = "SELECT basename_id ,pos,flag, channel, read_id,tracking_id.basename,file_path FROM caller_basecalled_template_%s inner join align_sam_basecalled_2d using (basename_id) inner join config_general using (basename_id) inner join tracking_id using (basename_id) group by basename_id " %(number)
	sql = "SELECT basename_id ,channel, read_id,tracking_id.basename,file_path FROM caller_basecalled_template_%s inner join config_general using (basename_id) inner join tracking_id using (basename_id) group by basename_id " %(number)
	#sql = "SELECT basename_id , channel, read_id FROM caller_basecalled_template_%s inner join config_general using (basename_id) where basename_id not in (select basename_id from align_sam_basecalled_template) group by basename_id" %(number)
	print sql	cursor.execute(sql)
					
	
	
	basenameids = cursor.fetchall()
	
	#print basenameids
	
	for ids in basenameids:
		#print ids[4]
		elements = ids[4].split('/')
		readlookup[elements[-1]]=elements[-2]
		readfilelookup[elements[-1]]=ids[4]
		#print type(elements)
		#print elements
		#print elements[-1]
		#print elements[-2]
#		sql = "SELECT count(*) FROM caller_basecalled_template_%s where basename_id = %s" % (number,ids[0])
		#print sql
		#sql = "SELECT mean FROM caller_basecalled_template where basename_id = %s " % (ids[0]) 
#		cursor.execute(sql)
#		length = cursor.fetchall()
		#print length[0][0]
#		templatelength[elements[-1]]=length[0][0]
		
#		sql = "SELECT count(*) FROM caller_basecalled_complement_%s where basename_id = %s" % (number,ids[0])
		#print sql
		#sql = "SELECT mean FROM caller_basecalled_template where basename_id = %s " % (ids[0]) 
#		cursor.execute(sql)
#		length = cursor.fetchall()
		#print length[0][0]
#		complementlength[elements[-1]]=length[0][0]		
##		meanlist = list()
##		for mean in means:
##			meanlist.append(float(mean[0]))
		
		#print len(meanlist)
		#print squiggle_search3(meanlist,kmerhash2,ids[3],ids[4],seqlengths),
		#print squiggle_search(meanlist,kmerhash2,ids[3],ids[4],seqlengths),
#		matchseq,dist,orientation,matchstart,amplicon = squiggle_search2(meanlist,kmerhash2,ids[3],ids[4],seqlengths)
#		print matchseq,dist,orientation,matchstart,amplicon,
#		x=np.array(meanlist)
#		Xprime = sklearn.preprocessing.scale(x, axis=0, with_mean=True, with_std=True, copy=True)
#		result_sum=dict()
#		scale_result_sum=dict()
		
#		if not os.path.exists('symsplit'):
#			os.makedirs('symsplit')
#		destdir = os.path.join('symsplit',str(amplicon),'downloads')    	
#		if not os.path.exists(destdir):
#			os.makedirs(destdir)
		
#		print ids[6],
#		sourcefile = ids[6]

#		filename = ids[5]+'.fast5'
#		destfile = os.path.join(destdir,filename)
#		print "sourcefile is:",sourcefile
#		print "destfile is:",destfile
		#os.symlink(sourcefile, destfile)
		
#		for seqid in kmerhash:
#			y=np.array(kmerhash[seqid]["F"])
#			z=np.array(kmerhash[seqid]["R"])
			
#			Yprime = np.array(kmerhash[seqid]["Fprime"])
#			Zprime = np.array(kmerhash[seqid]["Rprime"])
			
#			dist, cost, path = mlpy.dtw_subsequence(Xprime,Yprime)
#			secondarray=path[1]
#			scale_result_sum[dist]=dict()
#			scale_result_sum[dist][seqid]=dict()
#			scale_result_sum[dist][seqid]["F"]=list()
#			scale_result_sum[dist][seqid]["F"].append(secondarray[0])
			
#			dist, cost, path = mlpy.dtw_subsequence(Xprime,Zprime)
#			secondarray=path[1]
#			scale_result_sum[dist]=dict()
#			scale_result_sum[dist][seqid]=dict()
#			scale_result_sum[dist][seqid]["R"]=list()
#			scale_result_sum[dist][seqid]["R"].append(secondarray[0])
			
#			print (ids[0],ids[1],ids[2],ids[3],ids[4]),
#			#print (ids[0],ids[1],ids[2]),

#			for key in sorted(scale_result_sum):
#			    for seqid in scale_result_sum[key]:
#			    	for direct in scale_result_sum[key][seqid]:
#			    		for pos in scale_result_sum[key][seqid][direct]:
#			    			if not 25000<pos<30000:
#			    				print "SEQID",(seqid),
#			    				print "POS",(pos),
#			    				print "DIRECT",(direct),
#			    				print "KEY",key
#			    			else:
#			    				print "SEQID",(seqid),
#			    				print "POS",(pos),
#			    				print "DIRECT",(direct),
#			    				print "KEY",key
#			    break
		

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
	#print type(hdf)
	#for hdfelement in hdf:
	#	print hdfelement
	#print hdf['Sequences']
	#for hdfsequenceelement in hdf['Sequences']:
	#	print "\t",hdfsequenceelement
	for read in hdf['Analyses']['EventDetection_000']['Reads']:
	#	print "Read:",read
	#	print type(read)
		#for bit in read:
		#	print bit
	#	print hdf['Analyses']['EventDetection_000']['Reads'][read]
	#	for bit in hdf['Analyses']['EventDetection_000']['Reads'][read]:
	#		print bit
	#	print hdf['Analyses']['EventDetection_000']['Reads'][read]['Events']
		events = hdf['Analyses']['EventDetection_000']['Reads'][read]['Events'][()]
		event_collection=list()
		for event in events:
			#print event[2]
			event_collection.append(event[2])
		print len(event_collection)
		#print event_collection[50:250]
		matchseq,matchdist,orientation,matchstart,amplicon = squiggle_search2(event_collection[50:550],kmerhash2,'channel_id','read_id',seqlengths)
		print matchseq,matchdist,orientation,matchstart,amplicon
		
		matchseq2,matchdist2,orientation2,matchstart2,amplicon2 = squiggle_search2(event_collection[-550:-50],kmerhash3,'channel_id','read_id',seqlengths)
		print matchseq2,matchdist2,orientation2,matchstart2,amplicon2
		
		###Given this information we can determine the approximate number of expected events from the read.
		print len(kmerhash2[matchseq][amplicon][orientation])
		###An ideal read should be at least twice as many events as the length of the sequence so...
		testlen = len(kmerhash2[matchseq][amplicon][orientation])
		print "Amplicon half length", testlen
		print "Amplicon True length", (2 * testlen)
		print "Query length", len(event_collection)
		###Test a custom array
		
		if (amplicon2 == amplicon) and (orientation != orientation2):
			print "**************************BINGO***********************"
			if (amplicon not in readprediction):					readprediction[amplicon]=dict()
					
			if ("0" not in readprediction[amplicon]):
				readprediction[amplicon]["0"]=dict()
				
			if (filename not in readprediction[amplicon]["0"]):
				readprediction[amplicon]["0"][filename]=dict()
					
			readprediction[amplicon]["0"][filename]["name"]=filename
			readprediction[amplicon]["0"][filename]["matchdistance"]=matchdist

			break
		hairpin2=np.array([72.29113922,73.52673668,73.55946867,67.93627677,70.78341073,69.77046482,59.71001375,66.88501831,107.0488161,113.3369954,117.0188287,123.0165146,124.2285211,124.9807286,125.1277794,123.500317,124.180417,121.7045039,119.8741114,108.5086081,74.11850918,71.26983532,72.25233459,53.74366443,53.8004576,54.82023287,54.88335392,55.01291259,55.17210148,52.3131189,55.09850506,53.61134266,53.8489204,52.46566847,53.51238622,67.85040665,76.43446973,94.73609366,92.39555688,95.32671653,105.49724,112.8739329,112.0287627,118.7370562,117.5573792,119.1706588,122.5249539,121.5482435,122.3349303,123.0971417,120.6362789,120.5928762,121.932797,120.9672013,119.6665197,118.0263136,117.8794998,117.0123523,117.7958905,104.752265,81.78663782,76.73602154,84.10819387,76.19999462,70.38227417])
		hairpin3 = sklearn.preprocessing.scale(hairpin2, axis=0, with_mean=True, with_std=True, copy=True)
		dist,cost,path=checkhairpin(event_collection,hairpin3)
		print dist,path[1][0],path[1][-1],(path[1][-1]-path[1][1])
		filetocheck = os.path.split(filename)
		try:
			print "Template Event Guess:", (path[1][0]-1),"True:",templatelength[filetocheck[1]]
			print "Complement Event Guess:", (len(event_collection)-path[1][-1]+1),"True:",complementlength[filetocheck[1]]
		except Exception, err:
			print "Problem determining read lengths"
		
		
		#if (len(event_collection)-path[1][-1]+1)>(path[1][0]-1):
		if ( 0.5 <  ((path[1][0]-1)/(len(event_collection)-path[1][-1]+1)) < 2) or (len(event_collection)-path[1][-1]+1)>(path[1][0]-1) or ((8 * testlen) > len(event_collection) > (2 * testlen)):
			###If this criterion is satisfied we are going to put this read in tier 1 - i.e the best quality
			if not os.path.exists('prefiltered'):
				os.makedirs('prefiltered')
			destdir = os.path.join('prefiltered','downloads')    	
			if not os.path.exists(destdir):
				os.makedirs(destdir)
			filetocheck = os.path.split(filename)
			try:
				sourcefile = readfilelookup[filetocheck[1]]
				destfile = os.path.join(destdir,filetocheck[1])
			
				#sourcefile = filename
				#filename = ids[5]+'.fast5'
				#destfile = os.path.join(destdir,filename)
				#print "sourcefile is:",sourcefile
				#print "destfile is:",destfile
				#try:
				#	os.symlink(sourcefile, destfile)
				#except Exception, err:
				#	print "File Copy Failed",err
			except Exception, err:
				print "Weird bug I don't GROK"
			print "THIS IS REALLY GOOD"
			reallygood += 1
			
			if (8 * testlen) > len(event_collection) > (2 * testlen):
			
				if (amplicon not in readprediction):					readprediction[amplicon]=dict()
					
				if ("1" not in readprediction[amplicon]):
					readprediction[amplicon]["1"]=dict()
					
				if (filename not in readprediction[amplicon]["1"]):
					readprediction[amplicon]["1"][filename]=dict()
					
				readprediction[amplicon]["1"][filename]["name"]=filename
				readprediction[amplicon]["1"][filename]["matchdistance"]=matchdist
			else:
				if (amplicon not in readprediction):					readprediction[amplicon]=dict()
					
				if ("2" not in readprediction[amplicon]):
					readprediction[amplicon]["2"]=dict()
					
				if (filename not in readprediction[amplicon]["2"]):
					readprediction[amplicon]["2"][filename]=dict()
					
				readprediction[amplicon]["2"][filename]["name"]=filename
				readprediction[amplicon]["2"][filename]["matchdistance"]=matchdist
			try:
				#print readlookup[filetocheck[1]]
				if readlookup[filetocheck[1]] == "pass":
					reallytruegood += 1
				else:
					reallymismatch += 1
			except Exception, err:
				#print "Not basecalled", err
				reallynotbasecalled += 1
			
		else:
			if (8 * testlen) > len(event_collection) > (2 * testlen):
				if (amplicon not in readprediction):					readprediction[amplicon]=dict()
					
				if ("3" not in readprediction[amplicon]):
					readprediction[amplicon]["3"]=dict()
					
				if (filename not in readprediction[amplicon]["3"]):
					readprediction[amplicon]["3"][filename]=dict()
					
				readprediction[amplicon]["3"][filename]["name"]=filename
				readprediction[amplicon]["3"][filename]["matchdistance"]=matchdist
				#### A read here is in tier 2 - we don't know if it is any good or not...
				print "Looking promising"
				try:
				#print readlookup[filetocheck[1]]
					if readlookup[filetocheck[1]] == "pass":
						truegood += 1
					else:
						mismatch += 1
				except Exception, err:
					#print "Not basecalled", err
					notbasecalled += 1
			else:
				if (amplicon not in readprediction):					readprediction[amplicon]=dict()
					
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
				try:
					#print readlookup[filetocheck[1]]
					if readlookup[filetocheck[1]] == "fail":
						truebad += 1
					else:
						mismatch += 1
				except Exception, err:
					#print "Not basecalled", err
						notbasecalled += 1

				
				
			print "THIS DOESN'T SATISFY ME"
			reallybad += 1
		
#		if (8 * testlen) > len(event_collection) > (2 * testlen):
#			print "Looking promising"
#			goodcount+=1
#			###If the read looks promising we need to look for the hairpin.
#			dist,cost,path=checkhairpin(event_collection,hairpin3)
#			print dist,path[1][0],path[1][-1],(path[1][-1]-path[1][1])
#			#print path[0]
#			#print path[1]
#			### Now check and see if the complement strand is longer than the template...
#			## The start of the hairpin is path[1][0] above. The end of the hairpin is path[1][-1] - we just want to make sure that the complement length (defined as len(event_collection)-path[1][-1]+1) is longer than the template length (defined as (path[1][0]-1))
			
				
			
#			filetocheck = os.path.split(filename)
#			try:
#				#print readlookup[filetocheck[1]]
#				if readlookup[filetocheck[1]] == "pass":
#					truegood += 1
#				else:
#					mismatch += 1
#			except Exception, err:
#				#print "Not basecalled", err
#					notbasecalled += 1
#		else:
#			print "**********BAD READ************"
#			badcount+=1
#			###Even though the read looks bad lets check it anyway!	
#			dist,cost,path=checkhairpin(event_collection,hairpin3)
#			print dist,path[1][0],path[1][-1],(path[1][-1]-path[1][1])
#			filetocheck = os.path.split(filename)
#			try:
#				#print readlookup[filetocheck[1]]
#				if readlookup[filetocheck[1]] == "fail":
#					truebad += 1
#				else:
#					mismatch += 1
#			except Exception, err:
				#print "Not basecalled", err
#					notbasecalled += 1
		
		### We want to know if this read really is good or bad as defined by nanopore - we already have this in the database so we can check from there!
		filetocheck = os.path.split(filename)
		try:
			print readlookup[filetocheck[1]]
		except Exception, err:
			print "Not basecalled", err
	hdf.close()
	
#print "We have "+str(goodcount)+" potential reads and "+str(badcount)+" bad reads."
#print "We have "+str(truegood)+" true good reads, "+str(truebad)+" true bad reads, "+str(mismatch)+" mismatched reads and "+str(notbasecalled)+" notbasecalled reads."
#print "We have "+str(reallygood)+" really good reads and "+str(reallytruegood)+" really true good and "+str(reallymismatch)+" true mismatches and "+str(reallynotbasecalled)+" really not basecalled reads and "+str(reallybad)+" really bad reads."


#print readprediction

#for amplicon in readprediction:
#	for key in readprediction[amplicon]:
#		print amplicon,key

for amplicon in readprediction:
	print amplicon
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
					if not os.path.exists('prefiltered'):
						os.makedirs('prefiltered')
					destdir = os.path.join('prefiltered','downloads')    	
					if not os.path.exists(destdir):
						os.makedirs(destdir)
					try:
						filetocheck = os.path.split(read)
						subdir=readlookup[filetocheck[1]]
						sourcefile = os.path.join(filetocheck[0],'downloads',subdir,filetocheck[1])
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
		print "No reads of class 1"
	try:
		if (len(readprediction[amplicon]["1"].keys())>0):
			print len(readprediction[amplicon]["1"].keys())
			if (counter < args.depth):			
				ordered1 = OrderedDict(sorted(readprediction[amplicon]["1"].iteritems(), key=lambda x: x[1]['matchdistance']))
				#print ordered1
				#ordered1 = readprediction[amplicon]["1"].iteritems()
				for read in ordered1:
					print read, ordered1[read]["matchdistance"]
					if not os.path.exists('prefiltered'):
						os.makedirs('prefiltered')
					destdir = os.path.join('prefiltered','downloads')    	
					if not os.path.exists(destdir):
						os.makedirs(destdir)
					try:
						filetocheck = os.path.split(read)
						subdir=readlookup[filetocheck[1]]
						sourcefile = os.path.join(filetocheck[0],'downloads',subdir,filetocheck[1])
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
					if not os.path.exists('prefiltered'):
						os.makedirs('prefiltered')
					destdir = os.path.join('prefiltered','downloads')    	
					if not os.path.exists(destdir):
						os.makedirs(destdir)
					try:
						filetocheck = os.path.split(read)
						subdir=readlookup[filetocheck[1]]
						sourcefile = os.path.join(filetocheck[0],'downloads',subdir,filetocheck[1])
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
					if not os.path.exists('prefiltered'):
						os.makedirs('prefiltered')
					destdir = os.path.join('prefiltered','downloads')    	
					if not os.path.exists(destdir):
						os.makedirs(destdir)
					try:
						filetocheck = os.path.split(read)
						subdir=readlookup[filetocheck[1]]
						sourcefile = os.path.join(filetocheck[0],'downloads',subdir,filetocheck[1])
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
					if not os.path.exists('prefiltered'):
						os.makedirs('prefiltered')
					destdir = os.path.join('prefiltered','downloads')    	
					if not os.path.exists(destdir):
						os.makedirs(destdir)
					try:
						filetocheck = os.path.split(read)
						subdir=readlookup[filetocheck[1]]
						sourcefile = os.path.join(filetocheck[0],'downloads',subdir,filetocheck[1])
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