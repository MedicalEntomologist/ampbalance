#!/usr/bin/env python
import sys, os, re
import time
import threading, thread
from Bio import SeqIO
from StringIO import StringIO
#import MySQLdb
import string
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
import psutil

lenRef = 0


parser = configargparse.ArgParser(description='DTWmap - a squiggle based read aligner. Developed by Matt Loose @mattloose or matt.loose@nottingham.ac.uk for help!')
parser.add('-fasta', '--reference_fasta_file', type=str, dest='fasta', required=True, default=None, help="The fasta format file for the reference sequence for your organism.")
parser.add('-w', '--watch-dir', type=str, required=True, default=None, help="The path to the folder containing the downloads directory with fast5 reads to analyse - e.g. C:\data\minion\downloads (for windows).", dest='watchdir')
args = parser.parse_args()

###########################################################
def memory_usage_psutil():
    process = psutil.Process(os.getpid())
    pc = round(process.memory_percent(),2)
    print "Mem Used %s: " % (str(pc))

###########################################################
def make_hdf5_object_attr_hash(hdf5object, fields):
	att_hash=dict()
	for field in fields:
		if (field in hdf5object.attrs.keys() ):
			#print "filed: ",field (args.ref_fasta is not None), hdf5object.attrs[field]
			att_hash[field]=hdf5object.attrs[field]
	return att_hash

######################################################
def process_model_file(model_file):
	model_kmers = dict()
	with open(model_file, 'rb') as csv_file:
		reader = csv.reader(csv_file, delimiter="\t")
    		d = list(reader)
		#print d
		for r in range(1, len(d)):
			#print r, d[r]
			kmer = d[r][0]
			mean = d[r][1]
			#print r, kmer, mean
			model_kmers[kmer]=mean
	return 	model_kmers



######################################################
def get_seq_len(ref_fasta):
	seqlens=dict()
	for record in SeqIO.parse(ref_fasta, 'fasta'):
		seq=record.seq
		seqlens[record.id]=len(seq)
	return seqlens


#######################################################################
def raw_squiggle_search2(squiggle,hashthang):
    result=[]

    for ref in hashthang:
        try:
            memory_usage_psutil()
            queryarray = sklearn.preprocessing.scale(np.array(squiggle),axis=0,with_mean=True,with_std=True,copy=True)
            dist, cost, path = mlpy.dtw_subsequence(queryarray,hashthang[ref]['Fprime'])
            memory_usage_psutil()
            result.append((dist,ref,"F",path[1][0],path[1][-1],path[0][0],path[0][-1],cost,path))
            dist, cost, path = mlpy.dtw_subsequence(queryarray,hashthang[ref]['Rprime'])
            result.append((dist,ref,"R",(len(hashthang[ref]['Rprime'])-path[1][-1]),(len(hashthang[ref]['Rprime'])-path[1][0]),path[0][0],path[0][-1],cost,path))
            memory_usage_psutil()
        except Exception,err:
            print "Warp Fail"

	return sorted(result,key=lambda result: result[0])[0][1],sorted(result,key=lambda result: result[0])[0][0],sorted(result,key=lambda result: result[0])[0][2],sorted(result,key=lambda result: result[0])[0][3],sorted(result,key=lambda result: result[0])[0][4],sorted(result,key=lambda result: result[0])[0][5],sorted(result,key=lambda result: result[0])[0][6],sorted(result,key=lambda result: result[0])[0][7],sorted(result,key=lambda result: result[0])[0][8]



######################################################

'''
def reverseIndexVals(x): 
	if x==None: return None
	else: return lenRef - x
'''


def process_ref_fasta_raw(ref_fasta,model_kmer_means):
	global lenRef
	print "processing the reference fasta."
	kmer_len=7
	kmer_means=dict()
	for record in SeqIO.parse(ref_fasta, 'fasta'):
		kmer_means[record.id]=dict()
		kmer_means[record.id]["F"]=list()
		kmer_means[record.id]["R"]=list()
		kmer_means[record.id]["Fprime"]=list()
		kmer_means[record.id]["Rprime"]=list()
		print "ID", record.id
		print "length", len(record.seq)
		print "FORWARD STRAND"

		seq = record.seq
		lenRef = len(seq)
		for x in range(len(seq)+1-kmer_len):
			kmer = str(seq[x:x+kmer_len])
			kmer_means[record.id]["F"].append(float(model_kmer_means[kmer]))
			#if model_kmer_means[kmer]:
				#print x, kmer, model_kmer_means[kmer]

		print "REVERSE STRAND"
		seq = revcomp = record.seq.reverse_complement()
		for x in range(len(seq)+1-kmer_len):
			kmer = str(seq[x:x+kmer_len])
			kmer_means[record.id]["R"].append(float(model_kmer_means[kmer]))

		kmer_means[record.id]["Fprime"]=sklearn.preprocessing.scale(kmer_means[record.id]["F"], axis=0, with_mean=True, with_std=True, copy=True)
		kmer_means[record.id]["Rprime"]=sklearn.preprocessing.scale(kmer_means[record.id]["R"], axis=0, with_mean=True, with_std=True, copy=True)
	return kmer_means
#######################################################################

def processOutput(path, dist):
	f,t = path
	dat  = zip(f,t)
	for (a,b) in dat: h.write(toCSV(dat))

def toCSV(xs): return (','.join(map(str,xs)))+'\n'



fasta_file = args.fasta

model_file_template = "template_a68_2_final_model_501.model"
model_file_complement = "complement_a70_final_model_101.model"
model_kmer_means_template=process_model_file(model_file_template)
model_kmer_means_complement=process_model_file(model_file_complement)

kmerhashT = process_ref_fasta_raw(fasta_file,model_kmer_means_template)
kmerhashC = process_ref_fasta_raw(fasta_file,model_kmer_means_complement)

seqlengths = get_seq_len(fasta_file)


#sys.exit()


###A dictionary to store our results in

readprediction=dict()

print "Now we are going to try and open the raw reads and do the same as we have done above..."
for filename in glob.glob(os.path.join(args.watchdir, '*.fast5')):
    print filename
    hdf = h5py.File(filename, 'r')
    #try:
    for read in hdf['Analyses']['EventDetection_000']['Reads']:
        events = hdf['Analyses']['EventDetection_000']['Reads'][read]['Events'][()]
        event_collection=list()
        for event in events:
            event_collection.append(event[2])


        read_id_fields = ['duration','hairpin_found','hairpin_event_index','read_number','scaling_used','start_mux','start_time',]
        #print "!!!!!!READ IS",read
        read_info_hash =  make_hdf5_object_attr_hash(hdf['Analyses/EventDetection_000/Reads/'+read],read_id_fields)

        #print read_info_hash['hairpin_found']
        if read_info_hash['hairpin_found']==1:
            print "!!!!!!!!!!!Hairpin Found!!!!!!!!!!!"
            print "Template Length:", len(event_collection[0:read_info_hash['hairpin_event_index']])
            print "Complement Length:", len(event_collection[read_info_hash['hairpin_event_index']:len(event_collection)])

            #if len(event_collection[0:read_info_hash['hairpin_event_index']]) >= (0.5 * args.length) and len(event_collection[0:read_info_hash['hairpin_event_index']]) <= (1.5 * args.length) and len(event_collection[read_info_hash['hairpin_event_index']:len(event_collection)]) >= (0.5 * args.length) and len(event_collection[read_info_hash['hairpin_event_index']:len(event_collection)]) <= (1.5 * args.length):
            try:
                (seqmatchnameT,distanceT,frT,rsT,reT,qsT,qeT,costT,pathT) = raw_squiggle_search2(event_collection[0:read_info_hash['hairpin_event_index']],kmerhashT)
                print "Warp 1 Complete"
            except Exception,err:
                print "A time warp failed:", err
	#		print (seqmatchnameT,distanceT,frT,rsT,reT,qsT,qeT)
			#print squiggle_search2(meansquiggle[0:read_info_hash['hairpin_event_index']],kmerhash)
            try:
                (seqmatchnameC,distanceC,frC,rsC,reC,qsC,qeC,costC,pathC) = raw_squiggle_search2(event_collection[read_info_hash['hairpin_event_index']:len(event_collection)],kmerhashC)
                print "Warp 2 Complete"
            except Exception,err:
                print "A time warp failed:", err
            #print (seqmatchnameC,distanceC,frC,rsC,reC,qsC,qeC)
			### If the forward and reverse reads map appropriately and overlap to the reference we upload template,complement and 2d. But what coordinate do we give for the 2D? Perhaps the overlap region?
            if (seqmatchnameC==seqmatchnameT and frT != frC and reC >= rsT and rsC <= reT):
                print "Good Candidate"
                if (rsT < rsC):
                    start = rsT
                else:
                    start = rsC
                if (reT > reC):
                    end = reT
                else:
                    end = reC
                #target = 3.19
#                amplicon, value = min(ampstartdict.items(), key=lambda (_, v): abs(v - start))
                #print min(ampstarttuple, key=lambda y:ampstarttuple[y]-start)
#                print amplicon, value
#                key2, value2 = min(ampenddict.items(), key=lambda (_, v): abs(v - end))
                #print min(ampstarttuple, key=lambda y:ampstarttuple[y]-start)
#                print key2, value2
                print "CostT"
                print costT
                print "CostC"
                print costC
		
                print "PathT"
                print pathT
                print "PathC"
                print pathC
            else:
                print "We didn't satisfy something here"
                print "Template",frT,rsT,reT
                print "Complement",frC,rsC,reC
		print "CostT"
                print costT
                print "CostC"
                print costC
                print "PathT"
                print pathT
        else:
            print "!!!!!!!!!!!Hairpin Not Found!!!!!!!!!!!"

	print lenRef
	'''
	h= open('out.csv','w')
	processOutput(pathT, distanceT)
	processOutput(pathC, distanceC)
	h.close()
	'''
	h= open('out2D.csv','w')

	indT, indRT = pathT
	indC, indRC = pathC
	
	if frT=='F': 
		meansToUseT = 'Fprime'
		datT = zip(indRT, indRT)
	else: 
		meansToUseT = 'Rprime'
		# invertNumbering...
		datT = zip(lenRef - indRT, indT)

	if frC=='F': 
		meansToUseC = 'Fprime'
		datC = zip(indRC, indC)
	else: 
		meansToUseC = 'Rprime'
		# invertNumbering...
		datC = zip(lenRef - indRC, indC)


	for refPos in xrange(lenRef):
		t_vals = filter(lambda (k,_): k==refPos, datT)
		c_vals = filter(lambda (k,_): k==refPos, datC)
		if t_vals==[]: t_vals=[(refPos,None)]
		if c_vals==[]: c_vals=[(refPos,None)]

		for _,t in t_vals:
		  for _,c in c_vals:
		    try: 
			ref_kmermeanT = \
			    kmerhashT[seqmatchnameT][meansToUseT][refPos] 
		    except: 
			ref_kmermeanT = None
		    try: ref_kmermeanC = \
			    kmerhashC[seqmatchnameC][meansToUseC][refPos] 
		    except: 
			ref_kmermeanC = None
		    dat = refPos, t, c, ref_kmermeanT, ref_kmermeanC
		    #print dat
		    h.write(toCSV(dat))


	h.close()
	exit()





     #          for position in pathT[0]:
     #              print position, pathT[1][position]
                
            #else:
            #    print "Reads too long"
    #except Exception, err:
    #    print "Sorry - couldn't process read ",read
    #    print "Any usefull error message is ",err

    hdf.close()



