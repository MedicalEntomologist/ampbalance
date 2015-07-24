#!/usr/bin/env python
import sys, os, re
import time
import threading, thread
from Bio import SeqIO
from StringIO import StringIO
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

from dtw2sam import AlignedEvent
from dtw2sam import DtwToSam
from dtw2sam import DtwAlignment

lenRef = 0

parser = configargparse.ArgParser(description='DTWmap - a squiggle based read aligner. Developed by Matt Loose @mattloose or matt.loose@nottingham.ac.uk for help!')
parser.add('-fasta', '--reference_fasta_file', type=str, dest='fasta', required=True, default=None, help="The fasta format file for the reference sequence for your organism.")
parser.add('-w', '--watch-dir', type=str, required=True, default=None, help="The path to the folder containing the downloads directory with fast5 reads to analyse - e.g. C:\data\minion\downloads (for windows).", dest='watchdir')
parser.add('-o', '--output-sam', type=str, required=True, default='out.sam', help="Output filename for sam file.", dest='outsam')
parser.add('-v', '--verbose-true', action='store_true', help="Print detailed messages while processing files.", default=False, dest='verbose')
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
        reader = csv.reader(csv_file, delimiter=",")
        d = list(reader)
        #print d
        for r in range(1, len(d)):
            #print r, d[r]
            kmer = d[r][0]
            mean = d[r][1]
            #print r, kmer, mean
            model_kmers[kmer]=mean
    return  model_kmers



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
            if (args.verbose is True):
                memory_usage_psutil()
            queryarray = sklearn.preprocessing.scale(np.array(squiggle),axis=0,with_mean=True,with_std=True,copy=True)
            dist, cost, path = mlpy.dtw_subsequence(queryarray,hashthang[ref]['Fprime'])
            if (args.verbose is True):
                memory_usage_psutil()
            result.append((dist,ref,"F",path[1][0],path[1][-1],path[0][0],path[0][-1],cost,path))
            dist, cost, path = mlpy.dtw_subsequence(queryarray,hashthang[ref]['Rprime'])
            result.append((dist,ref,"R",(len(hashthang[ref]['Rprime'])-path[1][-1]),(len(hashthang[ref]['Rprime'])-path[1][0]),path[0][0],path[0][-1],cost,path))
            if (args.verbose is True):
                memory_usage_psutil()
        except Exception,err:
            print "Warp Fail"

    return sorted(result,key=lambda result: result[0])[0][1],sorted(result,key=lambda result: result[0])[0][0],sorted(result,key=lambda result: result[0])[0][2],sorted(result,key=lambda result: result[0])[0][3],sorted(result,key=lambda result: result[0])[0][4],sorted(result,key=lambda result: result[0])[0][5],sorted(result,key=lambda result: result[0])[0][6],sorted(result,key=lambda result: result[0])[0][7],sorted(result,key=lambda result: result[0])[0][8]



######################################################



def process_ref_fasta_raw(ref_fasta,model_kmer_means):
    print "processing the reference fasta."
    global lenRef
    kmer_len=5
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

model_file_template = "template_model_5.model"
model_file_complement = "complement_model_5.model"


model_kmer_means_template=process_model_file(model_file_template)
model_kmer_means_complement=process_model_file(model_file_complement)

kmerhashT = process_ref_fasta_raw(fasta_file,model_kmer_means_template)
kmerhashC = process_ref_fasta_raw(fasta_file,model_kmer_means_complement)

seqlengths = get_seq_len(fasta_file)



readprediction=dict()




#SamAlignments=[]
SamOutFile=DtwToSam(path_to_template_model_file=os.path.abspath(model_file_template),path_to_complement_model_file=os.path.abspath(model_file_complement),outfile=args.outsam,rname="EM_079517")


for filename in glob.glob(os.path.join(args.watchdir, '*.fast5')):
    print filename
    if (args.verbose is True):
        print os.path.abspath(filename)
    hdf = h5py.File(filename, 'r')
    #try:
    for read in hdf['Analyses']['EventDetection_000']['Reads']:
        events = hdf['Analyses']['EventDetection_000']['Reads'][read]['Events'][()]
        event_collection=list()
        for event in events:
            print event[2]
            #print type event[2]
            event_collection.append(event[2])


        read_id_fields = ['duration','hairpin_found','hairpin_event_index','read_number','scaling_used','start_mux','start_time','read_id']
        #print "!!!!!!READ IS",read
        read_info_hash =  make_hdf5_object_attr_hash(hdf['Analyses/EventDetection_000/Reads/'+read],read_id_fields)

        
        #print read_info_hash['hairpin_found']
        if read_info_hash['hairpin_found']==1:
            if (args.verbose is True):
                print "!!!!!!!!!!!Hairpin Found!!!!!!!!!!!"
                print "Template Length:", len(event_collection[0:read_info_hash['hairpin_event_index']])
                print "Complement Length:", len(event_collection[read_info_hash['hairpin_event_index']:len(event_collection)])

                print read_info_hash['read_id']
            ####Can we find a sudden shift in current - up and down... Lets look for a sudden shift down and declare that to be the start of the read...

            prev1 = 0
            prev2 = 0
            counter = 0
            Template_Start_Event = 0

            for event in event_collection[0:read_info_hash['hairpin_event_index']]:
                #print event,prev1,prev2
                if (prev1 - event > 40 ) or (prev2 - event > 40):
                    if (args.verbose is True):
                        print "Large Delta", counter
                    Template_Start_Event=counter+1
                    break
                prev2 = prev1
                prev1 = event
                counter +=1

            if (args.verbose is True):
                print "Template Start Event",Template_Start_Event
            
            try:
                (seqmatchnameT,distanceT,frT,rsT,reT,qsT,qeT,costT,pathT) = raw_squiggle_search2(event_collection[Template_Start_Event:read_info_hash['hairpin_event_index']],kmerhashT)
                print "Warp 1 Complete"
            except Exception,err:
                print "A time warp failed:", err

            ####Can we find a sudden shift in current - up and down... Lets look for a sudden shift down and declare that to be the start of the read...

            prev1 = 0
            prev2 = 0
            counter = 0
            deltacount = 0
            Complement_Start_Event = 0

            for event in event_collection[read_info_hash['hairpin_event_index']:len(event_collection)]:
                #print event,prev1,prev2
                if (prev1 - event > 40 ) or (prev2 - event > 40):
                    if (args.verbose is True):
                        print "Large Delta", counter
                    deltacount += 1
                    Complement_Start_Event=counter+1
                    #break
                if deltacount >= 2:
                    break
                prev2 = prev1
                prev1 = event
                counter +=1

            #exit()
            if (args.verbose is True):
                print "Complement Start Event",Complement_Start_Event

            try:
                (seqmatchnameC,distanceC,frC,rsC,reC,qsC,qeC,costC,pathC) = raw_squiggle_search2(event_collection[read_info_hash['hairpin_event_index']+Complement_Start_Event:len(event_collection)],kmerhashC)
                print "Warp 2 Complete"
            except Exception,err:
                print "A time warp failed:", err

            read_mapping_quality = 0
            ### If the forward and reverse reads map appropriately and overlap to the reference we upload template,complement and 2d. But what coordinate do we give for the 2D? Perhaps the overlap region?
            if (seqmatchnameC==seqmatchnameT and frT != frC and reC >= rsT and rsC <= reT):
                print "Good Candidate"
                read_mapping_quality = 1
                if (rsT < rsC):
                    start = rsT
                else:
                    start = rsC
                if (reT > reC):
                    end = reT
                else:
                    end = reC

                if (args.verbose is True):
                    print "PathT"
                    print pathT
                    print "PathC"
                    print pathC
            else:
                print "We didn't satisfy something here"
                read_mapping_quality = 0
                if (args.verbose is True):
                    print "Template",frT,rsT,reT
                    print "Complement",frC,rsC,reC
            #        print "CostT"
            #        print costT
            #        print "CostC"
            #        print costC
                    print "PathT"
                    print pathT
                    print "PathC"
                    print pathC

            if read_mapping_quality == 1:

                if (args.verbose is True):
                    print "PathT"
                f,t = pathT
                dat  = zip(f,t)

                alignedevents=[]
                for(a,b) in dat:
                    readpos = a+Template_Start_Event 
                    if (args.verbose is True):
                        print a+Template_Start_Event,
                    if (frT=="R"):
                        if (args.verbose is True):
                            print (lenRef-b),
                        refpos = (lenRef-b)
                    else:
                        if (args.verbose is True):
                            print b,
                        refpos = b
                    if (args.verbose is True):    
                        print distanceT,seqmatchnameT,frT,event_collection[a+Template_Start_Event],
                    if (frT=="F"):
                        if (args.verbose is True):
                            print "Forward Mapping",
                            print kmerhashT[seqmatchnameT]["F"][b],kmerhashT[seqmatchnameT]["Fprime"][b]
                    else:
                        if (args.verbose is True):
                            print "Reverse Mapping",
                            print kmerhashT[seqmatchnameT]["R"][b],kmerhashT[seqmatchnameT]["Rprime"][b]
                    alignedevents.append(AlignedEvent([readpos,refpos,distanceT,seqmatchnameT,frT]))
                
                leftClip=Template_Start_Event
                rightClip=len(event_collection)-read_info_hash['hairpin_event_index']
                #SamAlignments.append(DtwAlignment(alignedevents,filename = read_info_hash['read_id']+".template",path_to_read = os.path.abspath(filename), start_clipping=leftClip,end_clipping=rightClip))
                SamOutFile.convert(DtwAlignment(alignedevents,filename = read_info_hash['read_id']+".template",path_to_read = os.path.abspath(filename), start_clipping=leftClip,end_clipping=rightClip,))
                if (args.verbose is True):
                    print "pathC"
                f,t = pathC
                dat  = zip(f,t)
                alignedevents=[]
                for(a,b) in dat:
                    readpos = a+read_info_hash['hairpin_event_index']+Complement_Start_Event
                    if (args.verbose is True):
                        print a+read_info_hash['hairpin_event_index']+Complement_Start_Event,
                    if (frC=="R"):
                        if (args.verbose is True):
                            print (lenRef-b),
                        refpos = (lenRef-b)
                    else:
                        if (args.verbose is True):
                            print b,
                        refpos = b
                    if (args.verbose is True):
                        print distanceC,seqmatchnameC,frC,event_collection[a+read_info_hash['hairpin_event_index']+Complement_Start_Event],
                    if (frC=="F"):
                        if (args.verbose is True):
                            print "Forward Mapping",
                            print kmerhashC[seqmatchnameC]["F"][b],kmerhashC[seqmatchnameC]["Fprime"][b]
                    else:
                        if (args.verbose is True):
                            print "Reverse Mapping",
                            print kmerhashC[seqmatchnameC]["R"][b],kmerhashC[seqmatchnameC]["Rprime"][b]
                    leftClip=read_info_hash['hairpin_event_index']+Complement_Start_Event
                    rightClip=0
                    alignedevents.append(AlignedEvent([readpos,refpos,distanceC,seqmatchnameC,frC]))
                SamOutFile.convert(DtwAlignment(alignedevents,filename = read_info_hash['read_id']+".complement",path_to_read = os.path.abspath(filename), start_clipping=leftClip,end_clipping=rightClip))
                        
        else:
            print "!!!!!!!!!!!Hairpin Not Found!!!!!!!!!!!"
    
    hdf.close()

SamOutFile.close()

