#!/usr/bin/env python

#!Python
#Usage:
#BP.Merge.Multi.Sample.py --workdir workdir
#For debug use only
#command='SV.Merge.Multi.Sample.py  --workdir /scratch/remills_flux/xuefzhao/SV_discovery_index/download/SVelter.CommonBPs'
#sys.argv=command.split()
import os
import re
import getopt
import sys
import numpy
import scipy
import random
import math
import numpy as np
from numpy import vstack,array
from numpy.random import rand
from scipy.cluster.vq import vq, kmeans, whiten
from scipy import stats
from scipy.stats import linregress
from sklearn import cluster
from scipy.spatial import distance
import sklearn.datasets
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
def path_modify(in_path):
   if not in_path[-1]=='/':
      in_path+='/'
   return in_path

def path_mkdir(in_path):
   if not os.path.isdir(in_path):
      os.system(r'''mkdir %s'''%(in_path))

def BP_in_path_readin(workdir,sample_file):
	fin=open(sample_file)
	samples=[]
	for line in fin:
		pin=line.strip().split()
		if len(pin)>0:
			samples.append(pin[0])
	fin.close()
	out_path=[]
	for x in samples:
		out_path.append(workdir+'BreakPoints.'+x.split('/')[-1])
	return out_path

def BP_info_readin(workdir,sample_file):
	bp_paths=BP_in_path_readin(workdir,sample_file)
	bp_paths=[path_modify(i) for i in bp_paths]
	BP_info={}
	LN_info={}
	for k1 in bp_paths:
		key1=k1.split('/')[-2]
		for k2 in os.listdir(k1):
			if k2.split('.')[-1]=='SPs':
				chromo=k2.replace('.'.join(k1.split('/')[-2].split('.')[1:-1])+'.','').split('.')[0]
				SP_temp=SP_info_readin(k1+k2,chromo)
				if not chromo in BP_info.keys():
					BP_info[chromo]={}
				if not key1 in BP_info[chromo].keys():
					BP_info[chromo][key1]=[]
				BP_info[chromo][key1]=SP_temp
			if k2.split('.')[-1]=='LNs':
				chromo=k2.replace('.'.join(k1.split('/')[-2].split('.')[1:-1])+'.','').split('.')[0]
				LN_temp=LN_info_readin(k1+k2,chromo)
				if not chromo in LN_info.keys():
					LN_info[chromo]={}
				if not key1 in LN_info[chromo].keys():
					LN_info[chromo][key1]=[]
				LN_info[chromo][key1]=LN_temp
	for k1 in BP_info.keys():
		BP_info[k1]['new']={}
		for k2 in BP_info[k1].keys():
			if not k2=='new':
				for k3 in BP_info[k1][k2]:
					if not k3 in BP_info[k1]['new'].keys():
						BP_info[k1]['new'][k3]=1
					else:
						BP_info[k1]['new'][k3]+=1
		for k3 in BP_info[k1]['new'].keys():
			if BP_info[k1]['new'][k3]==1:
				del BP_info[k1]['new'][k3]
		temp=sorted(BP_info[k1]['new'].keys())
		for x1 in temp:
			for x2 in BP_info[k1].keys():
				if not x2=='new':
					if not x1 in BP_info[k1][x2]:
						BP_info[k1][x2].append(x1)
		LN_info_temp={}
		for k2 in LN_info[k1].keys():
			for k3 in LN_info[k1][k2]:
				if not '_'.join([str(i) for i in k3]) in LN_info_temp.keys():
					LN_info_temp['_'.join([str(i) for i in k3])]=1
				else:
					LN_info_temp['_'.join([str(i) for i in k3])]+=1
		for x1 in LN_info_temp.keys():
			if LN_info_temp[x1]==2:
				for x2 in LN_info[k1].keys():
					if not [int(i) for i in x1.split('_')] in LN_info[k1][x2]:
						LN_info[k1][x2].append([int(i) for i in x1.split('_')])
	for k1 in BP_info.keys():
		if 'new' in BP_info[k1].keys():
			del BP_info[k1]['new']
	return [BP_info,LN_info]

def SP_info_readin(SP_file,chromo):
	if not chromo in SP_info_rec.keys():
		SP_info_rec[chromo]={}
	out=[]
	fin=open(SP_file)
	for line in fin:
		pin=line.strip().split()
		if int(pin[1])>1:
			out.append(int(pin[0]))
		if not int(pin[0]) in SP_info_rec[chromo].keys():
			SP_info_rec[chromo][int(pin[0])]=int(pin[1])
		else:
			SP_info_rec[chromo][int(pin[0])]=max([int(pin[1]),SP_info_rec[chromo][int(pin[0])]])
	fin.close()
	out.sort()
	return out

def LN_info_readin(LN_file,chromo):
	if not chromo in SP_info_rec.keys():
		SP_info_rec[chromo]={}
	out={}
	fin=open(LN_file)
	for line in fin:
		pin=line.strip().split()
		if not int(pin[0]) in SP_info_rec[chromo].keys():
			SP_info_rec[chromo][int(pin[0])]=int(pin[1])
		else:
			SP_info_rec[chromo][int(pin[0])]=max([int(pin[1]),SP_info_rec[chromo][int(pin[0])]])
		if not int(pin[2]) in SP_info_rec[chromo].keys():
			SP_info_rec[chromo][int(pin[2])]=int(pin[3])
		else:
			SP_info_rec[chromo][int(pin[2])]=max([int(pin[3]),SP_info_rec[chromo][int(pin[2])]])
		pin2=sorted([int(pin[0]),int(pin[2])])
		if not pin2[0] in out.keys():
			out[pin2[0]]=[]
		if not pin2[1] in out[pin2[0]]:
			out[pin2[0]].append(pin2[1])
	fin.close()
	out2=[]
	for k1 in sorted(out.keys()):
		for k2 in sorted(out[k1]):
			out2.append([k1,k2])
	return out2

def SP_file_readin(workdir,sample_file):
	SP_files={}
	LN_files={}
	bp_paths=BP_in_path_readin(workdir,sample_file)
	for k1 in bp_paths:
		k1+='/'
		SP_files[k1.split('/')[-2]]={}
		LN_files[k1.split('/')[-2]]={}
		for k2 in os.listdir(k1):
			if k2.split('.')[-1]=='SPs':
				chromo=k2.replace('.'.join(k1.split('/')[-2].split('.')[1:-1])+'.','').split('.')[0]
				if not chromo in SP_files[k1.split('/')[-2]].keys():
					SP_files[k1.split('/')[-2]][chromo]=[]
				SP_files[k1.split('/')[-2]][chromo].append(k2)
			if k2.split('.')[-1]=='LNs':
				chromo=k2.replace('.'.join(k1.split('/')[-2].split('.')[1:-1])+'.','').split('.')[0]
				if not chromo in LN_files[k1.split('/')[-2]].keys():
					LN_files[k1.split('/')[-2]][chromo]=[]
				LN_files[k1.split('/')[-2]][chromo].append(k2)					
	return [SP_files,LN_files]

def SP_info_write(SP_info):
	for k1 in SP_info.keys():
		for k2 in SP_info[k1].keys():
			fileout=workdir+k2+'/'+files_names[0][k2][k1][0]+'.new'
			fo=open(fileout,'w')
			for ka in sorted(SP_info[k1][k2]):
				print >>fo, ' '.join([str(i) for i in [ka,SP_info_rec[k1][ka]]])
			fo.close()

def LN_info_write(LN_info):
	for k1 in LN_info.keys():
		for k2 in LN_info[k1].keys():
			fileout=workdir+k2+'/'+files_names[1][k2][k1][0]+'.new'
			fo=open(fileout,'w')
			for ka in LN_info[k1][k2]:
				print >>fo, ' '.join([str(i) for i in [ka[0],SP_info_rec[k1][ka[0]],ka[1],SP_info_rec[k1][ka[1]]]])
			fo.close()

def sample_name_readin():
	fin=open(sample_file)
	global sample_name
	sample_name=[]
	for line in fin:
		pin=line.strip().split()
		if not pin=='':
			sample_name+=pin
	fin.close()

def Apply_SVelter_BPIntegrate():
	for k1 in sample_name:
		os.system('SVelter.py BPIntegrate --workdir %s --sample %s --batch %s'%(workdir, k1,'0'))

def bp_files_readin():
	#this function readin only bp clusters on the same chromosome for now. will be modified later for multi-chrom events
	bp_paths=[]
	for k1 in sample_name:
		bp_paths.append(workdir+'bp_files.'+k1.split('/')[-1])
	global bp_files_names
	bp_files_names={}
	bp_info={}
	for k1 in bp_paths:
		temp_key1='.'.join(k1.split('/')[-1].split('.')[1:])
		bp_info[temp_key1]={}
		bp_files_names[temp_key1]={}
		k1=path_modify(k1)
		for k2 in os.listdir(k1):
			temp_path2=k1+k2+'/'
			for k3 in os.listdir(temp_path2):
				temp_path3=temp_path2+k3+'/'
				for k4 in os.listdir(temp_path3):
					if not k4.split('.')[-2] in bp_files_names[temp_key1].keys():
						if not k4.split('.')[-2]=='LN':
							bp_files_names[temp_key1][k4.split('.')[-2]]=temp_path3+k4
					fin=open(temp_path3+k4)
					for line in fin:
						pin=line.strip().split()
						if len(pin)>1:
							chromo_num=chromosome_num_check(pin)
							if chromo_num==1:
								if not pin[0] in bp_info[temp_key1].keys():
									bp_info[temp_key1][pin[0]]={}
								if not int(pin[1]) in bp_info[temp_key1][pin[0]].keys():
									bp_info[temp_key1][pin[0]][int(pin[1])]=[]
								bp_info[temp_key1][pin[0]][int(pin[1])].append([int(i) for i in pin[1:]])
					fin.close()
	return bp_info

def bp_info_reorder():
	bp_info=bp_files_readin()
	out={}
	for k1 in bp_info.keys():
		out[k1]={}
		for k2 in bp_info[k1].keys():
			out[k1][k2]=[]
			for k3 in sorted(bp_info[k1][k2].keys()):
				for k4 in bp_info[k1][k2][k3]:
					if not k4 in out[k1][k2]:
						out[k1][k2].append(k4)
	out2={}
	for k1 in bp_info.keys():
		for k2 in bp_info[k1].keys():
			if not k2 in out2.keys():
				out2[k2]={}
			for k3 in bp_info[k1][k2].keys():
				if not k3 in out2[k2].keys():
					out2[k2][k3]=[]
				for k4 in bp_info[k1][k2][k3]:
					if not k4 in out2[k2][k3]:
						out2[k2][k3].append(k4)
	return [out,out2]

def bp_files_compare():
	bp_info_all=bp_info_reorder()	
	bp_inf=bp_info_all[0]
	total_info=bp_info_all[1]
	common_info={}
	unique_info={}
	for k1 in total_info.keys():
		for k2 in total_info[k1].keys():
			for k3 in total_info[k1][k2]:
				temp_score=0
				for k4 in bp_inf.keys():
					if k1 in bp_inf[k4].keys():
						if k3 in bp_inf[k4][k1]:
							temp_score+=1
				if temp_score>1:
					if not k1 in common_info.keys():
						common_info[k1]=[]
					common_info[k1].append(k3)
	for k1 in bp_inf.keys():
		unique_info[k1]={}
		for k2 in bp_inf[k1].keys():
			for k3 in bp_inf[k1][k2]:
				if k2 in common_info.keys():
					if k3 in common_info[k2]:
						continue
					else:
						if not k2 in unique_info[k1].keys():
							unique_info[k1][k2]=[]
						unique_info[k1][k2].append(k3)
	return [common_info,unique_info]

def bp_files_write():
	write_info=bp_files_compare()
	common_info=write_info[0]
	unique_info=write_info[1]
	for k1 in common_info.keys():
		for k2 in bp_files_names.keys():
			if k1 in bp_files_names[k2].keys():
				k4=bp_files_names[k2][k1]
				fo=open('.'.join(k4.split('.')[:-1]+['integrated.txt']),'w')
				for k5 in common_info[k1]:
					print >>fo, ' '.join([str(i) for i in [k1]+k5])
					print >>fo, ' '
					print ' '.join([str(i) for i in [k1]+k5])
				fo.close()
	for k1 in unique_info.keys():
		for k2 in unique_info[k1].keys():
			if k1 in bp_files_names.keys():
				if k2 in bp_files_names[k1].keys():
					k4=bp_files_names[k1][k2]
					if os.path.isfile('.'.join(k4.split('.')[:-1]+['integrated.txt'])):
						fo=open('.'.join(k4.split('.')[:-1]+['integrated.txt']),'a')
					else:
						fo=open('.'.join(k4.split('.')[:-1]+['integrated.txt']),'w')
					for k3 in unique_info[k1][k2]:
						print >>fo, ' '.join([str(i) for i in [k2]+k3])
						print >>fo, ' '
						print ' '.join([str(i) for i in [k2]+k3])
					fo.close()

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def chromosome_names_readin(chromosome_name):
	ref_in=chromosome_name
	chromosome_names=[]
	fin=open(ref_in+'.fai')
	for line in fin:
		pin=line.strip().split()
		chromosome_names.append(pin[0])
	return chromosome_names

def chromosome_num_check(pin):
	temp_num=0
	for x in pin:
		if x in chromosome_names:
			temp_num+=1
		elif is_number(x)==False:
			temp_num+=1
	return temp_num

def read_sample_in(workdir):
	sample_name={}
	for k1 in os.listdir(workdir):
		if k1.split('.')[0]=='BreakPoints' and os.path.isdir(workdir+k1):
			key_name='.'.join(k1.split('.')[1:])
			if not key_name in sample_name.keys():
				sample_name[key_name]={}
			for k2 in os.listdir(workdir+k1):
				if k2.split('.')[-1] in ['LNs','SPs']:
					key_chr=k2.split('.')[-6]
					if not key_chr in sample_name[key_name].keys():
						sample_name[key_name][key_chr]=[]
					sample_name[key_name][key_chr].append(workdir+k1+'/'+k2)
	for k1 in sample_name.keys():
		for k2 in sample_name[k1].keys():
			sample_name[k1][k2].sort()
	return sample_name

def read_info_in(workdir,chromo_name):
	out={}
	sample_name=read_sample_in(workdir)
	for k1 in sample_name.keys():
		out[k1]={}
		for k2 in sample_name[k1].keys():
			if k2==chromo_name:
				print [k1,k2]
				out[k1][k2]={}
				temp_SP={}
				for k3 in sample_name[k1][k2]:
					out[k1][k2][k3]={}
					if k3.split('.')[-1]=='SPs':
						ka=k3
						fin=open(k3)
						for line in fin:
							pin=line.strip().split()
							if not int(pin[0]) in out[k1][k2][k3].keys():
								out[k1][k2][k3][int(pin[0])]=0
							if int(pin[1])>out[k1][k2][k3][int(pin[0])]:
								out[k1][k2][k3][int(pin[0])]=int(pin[1])
						fin.close()
					elif k3.split('.')[-1]=='LNs':
						kb=k3
						fin=open(k3)
						for line in fin:
							pin=line.strip().split()
							if pin[0]==pin[2]:
								if not int(pin[0]) in temp_SP.keys():
									temp_SP[int(pin[0])]=0
								if int(pin[1]) > temp_SP[int(pin[0])]:
									temp_SP[int(pin[0])]=int(pin[1])
							else:
								if not int(pin[0]) in out[k1][k2][k3].keys():
									out[k1][k2][k3][int(pin[0])]=[]
								if not int(pin[2]) in out[k1][k2][k3][int(pin[0])]:
									out[k1][k2][k3][int(pin[0])].append(int(pin[2]))
								if not int(pin[0]) in temp_SP.keys():
									temp_SP[int(pin[0])]=0
								if int(pin[1]) > temp_SP[int(pin[0])]:
									temp_SP[int(pin[0])]=int(pin[1])
								if not int(pin[2]) in temp_SP.keys():
									temp_SP[int(pin[2])]=0
								if int(pin[3]) > temp_SP[int(pin[2])]:
									temp_SP[int(pin[2])]=int(pin[3])
						fin.close()
				for x in temp_SP.keys():
					if not x in out[k1][k2][ka].keys():
						out[k1][k2][ka][x]=0
					if temp_SP[x]>out[k1][k2][ka][x]:
						out[k1][k2][ka][x]=temp_SP[x]
	out2={}
	for k1 in out.keys():
		for k2 in out[k1].keys():
			if k2 == chromo_name:
				if not k2 in out2.keys():
					out2[k2]={}
				if not k1 in out2[k2].keys():
					out2[k2][k1]={}
				out2[k2][k1]=out[k1][k2]
	return out2

def sample_index(sample_names):
	out={}
	rec=0
	for k1 in sample_names:
		if not k1 in out.keys():
			rec+=1
			out[k1]=rec
	return out

def hash_reverse(hash):
	out={}
	for x in hash.keys():
		out[hash[x]]=x
	return out

def cluster_subsplit(list,dis_cff):
	#eg of list: [247042558, 247042569, 247042616]
	if max(list)-min(list)>dis_cff:
		dis_list=[]
		for x in range(len(list)-1):
			dis_list.append(list[x+1]-list[x])
		dis_index=dis_list.index(max(dis_list))
		out=[list[:dis_index+1],list[dis_index+1:]]
		return out
	else:
		return [list]

def culster_split_recursive(list,dis_cff):
	out=[]
	temp=cluster_subsplit(list,dis_cff)
	if len(temp)==1:
		out=temp
	else:
		for i in temp:
			out+=culster_split_recursive(i,dis_cff)
	return out

def cluster_num(all_SP,dis_cff):
	all_keys=sorted(all_SP.keys())
	out_keys=[[]]
	for k1 in all_keys:
		if out_keys[-1]==[]:
			out_keys[-1].append(k1)
		else:
			if k1-out_keys[-1][-1]<dis_cff:
				out_keys[-1].append(k1)
			else:
				out_keys.append([k1])
	out_keys_2=[]
	for k1 in out_keys:
		out_keys_2+=culster_split_recursive(k1,dis_cff)
	return out_keys_2

def pick_largest_support(all_SP_item):
	#eg pf all_SP_item:[[0, 2]] 
	out=0
	for x in all_SP_item:
		if x[0]>out:
			out=x[0]
	return out

def pick_unique_from_list(list):
	out=[]
	for x in list:
		if not x in out:
			out.append(x)
	return out

def group_SP(clu_SP,all_SP):
	out=[]
	out2=[]
	for x in clu_SP:
		temp=[]
		temp2=[]
		for y in x:
			temp.append(pick_largest_support(all_SP[y]))
			temp2+=[i[1] for i in all_SP[y]]
		out.append(x[temp.index(max(temp))])
		out2.append(pick_unique_from_list(temp2))
	return [out,out2]

def ref_region_readin(chromo,start,end):
	fref=os.popen(r'''samtools faidx %s %s:%d-%d'''%(reference,chromo,int(start),int(end)))
	fref.readline().strip().split()
	seq=''
	while True:
			pref=fref.readline().strip().split()
			if not pref: break
			seq+=pref[0]
	fref.close()
	return seq

def formalize_ref_seq(seq):
	#transfer R,Y,S,W,K,M,B,D,H,V,.,- to N  
	seq2=seq
	if 'R' in seq:
		seq2=seq2.replace('R','N')
	if 'Y' in seq:
		seq2=seq2.replace('Y','N')
	if 'S' in seq:
		seq2=seq2.replace('S','N')
	if 'W' in seq:
		seq2=seq2.replace('W','N')
	if 'K' in seq:
		seq2=seq2.replace('K','N')
	if 'M' in seq:
		seq2=seq2.replace('M','N')
	if 'B' in seq:
		seq2=seq2.replace('B','N')
	if 'D' in seq:
		seq2=seq2.replace('D','N')
	if 'H' in seq:
		seq2=seq2.replace('H','N')
	if 'V' in seq:
		seq2=seq2.replace('V','N')
	if '.' in seq:
		seq2=seq2.replace('.','N')
	if '-' in seq:
		seq2=seq2.replace('-','N')
	return seq2

def dotdata(kmerlen,seq1, seq2):
	if len(sys.argv) < 4:
		sys.exit(1)
	nth_base = 1
	inversions = True
	seq1a=formalize_ref_seq(seq1)
	seq2a=formalize_ref_seq(seq2)
	hits = kmerhits(seq1a, seq2a, kmerlen, nth_base, inversions)
	return hits

def kmerhits(seq1, seq2, kmerlen, nth_base=1, inversions=False):
	# hash table for finding hits
	lookup = {}
	# store sequence hashes in hash table
	#print "hashing seq1..."
	seq1len = len(seq1)
	for i in xrange(seq1len - kmerlen + 1):
		key = seq1[i:i+kmerlen]
		for subkey in subkeys(key, nth_base, inversions):
			lookup.setdefault(subkey, []).append(i)
	# match every nth base by 
	# look up hashes in hash table
	#print "hashing seq2..."
	hits = []
	for i in xrange(len(seq2) - kmerlen + 1):
		key = seq2[i:i+kmerlen]
		# only need to specify inversions for one seq
		for subkey in subkeys(key, nth_base, False):
			subhits = lookup.get(subkey, [])
			if subhits != []:
				# store hits to hits list
				for hit in subhits:
					hits.append((i, hit))
				# break out of loop to avoid doubly counting
				# exact matches
				break
	return hits

def check_ref_qual(chromo,start,end):
	ref_seq=ref_region_readin(chromo,start,end)
	hits=dotdata(window_size,ref_seq,ref_seq)
	region_QC=qual_check_repetitive_region(hits)
	return region_QC

def check_ref_qual_write(chromo,start,end):
	ref_seq=ref_region_readin(chromo,start,end)
	hits=dotdata(window_size,ref_seq,ref_seq)
	fo=open('/nfs/remills-data/xuefzhao/test.txt','w')
	for x in hits:
		print >>fo, ' '.join([str(i) for i in x])
	fo.close()

def qual_check_repetitive_region(dotdata_qual_check):
	#qual_check_R_2(dotdata_qual_check,txt_file)  #not necessary for validator
	diagnal=0
	other=[[],[]]
	for x in dotdata_qual_check:
		if x[0]==x[1]:
			diagnal+=1 
		else:
			if x[0]>x[1]:
				other[0].append(x[0])
				other[1].append(x[1])
	if float(len(other[0]))/float(len(dotdata_qual_check))>0.1 and float(len(other[0]))/float(len(dotdata_qual_check))<0.5:
		other_cluster=X_means_cluster_reformat(other)
		range_cluster=cluster_range_decide(other_cluster)
		size_cluster=cluster_size_decide(range_cluster)
	else:
		size_cluster=[0]
	return [float(diagnal)/float(len(dotdata_qual_check)),size_cluster]

def classify_SP_based_on_samples(clu_SP,all_SP):
	gro_SP=group_SP(clu_SP,all_SP)
	out={}
	for x in range(len(gro_SP[1])):
		if not len(gro_SP[1][x]) in out.keys():
			out[len(gro_SP[1][x])]=[]
		out[len(gro_SP[1][x])].append(gro_SP[0][x])
	for x in sorted(out.keys()):
		print [x,len(out[x])]
	return out

def reverse(seq):
	return seq[::-1]

def subkeys(key, nth_base, inversions):
	subkeys_info = []
	keylen = len(key)
	# speed tip from:
	# http://wiki.python.org/moin/PythonSpeed/PerformanceTips#String_Concatenation
	if nth_base == 1:
		subkeys_info = [key]
	elif nth_base != 0:
		for k in range(nth_base):
			substr_list = [key[j] for j in range(keylen) if (j % nth_base == k)]
			subkeys_info.append("".join(substr_list))
	else:
		# nth_base = 0 is a special case for third base mismatches
		# for every codon, only include the first 2 bases in the hash
		subkeys_info = ["".join([key[i] for i in range(len(key)) if i % 3 != 2])]
	if inversions:
		for i in range(len(subkeys_info)):
			subkeys_info.append("".join([invert_base[c] for c in reversed(subkeys_info[i])]))
	return subkeys_info

def k_means_cluster(data_list):
	#print data_list
	array_diagnal=array([[data_list[0][x],data_list[1][x]] for x in range(len(data_list[0]))])
	ks = range(1,min([5,len(data_list[0])+1]))
	KMeans = [cluster.KMeans(n_clusters = i, init="k-means++").fit(array_diagnal) for i in ks]
	BIC = [compute_bic(kmeansi,array_diagnal) for kmeansi in KMeans]
	ks_picked=ks[BIC.index(max(BIC))]
	if ks_picked==1:
		return [data_list]
	else:
		out=[]
		std_rec=[scipy.std(data_list[0]),scipy.std(data_list[1])]
		whitened = whiten(array_diagnal)
		centroids, distortion=kmeans(whitened,ks_picked)
		idx,_= vq(whitened,centroids)
		for x in range(ks_picked):
			group1=[[int(i) for i in array_diagnal[idx==x,0]],[int(i) for i in array_diagnal[idx==x,1]]]
			out.append(group1)
		return out

def k_means_cluster(data_list):
	if max(data_list[0])-min(data_list[0])>10 and max(data_list[1])-min(data_list[1])>10:
		array_diagnal=array([[data_list[0][x],data_list[1][x]] for x in range(len(data_list[0]))])
		ks = range(1,min([5,len(data_list[0])+1]))
		KMeans = [cluster.KMeans(n_clusters = i, init="k-means++").fit(array_diagnal) for i in ks]
		KMeans_predict=[cluster.KMeans(n_clusters = i, init="k-means++").fit_predict(array_diagnal) for i in ks]
		BIC=[]
		BIC_rec=[]
		for x in ks:
			if KMeans_predict[x-1].max()<x-1: continue
			else:
				BIC_i=compute_bic(KMeans[x-1],array_diagnal)
				if abs(BIC_i)<10**8:
					BIC.append(BIC_i)
					BIC_rec.append(x)
		#BIC = [compute_bic(kmeansi,array_diagnal) for kmeansi in KMeans]
		#ks_picked=ks[BIC.index(max(BIC))]
		ks_picked=BIC_rec[BIC.index(max(BIC))]
		if ks_picked==1:
			return [data_list]
		else:
			out=[]
			std_rec=[scipy.std(data_list[0]),scipy.std(data_list[1])]
			whitened = whiten(array_diagnal)
			centroids, distortion=kmeans(whitened,ks_picked)
			idx,_= vq(whitened,centroids)
			for x in range(ks_picked):
				group1=[[int(i) for i in array_diagnal[idx==x,0]],[int(i) for i in array_diagnal[idx==x,1]]]
				out.append(group1)
			return out
	else:
		return [data_list]

def k_means_cluster_Predict(data_list,info):
	array_diagnal=array([[data_list[0][x],data_list[1][x]] for x in range(len(data_list[0]))])
	ks = range(1,len(info))
	KMeans = [cluster.KMeans(n_clusters = i, init="k-means++").fit(array_diagnal) for i in ks]
	BIC = [compute_bic(kmeansi,array_diagnal) for kmeansi in KMeans]
	ks_picked=ks[BIC.index(max(BIC))]
	if ks_picked==1:
		return [data_list]
	else:
		out=[]
		std_rec=[scipy.std(data_list[0]),scipy.std(data_list[1])]
		whitened = whiten(array_diagnal)
		centroids, distortion=kmeans(whitened,ks_picked)
		idx,_= vq(whitened,centroids)
		for x in range(ks_picked):
			group1=[[int(i) for i in array_diagnal[idx==x,0]],[int(i) for i in array_diagnal[idx==x,1]]]
			out.append(group1)
		return out

def X_means_cluster(data_list):
    temp_result=[i for i in k_means_cluster(data_list) if not i==[[],[]]]
    if temp_result==[data_list]:
        return temp_result[0]
    else:
        out=[]
        for i in temp_result:
            out+=X_means_cluster(i)
        return out

def X_means_cluster_reformat(data_list):
	out=X_means_cluster(data_list)
	out2=[]
	for y in range(len(out)/2):
		out2.append([out[2*y],out[2*y+1]])
	return out2

def compute_bic(kmeans,X):
	"""
	Computes the BIC metric for a given clusters
	Parameters:
	-----------------------------------------
	kmeans:  List of clustering object from scikit learn
	X     :  multidimension np array of data points
	Returns:
	-----------------------------------------
	BIC value
	"""
	# assign centers and labels
	centers = [kmeans.cluster_centers_]
	labels  = kmeans.labels_
	#number of clusters
	m = kmeans.n_clusters
	# size of the clusters
	n = np.bincount(labels)
	#size of data set
	N, d = X.shape
	#compute variance for all clusters beforehand
	cl_var=[]
	for i in xrange(m):
		if not n[i] - m==0:
			cl_var.append((1.0 / (n[i] - m)) * sum(distance.cdist(X[np.where(labels == i)], [centers[0][i]], 'euclidean')**2))
		else:
			cl_var.append(float(10**20) * sum(distance.cdist(X[np.where(labels == i)], [centers[0][i]], 'euclidean')**2))
	const_term = 0.5 * m * np.log10(N)
	BIC = np.sum([n[i] * np.log10(n[i]) -
	       n[i] * np.log10(N) -
	     ((n[i] * d) / 2) * np.log10(2*np.pi) -
	      (n[i] / 2) * np.log10(cl_var[i]) -
	     ((n[i] - m) / 2) for i in xrange(m)]) - const_term
	return(BIC)

def cluster_range_decide(other_cluster):
	out=[]
	for x in other_cluster:
		out.append([])
		out[-1].append([min(x[0]),max(x[0])])
		out[-1].append([min(x[1]),max(x[1])])
	return out

def cluster_simple_dis(y):
	temp=tranform_diagnal_to_horizonal(y)
	temp_hash=list_to_hash(temp)
	temp[0].sort()
	temp2=[temp[0][i+1]-temp[0][i] for i in range(len(temp[0])-1)]
	temp2_index=[i for i in range(len(temp2)) if temp2[i]>simple_dis_cff]
	if not temp2_index==[]:
		temp2_index=[-1]+temp2_index
		temp2_new=[temp[0][temp2_index[i]+1:temp2_index[i+1]+1] for i in range(len(temp2_index)-1)]
		temp2_new.append(temp[0][(temp2_index[-1]+1):])
		temp3_new=[]
		for y in temp2_new:
			temp3_new.append([])
			y2=list_unique(y)
			for z in y2:
				temp3_new[-1]+=temp_hash[z]
		temp4_new=[tranform_horizonal_to_diagnal([temp2_new[i],temp3_new[i]]) for i in range(len(temp2_new))]
		return temp4_new
	else:
		return [y]

def cluter_to_diagnal(out):
    #search for clusters according to diagnal. based on dis of a dot to diagnal:
    out2=tranform_diagnal_to_horizonal(out)
    out2_hash1={}
    for k1 in range(len(out2[1])):
        if not out2[1][k1] in out2_hash1.keys():
            out2_hash1[out2[1][k1]]=[]
        out2_hash1[out2[1][k1]].append(out2[0][k1])
    cluster_a=cluster_numbers(out2_hash1.keys(),clu_dis_cff)
    cluster_b=[]
    for x in cluster_a:
        cluster_b.append([])
        for y in x:
            cluster_b[-1]+=out2_hash1[y]
    cluster_2_a=[]
    cluster_2_b=[]
    cluster_2_rest=[[],[]]
    for x in range(len(cluster_b)):
        if max(cluster_b[x])-min(cluster_b[x])>lenght_cff and len(cluster_b[x])>dots_num_cff:
            cluster_2_a.append([])
            for y in cluster_a[x]:
                cluster_2_a[-1]+=[y for i in out2_hash1[y]]
            cluster_2_b.append(cluster_b[x])
        else:
            cluster_2_rest[0]+=cluster_a[x]
            cluster_2_rest[1]+=cluster_b[x]
    diagnal_segs=[]
    for x in range(len(cluster_2_a)):
        diagnal_segs.append(tranform_horizonal_to_diagnal([cluster_2_b[x],cluster_2_a[x]]))

def cluster_numbers(list,dis_cff):
    out=[[]]
    for k1 in sorted(list):
        if out[-1]==[]:
            out[-1].append(k1)
        else:
            if k1-out[-1][-1]<dis_cff:
                out[-1].append(k1)
            else:
                out.append([k1])
    return out

def cluster_subgroup(cluster_a,point_dis_cff):
    out=[[]]
    for x in sorted(cluster_a):
        if out[-1]==[]:
            out[-1].append(x)
        else:
            if x-out[-1][-1]<point_dis_cff:
                out[-1].append(x)
            else:
                out.append([x])
    return out

def cluster_check(cluster_a,cluster_b):
    out=[[],[]]
    rec=-1
    for x in cluster_b:
        rec+=1
        temp_out=cluster_subgroup(x,point_dis_cff)
        if len(temp_out)==1:
            out[0].append(cluster_a[rec])
            out[0].append()

def cluster_dis_to_diagnal(out):
	out1=tranform_diagnal_to_distance(out)
	out1_hash=list_to_hash(out1)
	cluster_a=cluster_numbers(out1[0],clu_dis_cff)
	out_clustered=[]
	out_left=[[],[]]
	for x in cluster_a:
		if len(x)>dots_num_cff:
			x2=list_unique(x)
			out_clustered.append([[],[]])
			for y in x2:
				for z in out1_hash[y]:
					out_clustered[-1][0].append(out[0][z])
					out_clustered[-1][1].append(out[1][z])
		else:
			x2=list_unique(x)
			for y in x2:
				for z in out1_hash[y]:
					out_left[0].append(out[0][z])
					out_left[1].append(out[1][z])
	out2=tranform_anti_diagnal_to_distance(out_left)			
	out2_hash=list_to_hash(out2)
	cluster_b=cluster_numbers(out2[0],clu_dis_cff)
	out_anti_diag=[]
	for x in cluster_b:
		if len(x)>dots_num_cff:
			x2=list_unique(x)
			out_anti_diag.append([[],[]])
			for y in x2:
				for z in out2_hash[y]:
					out_anti_diag[-1][0].append(out_left[0][z])
					out_anti_diag[-1][1].append(out_left[1][z])
	return [out_clustered,out_anti_diag]

def cluster_dots_based_simplely_on_dis(data_list,simple_dis_cff):
	#subgroups dots based on their distance on x-axil first, and then on y-axil. distance >50 would be used as the cutoff.
	out=[]
	for x in data_list:
		y=k_means_cluster_Predict(x,info)
		out+=y
	return out

def cluster_range_decide(other_cluster):
	out=[]
	for x in other_cluster:
		out.append([])
		out[-1].append([min(x[0]),max(x[0])])
		out[-1].append([min(x[1]),max(x[1])])
	return out

def cluster_size_decide(range_cluster):
	out=[]
	for x in range_cluster:
		size=(x[0][1]-x[0][0])*(x[1][1]-x[1][0])
		out.append(numpy.sqrt(size))
	return out

global window_size
window_size=10
global invert_base
invert_base = { 'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C','N' : 'N','a' : 't', 't' : 'a', 'c' : 'g', 'g' : 'c','n' : 'n'}
def classify_LN_based_on_samples(clu_LN,all_LN):
	out=[]
	for x in clu_LN:
		temp=[]
		for y in x:
			sample_names=[name_index_reverse[i] for i in all_LN[y]]
			for sample_name_i in sample_names:
				for keys_3 in BP_info[chromo][sample_name_i].keys():
					if 'LNs' in keys_3:
						if y in BP_info[chromo][sample_name_i][keys_3].keys():
							temp+=BP_info[chromo][sample_name_i][keys_3][y]
		out.append([int(numpy.median(x)),int(numpy.median(temp))])
	out2=[out[0]]
	for x in out[1:]:
		if not x[1]-x[0]>10**6:
			if x[0]<max(out2[-1])+100:
				out2[-1]+=x
			else:
				out2.append(x)
	out3=[]
	for x in out2:
		x.sort()
		x2=pick_unique_from_list(x)
		#y=[x2[i+1]-x2[i] for i in range(len(x2)-1)]
		if not x2 in out3:
			out3.append(x2)
	return out3

def Associate_SP_LN(class_SP,class_LN):
	temp_all_SP=[]
	temp_all_LN=[]
	for x in class_SP.keys():
		temp_all_SP+=class_SP[x]
	for x in class_LN:
		if len(x)==1:
			temp_all_SP+=x
		else:
			temp_all_LN.append(x)
	temp_all_SP.sort()
	rec_SP=0
	rec_LN=0
	included_SP=[]
	while True:
		if rec_SP==len(temp_all_SP) or rec_LN==len(temp_all_LN): break
		if temp_all_SP[rec_SP]>temp_all_LN[rec_LN][-1]+1000:
			rec_LN+=1
		elif temp_all_SP[rec_SP]>temp_all_LN[rec_LN][0]-1001 and temp_all_SP[rec_SP]<temp_all_LN[rec_LN][-1]+1001:
			temp_all_LN[rec_LN].append(temp_all_SP[rec_SP])
			included_SP.append(temp_all_SP[rec_SP])
			rec_SP+=1
		elif temp_all_SP[rec_SP]<temp_all_LN[rec_LN][0]-1000:
			rec_SP+=1
	temp_out_LN=[]
	for x in temp_all_LN:
		x.sort()
		x2=pick_unique_from_list(x)
		y=[x2[i+1]-x2[i] for i in range(len(x2)-1)]
		x3=[[x2[0]]]
		for z in range(len(y)):
			if y[z]<50:
				x3[-1].append(x2[z+1])
			else:
				x3.append([x2[z+1]])
		x4=[]
		for z in x3:
			if len(z)>1:
				temp=[]
				for z2 in z:
					if z2 in included_SP:
						temp.append(z2)
				if temp==[]:
					z4=[z[int(numpy.median(range(len(z))))]]
				elif len(temp)==1:
					z4=temp
				else:
					z3=[abs(i-z[int(numpy.median(range(len(z))))]) for i in temp]
					z4=[temp[z3.index(min(z3))]]
			else:
				z4=z
			x4.append(z4)
		x5=[i[0] for i in x4]
		temp_out_LN.append(x5)
	temp_out_SP=[]
	for x in temp_all_SP:
		if not x in included_SP:
			temp_out_SP.append(x)
	temp_out_SP2=[i for i in subgroup_list_base_on_distance(temp_out_SP,1000) if len(i)>1]
	return temp_out_LN+temp_out_SP2

def QC_SP(clu_SP):
	out=[]
	for x in clu_SP:
		if len(x)==1 or max(x)-min(x)<100:
			temp_QC=check_ref_qual(chromo,x[0]-500,x[0]+500)
			if temp_QC[0]>0.4 or sum(temp_QC[1])/1000<0.3:
				out.append(x)
			else:
				print x+temp_QC
		else:
			temp_QC=check_ref_qual(chromo,x[0],x[-1])
			if temp_QC[0]>0.4 or sum(temp_QC[1])/float((x[-1]-x[0]))<0.3:
				out.append(x)
			else:
				print x+temp_QC
	return out

def subgroup_list_base_on_distance(temp_out_SP,dis):
	temp_gap_SP=[temp_out_SP[i+1]-temp_out_SP[i] for i in range(len(temp_out_SP)-1)]
	temp_out_SP2=[[temp_out_SP[0]]]
	for x in range(len(temp_gap_SP)):
		if temp_gap_SP[x]<dis:
			temp_out_SP2[-1].append(temp_out_SP[x+1])
		else:
			temp_out_SP2.append([temp_out_SP[x+1]])
	return temp_out_SP2

def Search_for_Common_BP(BP_info,dis_cff,num_cff):
	#dis_cff=50;num_cff=14
	out={}
	for k1 in BP_info.keys():
		all_SP={}
		all_LN={}
		global name_index
		name_index=sample_index(BP_info[k1].keys())
		global name_index_reverse
		name_index_reverse=hash_reverse(name_index)
		for k2 in BP_info[k1].keys():
			name_number=name_index[k2]
			for k3 in BP_info[k1][k2].keys():
				if k3.split('.')[-1]=='SPs':
					for k4 in BP_info[k1][k2][k3].keys():
						if not k4 in all_SP.keys():
							all_SP[k4]=[]
						all_SP[k4].append([BP_info[k1][k2][k3][k4],name_number])
				elif k3.split('.')[-1]=='LNs':
					for k4 in BP_info[k1][k2][k3].keys():
						if not k4 in all_LN.keys():
							all_LN[k4]=[]
						all_LN[k4].append(name_number)
						for k5 in BP_info[k1][k2][k3][k4]:
							if not k5 in all_LN.keys():
								all_LN[k5]=[]
							all_LN[k5].append(name_number)
					BP_info[k1][k2][k3+'.re']={}
					for k4 in  BP_info[k1][k2][k3].keys():
						for k5 in BP_info[k1][k2][k3][k4]:
							if not k5 in BP_info[k1][k2][k3+'.re'].keys():
								BP_info[k1][k2][k3+'.re'][k5]=[]
							BP_info[k1][k2][k3+'.re'][k5].append(k4)
		clu_SP=QC_SP(cluster_num(all_SP,dis_cff))
		class_SP=classify_SP_based_on_samples(clu_SP,all_SP)
		#class_SP was organized in has fasion, with each sub list containing SPs from 1 / 2 / 3 samples. could be useful for further analysis
		clu_LN=QC_SP(cluster_num(all_LN,dis_cff))
		class_LN=classify_LN_based_on_samples(clu_LN,all_LN)
		out[k1]=Associate_SP_LN(class_SP,class_LN)
	return out

def Write_Common_BP(Common_BPs):
	path_mkdir(workdir+'Common_BPs')
	for x in Common_BPs.keys():
		fo=open(workdir+'Common_BPs/Common.BP.'+x+'.bps','w')
		for x2 in Common_BPs[x]:
			print >>fo, ' '.join([str(i) for i in [x]+x2])
			print >>fo, ' '
		fo.close()

def svelter_file_names_readin(workdir):
	out=[]
	for x in os.listdir(workdir):
		if x.split('.')[-1]=='svelter':
			if not x in out:
				out.append(x)
	return out

def svelter_file_readin(workdir):
	global file_names
	file_names=svelter_file_names_readin(workdir)
	global file_index
	file_index=range(len(file_names))
	out={}
	for x in file_names:
		fin=open(workdir+x)
		pin=fin.readline().strip().split()
		for line in fin:
			pin=line.strip().split()
			if not pin[0] in out.keys():
				out[pin[0]]={}
			if not int(pin[1]) in out[pin[0]].keys():
				out[pin[0]][int(pin[1])]={}
			if not int(pin[2]) in out[pin[0]][int(pin[1])].keys():
				out[pin[0]][int(pin[1])][int(pin[2])]=[]
			out[pin[0]][int(pin[1])][int(pin[2])].append(pin+[file_index[file_names.index(x)]])
		fin.close()
	return out

def svelter_file_sample_collection(info_unit):
	#eg of info_unit: info[k1][k2][k3]=[['chr8', '93231103', '93231375', 'chr8:93231103:93231375', 'a/a', '/a', '-6.62007746853', 0], ['chr8', '93231103', '93231375', 'chr8:93231103:93231375', 'a/a', '/a', '-4.90683436002', 1], ['chr8', '93231103', '93231375', 'chr8:93231103:93231375', 'a/a', '/a', '-7.17113907948', 2]]
	out_index=[]
	for x in info_unit:
		if not x[-1] in out_index:
			out_index.append(x[-1])
	out_index.sort()
	return out_index
def svelter_file_interprete(workdir):
	info=svelter_file_readin(workdir)
	info_all={}
	for k1 in info.keys():
		for k2 in info[k1].keys():
			for k3 in info[k1][k2].keys():

				for length in range(len(file_names)+1)[1:]:
					if len(info[k1][k2][k3])==length:
						minus_length=len(file_names)-length
						if not minus_length in info_all.keys():
							info_all[minus_length]={}
						if not k1 in info_all[minus_length].keys():
							info_all[minus_length][k1]={}
						if not k2 in info_all[minus_length][k1].keys():
							info_all[minus_length][k1][k2]={}
						if not k3 in info_all[minus_length][k1][k2].keys():
							info_all[minus_length][k1][k2][k3]=[]
						info_all[minus_length][k1][k2][k3]=info[k1][k2][k3]
				if len(info[k1][k2][k3])>len(file_names):
					print info[k1][k2][k3]

	for k1 in info.keys():
		info_all[k1]={}
		for k2 in info[k1].keys():
			for k3 in info[k1][k2].keys():
				if len(info[k1][k2][k3])==len(file_names):
					if not k2 in info_all[k1].keys():
						info_all[k1][k2]={}
					if not k3 in info_all[k1][k2].keys():
						info_all[k1][k2][k3]=[]
					info_all[k1][k2][k3]=info[k1][k2][k3]
				else:
					for i in range(len(file_names))[1:]:

	svelter_file_write_common(info_all,'Common.SVs.minus.0.vcf')

def svelter_file_write_common(info_all,fileout):
	fo=open(workdir+fileout,'w')
	print >>fo, '\t'.join(['chr', 'start', 'end', 'bp_info', 'ref', 'alt', 'score','sample_name'])
	for k1 in info_all.keys():
		for k2 in info_all[k1].keys():
			for k3 in info_all[k1][k2].keys():
				for k4 in info_all[k1][k2][k3]:
					k4+=[file_names[file_index.index(k4[-1])]]
					print >>fo, '\t'.join([str(i) for i in k4[:-2]+[k4[-1]]])
	fo.close()

def main():
	opts,args=getopt.getopt(sys.argv[1:],'o:h:S:',['workdir=','chromosome=','sample-file=','reference=','exclude='])
	dict_opts=dict(opts)
	global workdir 
	workdir=path_modify(dict_opts['--workdir'])


	global reference
	reference=dict_opts['--reference']
	chromo_names=chromosome_names_readin(reference)
	global num_cff
	num_cff=14
	global dis_cff
	dis_cff=50
	global chromo
	for chromo in chromo_names:
		if chromo==dict_opts['--chromosome']:
			global BP_info
			BP_info=read_info_in(workdir,chromo) #read in BP info for certain chromosome
			global Common_BPs
			Common_BPs=Search_for_Common_BP(BP_info,dis_cff,num_cff)
			Write_Common_BP(Common_BPs)

main()