#!/usr/bin/env python

#!Python
#Usage:
#Compare.VCFs.py --file input.txt --reference ref.fa
#For debug use only
#command='Compare.VCFs.py --file /scratch/remills_flux/xuefzhao/SV_discovery_index/download/Yoruban/Yoruban.Trio.vcf.txt --reference /scratch/remills_flux/xuefzhao/reference/hg19_platinum/hg19_platinum.fa'
#sys.argv=command.split()
import os
import re
import getopt
import sys
import numpy
def sample_vcf_readin(vcfin):
	fin=open(vcfin)
	out_hash={}
	sv_signal=''
	for line in fin:
		pin=line.strip().split()
		if not pin[0][0]=='#':
			if chromosomes_check(pin)==1:
				if not ';'.join(pin[2].split(';')[:-1]) == sv_signal:
					sv_signal=';'.join(pin[2].split(';')[:-1])
					chrom=pin[0]
					start=int(pin[2].split(';')[1])
					end=int(pin[2].split(';')[-2])
					sv_type=pin[2].split(';')[-1]
					bps=pin[2].split(';')[1:-1]
					if not chrom in out_hash.keys():
						out_hash[chrom]={}
					if not start in out_hash[chrom].keys():
						out_hash[chrom][start]={}
					if not end in out_hash[chrom][start].keys():
						out_hash[chrom][start][end]=[]
				if not '[' in pin[4] and not ']' in pin[4]:
					sv_name=pin[4].replace('<','').replace('>','').replace(':','_')
					start_SV=int(pin[1])
					end_SV=end_pos_define(pin)
					out_hash[chrom][start][end].append([sv_name,start_SV,end_SV])
				else:
					out_hash[chrom][start][end].append(pin[4])
	return out_hash

def end_pos_define(pin):
	y=''
	for x in pin[7].split(';'):
		if 'END' in x:
			y=x.split('=')[1]
	if not y=='':
		return int(y)
	else:
		return 'Error'

def txt_file_readin(filein):
	fin=open(filein)
	samples={}
	for line in fin:
		pin=line.strip().split()
		vcf_hash=sample_vcf_readin(pin[0])
		samples[pin[0]]=vcf_hash
	fin.close()
	return samples

def BP_correction(vcf_hashes):
	label_hash={}
	rec=0
	for k1 in vcf_hashes.keys():
		rec+=1
		label_hash[k1]=rec
	total_SV_regions={}
	for k1 in vcf_hashes.keys():
		label=label_hash[k1]
		for k2 in vcf_hashes[k1].keys():
			if not k2 in total_SV_regions.keys():
				total_SV_regions[k2]={}
			for k3 in vcf_hashes[k1][k2].keys():
				if not k3 in total_SV_regions[k2].keys():
					total_SV_regions[k2][k3]={}
				for k4 in vcf_hashes[k1][k2][k3].keys():
					if not k4 in total_SV_regions[k2][k3].keys():
						total_SV_regions[k2][k3][k4]=[]
					total_SV_regions[k2][k3][k4].append([k3,k4,label])
	total_sv_list={}
	for k1 in total_SV_regions.keys():
		total_SV_regions[k1]=[]
		for k2 in sorted(total_SV_regions[k1].keys()):
			for k3 in sorted(total_SV_regions[k1][k2].keys()):
				total_SV_regions[k1]+=total_SV_regions[k1][k2][k3]

def chromosomes_readin(reference):
	ref_index=reference+'.fai'
	if not os.path.isfile(ref_index):
		print 'Error: reference not indexed!'
	else:
		out=[]
		fin=open(ref_index)
		for line in fin:
			pin=line.strip().split()
			out.append(pin[0])
		fin.close()
		return out

def chromosomes_check(pin):
	num_chr=0
	for x in pin[2].split(';'):
		if x in chromos:
			num_chr+=1
	return num_chr

def main():
	opts,args=getopt.getopt(sys.argv[1:],'o:h:S:',['file=','reference='])
	dict_opts=dict(opts)
	filein=dict_opts['--file']
	reference=dict_opts['--reference']
	global chromos
	chromos=chromosomes_readin(reference)
	vcf_hashes=txt_file_readin(filein)


