#!/usr/bin/env python

#!Python
#Usage:
#BP.Merge.Multi.Sample.py --workdir workdir
#For debug use only
#command='SV.Merge.Multi.Sample.py vcf --workdir /scratch/remills_flux/xuefzhao/SV_discovery_index/download/SVelter.CommonBPs'
#sys.argv=command.split()
import os
import re
import getopt
import sys
import numpy
import scipy
import random
def print_read_me():
	print 'SV.Merge.Multi.Sample.py V1.0 02-03-2016'
	print 'Contact: Xuefang Zhao (xuefzhao@umich.edu)'
	print ''
	print 'Usage: SV.Merge.Multi.Sample.py [Options] [Parameters]'
	print 'Options: '
	print '	svelter'
	print '	vcf'
	print ' bed'
	print 'Parameters:'
	print '	--workdir: workding directroy where all input files were placed, and also all output files will be written to'

if len(sys.argv)<2: 
	print_read_me()
else:
	def path_mkdir(path):
			if not os.path.isdir(path):
					os.system(r'''mkdir %s'''%(path))
	def path_modify(path):
		if not path[-1]=='/':
			path+='/'
		return path
	if sys.argv[1]=='svelter':
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
						sample_unit_index=svelter_file_sample_collection(info[k1][k2][k3])
						for length in range(len(file_names)+1)[1:]:
							if len(sample_unit_index)==length:
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
			for k1 in info_all.keys():
				fileout='Common.SVs.samples.minus.'+str(k1)+'.svelter'
				svelter_file_write_common(info_all[k1],fileout)
		def svelter_file_write_common(input_info,fileout):
			fo=open(workdir+fileout,'w')
			print >>fo, '\t'.join(['chr', 'start', 'end', 'bp_info', 'ref', 'alt', 'score','sample_name'])
			for k1 in input_info.keys():
				for k2 in input_info[k1].keys():
					for k3 in input_info[k1][k2].keys():
						for k4 in input_info[k1][k2][k3]:
							k4+=[file_names[file_index.index(k4[-1])]]
							print >>fo, '\t'.join([str(i) for i in k4[:-2]+[k4[-1]]])
			fo.close()
		def main():
			opts,args=getopt.getopt(sys.argv[2:],'o:h:S:',['workdir=','chromosome=','sample-file=','reference=','exclude='])
			dict_opts=dict(opts)
			global workdir 
			workdir=path_modify(dict_opts['--workdir'])
			svelter_file_interprete(workdir)
	elif sys.argv[1]=='vcf':
		def vcf_file_names_readin(workdir):
			out=[]
			for x in os.listdir(workdir):
				if x.split('.')[-1]=='vcf':
					if not x in out:
						out.append(x)
			return out
		def vcf_file_readin(workdir):
			global vcf_header
			vcf_header=[]
			global file_names
			file_names=vcf_file_names_readin(workdir)
			global file_index
			file_index=range(len(file_names))
			out={}
			for x in file_names:
				fin=open(workdir+x)
				for line in fin:
					pin=line.strip().split()
					if pin[0][0]=='#': 
						if vcf_header==[]:
							vcf_header.append(pin)
						if not vcf_header[-1][0]=='#CHROM':
							vcf_header.append(pin)
						continue
					if not pin[0] in out.keys():
						out[pin[0]]={}
					if not int(pin[1]) in out[pin[0]].keys():
						out[pin[0]][int(pin[1])]={}
					if not int(pin[2].split(';')[-2]) in out[pin[0]][int(pin[1])].keys():
						out[pin[0]][int(pin[1])][int(pin[2].split(';')[-2])]=[]
					out[pin[0]][int(pin[1])][int(pin[2].split(';')[-2])].append(pin+[file_index[file_names.index(x)]])
				fin.close()
			return out
		def vcf_file_sample_collection(info_unit):
			#eg of info_unit: info[k1][k2][k3]=[['chr8', '93231103', '93231375', 'chr8:93231103:93231375', 'a/a', '/a', '-6.62007746853', 0], ['chr8', '93231103', '93231375', 'chr8:93231103:93231375', 'a/a', '/a', '-4.90683436002', 1], ['chr8', '93231103', '93231375', 'chr8:93231103:93231375', 'a/a', '/a', '-7.17113907948', 2]]
			out_index=[]
			for x in info_unit:
				if not x[-1] in out_index:
					out_index.append(x[-1])
			out_index.sort()
			return out_index
		def vcf_file_interprete(workdir):
			info=vcf_file_readin(workdir)
			info_all={}
			for k1 in info.keys():
				for k2 in info[k1].keys():
					for k3 in info[k1][k2].keys():
						sample_unit_index=vcf_file_sample_collection(info[k1][k2][k3])
						for length in range(len(file_names)+1)[1:]:
							if len(sample_unit_index)==length:
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
			for k1 in info_all.keys():
				fileout='Common.SVs.samples.minus.'+str(k1)+'.vcf'
				vcf_file_write_common(info_all[k1],fileout)
		def vcf_file_write_common(input_info,fileout):
			fo=open(workdir+fileout,'w')
			for x in vcf_header[:-1]:
				print >>fo, '\t'.join(x)
			print >>fo, '\t'.join(vcf_header[-1][:-1]+['genotype','sample_name'])
			for k1 in input_info.keys():
				for k2 in input_info[k1].keys():
					for k3 in input_info[k1][k2].keys():
						for k4 in input_info[k1][k2][k3]:
							k4+=[file_names[file_index.index(k4[-1])]]
							print >>fo, '\t'.join([str(i) for i in k4[:-2]+[k4[-1]]])
			fo.close()
		def main():
			opts,args=getopt.getopt(sys.argv[2:],'o:h:S:',['workdir=','chromosome=','sample-file=','reference=','exclude='])
			dict_opts=dict(opts)
			global workdir 
			workdir=path_modify(dict_opts['--workdir'])
			vcf_file_interprete(workdir)

main()