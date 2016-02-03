#!/usr/bin/env python

#!Python
#Usage:
#BP.Merge.Trio.py --workdir workdir --sample-file sample-file
#For debug use only
#command='BP.Merge.Trio.py --workdir /scratch/remills_flux/xuefzhao/SV_discovery_index/download/ --sample-file /scratch/remills_flux/xuefzhao/SV_discovery_index/download/Trio_Info/Yoruban_File'
#sys.argv=command.split()
import os
import re
import getopt
import sys
import numpy
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

def chromosome_names_readin():
	global chromosome_names
	chromosome_names=[]
	for k1 in os.listdir(workdir):
		if 'reference' in k1 and os.path.isdir(workdir+k1):
			for k2 in os.listdir(workdir+k1):
				if k2.split('.')[-1]=='fa':
					ref_in=workdir+k1+'/'+k2
					fin=open(ref_in+'.fai')
					for line in fin:
						pin=line.strip().split()
						chromosome_names.append(pin[0])
	if chromosome_names==[]:
		print 'no reference file available under workdir:'+workdir

def chromosome_num_check(pin):
	temp_num=0
	for x in pin:
		if x in chromosome_names:
			temp_num+=1
		elif is_number(x)==False:
			temp_num+=1
	return temp_num

def main():
	opts,args=getopt.getopt(sys.argv[1:],'o:h:S:',['workdir=','sample-file=','reference=','exclude='])
	dict_opts=dict(opts)
	global sample_file
	sample_file=dict_opts['--sample-file']
	sample_name_readin()
	global workdir
	workdir=path_modify(dict_opts['--workdir'])
	chromosome_names_readin()
	global SP_info_rec
	SP_info_rec={}
	global files_names
	files_names=SP_file_readin(workdir,sample_file)
	All_Info=BP_info_readin(workdir,sample_file)
	SP_info_write(All_Info[0])
	LN_info_write(All_Info[1])
	Apply_SVelter_BPIntegrate()
	bp_files_write()

main()

