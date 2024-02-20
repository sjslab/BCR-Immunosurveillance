# MRDARCY developed by Rachael Bashford-Rogers (2018)
# at the University of Oxford and Univeristy of Cambridge
# E-mail: rbr1@well.ox.ac.uk

# If you find the methods in MRDARCY, please cite the following reference:
# Bashford-Rogers, R. J., et al. (2016). "Eye on the B-ALL: B-cell receptor 
# repertoires reveal persistence of numerous B-lymphoblastic leukemia 
# subclones from diagnosis to relapse." Leukemia 30(12): 2312-2321.

# Copyright (C) 2018  Dr Rachael Bashford-Rogers

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.


#!/usr/bin/python
import math
import sys
from collections import defaultdict
import os
#import commands
import sys
from operator import itemgetter
import networkx as nx
import numpy as np
#from itertools import izip
import copy
import os.path

def fasta_iterator(fh):
  while True:
    line = fh.readline()
    if line.startswith('>'): break  
  while True:
    header = line[1:-1].rstrip()
    sequence = fh.readline().rstrip()
    while True:
      line = fh.readline()
      if not line: break
      if line.startswith('>'): break
      sequence += line.rstrip()
    yield(header, sequence)
    if not line: return

def Cluster_i(seqfile, outfile, edge_lengths,cd_hit_directory):
  command= cd_hit_directory+"cd-hit -i "+seqfile+" -o "+outfile+" -c "+str(edge_lengths)+" -d 200 -T 10  -M 0 -AL 40 -n 5 -AS 30"
  os.system(command)
  return()

def Write_output(out, file):
  fh=open(file, "a")
  fh.write(out)
  fh.close()
  return ()

def CDR3_define_sequences(grouped, pat, dir, combined_sequences_file,header_file,CDR3_info_file,merge_cluster_file):
  for f in [CDR3_info_file,combined_sequences_file,merge_cluster_file]:
    fh=open(f,"w")
    fh.close()
  seqs = Functions.Tree()
  ids,ind,chains_use = [],0,{}
  fh=open(CDR3_info_file,"w")
  fh.close()
  out, ind = '#sample\tid\tv\tj\tCDR3 region\tV:J:CDR3\n',0
  cluster_codes,codes_clusters=Tree(),Tree()
  all_clusters,total_seq = Tree(),0
  all_seqs={}
  for info in grouped[pat]:
    sample = info[1]
    print (sample)
    seq_file = info[4]+"ORIENTATED_SEQUENCES/NETWORKS/Fully_reduced_"+sample+".fasta"
    annot_file = info[4]+"ORIENTATED_SEQUENCES/ANNOTATIONS/TMP/Annotation_"+sample+".txt"
    cluster_file = info[4]+"ORIENTATED_SEQUENCES/NETWORKS/Cluster_identities_"+sample+".txt"
    fh=open(seq_file,"r")
    for header,sequence in fasta_iterator(fh):
      all_seqs[header.split("__")[0]] = [header+"|"+sample, sequence]
    fh.close()
    fh=open(annot_file,"r")
    out, ind = '#sample\tid\tV\tJ\tCDR3\tcode\n',0
    codes, inv_codes={},Tree()
    for l in fh:
      if(l[0]!="#"):
        l=l.strip().split()
        if(len(l)>=18):
          v,j,cdr3 = l[1],l[13],l[16]
          #v,j = v.split("*")[0], j.split("*")[0]
          if(len(cdr3)<=7):cdr3 = "-"
          out=out+sample+"\t"+l[0]+"\t"+v+"\t"+j+"\t"+cdr3+"\t"+v+":"+j+":"+cdr3+"\n"
          if(cdr3 != "-"):
            codes[l[0].split("__")[0]] = v+":"+j+":"+cdr3
            inv_codes[v+":"+j+":"+cdr3][l[0].split("__")[0]].value = 1
          ind = ind+1
          if(ind>500):
            Write_output(out, CDR3_info_file)
            out, ind = '',0


def CDR3_defined_sequences_non_IMGT(pat, output_dir, seq_files, annot_files,samples,cd_hit_directory,tmp_file,edge_lengths,merge_cluster_file_count):
  all_seqs = Tree()
  outseq,indseq = '',0
  step = 1
  if(step==1):
    for f in [CDR3_info_file,combined_sequences_file]:
      fh=open(f,"w")
      fh.close()
    for i in range(0,len(samples)):
      sample = samples[i]
      fh=open(seq_files[i],"r")
      for header,sequence in fasta_iterator(fh):
        all_seqs[header.split("__")[0]] = [header+"|"+sample, sequence]
      fh.close()
      fh=open(annot_files[i],"r")
      out, ind = '#sample\tid\tV\tJ\tCDR3\tcode\n',0
      codes, inv_codes={},Tree()
      for l in fh:
        if(l[0]!="#"):
          l=l.strip().split()
          if(len(l)<18):
            print (len(l))
          if(len(l)>=18):
            v,j,cdr3 = l[1],l[13],l[16]
            v,j = v.split("*")[0], j.split("*")[0]
            if(len(cdr3)>=0*7):
              id = l[0].split("__")[0]
              cdr3 = all_seqs[id][1]
              out=out+sample+"\t"+id+"\t"+v+"\t"+j+"\t"+cdr3+"\t"+v+":"+j+":"+cdr3+"\n"
              outseq = outseq+">"+all_seqs[id][0]+"\n"+cdr3.upper()+"\n"
              codes[l[1].split("__")[0]] = v+":"+j+":"+cdr3
              inv_codes[v+":"+j+":"+cdr3][id].value = 1
              ind = ind+1
              if(ind>500):
                Write_output(out, CDR3_info_file)
                Write_output(outseq,combined_sequences_file)
                out, ind,outseq = '',0,''
      Write_output(out, CDR3_info_file)
      Write_output(outseq,combined_sequences_file)
      out, ind,outseq = '',0,''
      fh.close()
  ### cluster
  step = 21
  if(step==2):
    for f in [merge_cluster_file, merge_cluster_file_count]:
      fh=open(f,"w")
      fh.close()
    Cluster_i(combined_sequences_file, tmp_file+"clust",edge_lengths,cd_hit_directory)
    Get_clusters_merged (tmp_file, combined_sequences_file,merge_cluster_file,CDR3_info_file,merge_cluster_file_count)
  return()

def CDR3_defined_sequences_SC(pat, output_dir, seq_files, annot_files,samples,cd_hit_directory,tmp_file1,tmp_file2,edge_lengths,merge_cluster_file_count,combined_sequences_file1,combined_sequences_file2):
  Step1_combine_sequences(pat, output_dir, seq_files, annot_files,samples,cd_hit_directory,combined_sequences_file1,combined_sequences_file2)
  Step2_cluster_combined_sequences(pat, output_dir, seq_files, annot_files,samples,cd_hit_directory,tmp_file1,tmp_file2,edge_lengths,merge_cluster_file_count,combined_sequences_file1,combined_sequences_file2)
  return()

def Step1_combine_sequences(pat, output_dir, seq_files, annot_files,samples,cd_hit_directory,combined_sequences_file1,combined_sequences_file2):
  all_seqs = Tree()
  outseq,indseq = '',0
  step = 1
  if(step==1):
    if(os.path.exists(output_dir)==False):os.system("mkdir "+output_dir)
    for f in [CDR3_info_file,combined_sequences_file1, combined_sequences_file2]:
      fh=open(f,"w")
      fh.close()
    for i in range(0,len(samples)):
      #print samples[i]
      sample = samples[i]
      second_fasta = seq_files[i].replace("IGH","IGL")
      #print second_fasta
      fh=open(seq_files[i],"r")
      all_seqs,all_seqs2 = {},{}
      for header,sequence in fasta_iterator(fh):
        header = header.split("_contig")[0]
        all_seqs[header] = [header+"||"+sample, sequence]
      fh.close()
      fh=open(second_fasta,"r")
      for header,sequence in fasta_iterator(fh):
        header = header.split("_contig")[0]
        all_seqs2[header] = [header+"||"+sample, sequence]
      fh.close()
      fh=open(annot_files[i],"r")
      codes, inv_codes={},Tree()
      out, ind, outseq= '#sample\tid\tIGHV\tIGHJ\tIGH_CDR3\tcode\tisotype\tIGLV\tIGLJ\tl/k\treads_IGH\treads_IGL\tmm_IGHV\tmm_IGHJ\tmm_IGLV\tmm_IGLJ\n',0,''
      outseq2 = ''
      for l in fh:
        if(l[0]!="#"):
          l=l.strip().split()
          id,v,j,cdr3,iso = l[0],l[5],l[6],l[8],l[3]
          if(l[1]!='-'):
            if(id in all_seqs):
              if(id in all_seqs2):
                outseq2 = outseq2+">"+all_seqs2[id][0]+"\n"+all_seqs2[id][1]+"\n"
                vl,jl,isol,reads_IGH,reads_IGL = l[17],l[18],l[14],l[4],l[16]
                mm_v,mmj,mm_vl, mm_jl = l[11],l[12],l[23],l[24]
                seq = all_seqs[id][1]
                outseq = outseq+">"+all_seqs[id][0]+"\n"+seq+"\n"
                cod = v[0:5]+":"+j+":"+cdr3
                codes[id] = cod
                out=out+"\t".join(map(str,[sample, id, v,j,cdr3,cod, iso, vl,jl,isol,reads_IGH,reads_IGL,mm_v,mmj,mm_vl, mm_jl]))+"\n"
                inv_codes[v+":"+j+":"+cdr3][id].value = 1
                ind = ind+1
                if(ind>500):
                  Write_output(out, CDR3_info_file)
                  Write_output(outseq,combined_sequences_file1)
                  Write_output(outseq2,combined_sequences_file2)
                  out, ind,outseq,outseq2 = '',0,'',''
              #else:print "no IGL",id, sample
      fh.close()
      Write_output(out, CDR3_info_file)
      Write_output(outseq,combined_sequences_file1)
      Write_output(outseq2,combined_sequences_file2)
      out, ind,outseq,outseq2 = '',0,'',''
      fh.close()
  return()

def Step2_cluster_combined_sequences(pat, output_dir, seq_files, annot_files,samples,cd_hit_directory,tmp_file1,tmp_file2,edge_lengths,merge_cluster_file_count,combined_sequences_file1,combined_sequences_file2):
  ### cluster
  step = 2
  if(step==2):
    for f in [merge_cluster_file, merge_cluster_file_count]:
      fh=open(f,"w")
      fh.close()
    python  = "/usr/bin/python"
    command1 = python+" /users/immune-rep/mfj169/GENERAL_PROGRAMS/NETWORK_GENERATION/Network_generation_Single_cell_2.0.py "+output_dir+" "+pat+"_A "+combined_sequences_file1+" "+output_dir+"CDR3_information_"+pat+".txt"
    command2 = python+" /users/immune-rep/mfj169/GENERAL_PROGRAMS/NETWORK_GENERATION/Network_generation_Single_cell_2.0.py "+output_dir+" "+pat+"_B "+combined_sequences_file2+" "+output_dir+"CDR3_information_"+pat+".txt"
    os.system(command1)
    os.system(command2)
    ####Cluster_i(combined_sequences_file1, tmp_file1+"clust",edge_lengths,cd_hit_directory)
    ####Cluster_i(combined_sequences_file2, tmp_file2+"clust",edge_lengths,cd_hit_directory)
    cluster_file_1 = output_dir+"Cluster_identities_"+pat+"_A.txt"
    cluster_file_2 = output_dir+"Cluster_identities_"+pat+"_B.txt"
    Get_clusters_merged_SC (cluster_file_1, cluster_file_2, merge_cluster_file, merge_cluster_file_count)

    #Get_clusters_merged_SC (tmp_file1, tmp_file2, combined_sequences_file1, combined_sequences_file2,merge_cluster_file,CDR3_info_file,merge_cluster_file_count)
    #os.system("rm "+tmp_file1+"clust "+tmp_file1+"clust.bak.clstr "+tmp_file1+"clust.clstr")
    #os.system("rm "+tmp_file2+"clust "+tmp_file2+"clust.bak.clstr "+tmp_file2+"clust.clstr")
  return()

def Get_clusters (file):
  fh=open(file,"r")
  clusters,inv_cluster=Tree(),{}
  clust = '-1'
  for l in fh:
    l=l.strip().split()
    if(l[0].count(">")!=0):
        clust = l[1]
    else:
        id1 = l[2].replace("...","").replace(">","")
        id = id1.split("__")[0].split("|")[0]+"|"+id1.split("|")[len(id1.split("|"))-1]
        clusters["C"+clust][id].value=1
        inv_cluster[id]=["C"+clust,id1]
  fh.close()
  return(clusters,inv_cluster)

def CDR3_defined_sequences_IMGT(pat, output_dir, seq_files, annot_files,samples,cd_hit_directory,tmp_file,edge_lengths,merge_cluster_file_count,combined_sequences_file,summary_merge_clustering_count):
  step = 1
  if(step==1):
    for f in [CDR3_info_file,combined_sequences_file]:
      fh=open(f,"w")
      fh.close()
    for i in range(0,len(samples)):
      all_seqs = Tree()
      outseq,indseq = '',0
      sample = samples[i]
      CDR3_IMGT_merged_file = output_dir+"Cluster_CDR3_defined_identities_"+sample+".txt"
      fh=open(seq_files[i],"r")
      lens = []
      for header,sequence in fasta_iterator(fh):
        id = header.split("__")[0]
        all_seqs[id] = [header+"|"+sample, sequence]
        lens.append(len(sequence))
      fh.close()
      print (seq_files[i], sum(lens)*1.0/len(lens))
      fh=open(annot_files[i],"r")
      out, ind = '#sample\tid\tV\tJ\tCDR3\tcode\n',0
      codes, inv_codes={},Tree()
      passed, failed = 0,0
      for l in fh:
        if(l[0]!="#"):
          l=l.strip().split("\t")
          id = l[1].split("__")[0]
          if(l[2].count("unknown")==0 and id in all_seqs and l[2].count("No rearrangement found")==0):
            if(len(l)>=16):
              cdr3,v,j = l[15],l[3].split()[1],l[4].split()[1]
              v,j = v.split("*")[0], j.split("*")[0]
              passed = passed+1
              if(len(cdr3)>=7):
                out=out+sample+"\t"+id+"\t"+v+"\t"+j+"\t"+cdr3+"\t"+v+":"+j+":"+cdr3+"\n"
                seq1 = all_seqs[id][1]
                #seq1 = l[6]
                outseq = outseq+">"+all_seqs[id][0]+"\n"+seq1.upper()+"\n"
                codes[id] = v+":"+j+":"+cdr3
                inv_codes[v+":"+j+":"+cdr3][id].value = 1
                ind = ind+1
                if(ind>500):
                  Write_output(out, CDR3_info_file)
                  Write_output(outseq,combined_sequences_file)
                  out, ind,outseq = '',0,''
            else:failed = failed+1
          else:
            if(l[2].count("unknown")!=0):print (3) 
            elif(id not in all_seqs):a = 1
             # print id
             # print 4, id,seq_files[i]
      Write_output(out, CDR3_info_file)
      Write_output(outseq,combined_sequences_file)
      out,outseq = '',''
      missing = {}
      for id in all_seqs:
        if(id not in codes):
          missing[id] = 1
      Write_output(out, CDR3_info_file)
      Write_output(outseq,combined_sequences_file)
      print ("passed",passed, "all",len(all_seqs))
      Write_output(outseq,combined_sequences_file)
      out, ind,outseq = '',0,''
      fh.close()
      missing = {}
  ### cluster
  step = 2
  if(step==2):
    for f in [merge_cluster_file, merge_cluster_file_count]:
      fh=open(f,"w")
      fh.close()
    Cluster_i(combined_sequences_file, tmp_file+"clust",edge_lengths,cd_hit_directory)
    Get_clusters_merged (tmp_file, combined_sequences_file,merge_cluster_file,CDR3_info_file,merge_cluster_file_count)
    #os.system("rm "+tmp_file+"clust "+tmp_file+"clust.bak.clstr "+tmp_file+"clust.clstr")
  step = 3
  if(step==3):
    Summarise_clusters(combined_sequences_file,merge_cluster_file,summary_merge_clustering_count)
  return()

def Summarise_clusters(combined_sequences_file,merge_cluster_file,summary_merge_clustering_count):
  fh=open(summary_merge_clustering_count,"w")
  fh.close()
  fh=open(merge_cluster_file,"r")
  cluster = {}
  samples = {}
  for l in fh: 
    if(l[0]!="#"):
      l=l.strip().split()
      id,clust = l[1].split("|")[0],l[0]
      cluster[id] = clust
      samples[l[2]] = 1
  fh.close()
  print (len(cluster))
  fh=open(combined_sequences_file,"r")
  samples1 = []
  for s in samples: 
    samples1.append(s)
  samples1.sort()
  l = len(samples1)
  sample_ind, ind = {},0
  for s in samples1:
    sample_ind[s] = ind
    ind = ind+1
  cluster_size = {}
  for header,sequence in fasta_iterator(fh):
    id = header.split("__")[0]
    freq = sum(map(int, header.split("__")[1].split("|")[0].split("_")))
    sample = header.split("|")
    sample = sample[len(sample)-1]
    if(id in cluster):
      clust = cluster[id]
      if(clust in cluster_size): 
        x = copy.deepcopy(cluster_size[clust])
        x[sample_ind[sample]] = x[sample_ind[sample]]+freq
        cluster_size[clust] = x
      else:
        x = [0]*l
        x[sample_ind[sample]] = x[sample_ind[sample]]+freq
        cluster_size[clust] = x
  fh.close()
  out = '#clone\t'+"\t".join(samples1)+"\n"
  ind = 0
  for c in cluster_size:
    out=out+"\t".join(map(str, [c]+cluster_size[c]))+"\n"
    ind= ind+1
    if(ind>500):
      Write_output(out, summary_merge_clustering_count)
      out, ind = '',0
  Write_output(out, summary_merge_clustering_count)
  return()

def Get_joint_clusters_SC(file1, file2):
  fh=open(file1,"r")
  clusters,inv_cluster1,inv_cluster=Tree(),{},{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id,clust = l[2],"C"+l[1]
      inv_cluster1[id]=[clust,id]
  fh.close()
  fh=open(file2,"r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id,clust = l[2],"C"+l[1]
      if(id in inv_cluster1):
        c1 = inv_cluster1[id][0]+"_"+clust
        clusters[c1][id].value = 1
        inv_cluster[id]=[c1,id]
  del inv_cluster1
  return(clusters,inv_cluster)

def Get_clusters_SC (file1, file2):
  fh=open(file1,"r")
  clusters,inv_cluster1,inv_cluster=Tree(),{},{}
  for l in fh:
    l=l.strip().split()
    id1 = l[2].replace("...","").replace(">","")
    id = id1.split("__")[0].split("|")[0]+"|"+id1.split("|")[2]
    inv_cluster1[id]=["C"+l[0],id1]
  fh.close()
  fh=open(file2,"r")
  for l in fh:
    l=l.strip().split()
    id1 = l[2].replace("...","").replace(">","")
    id = id1.split("__")[0].split("|")[0]+"|"+id1.split("|")[2]
    if(id in inv_cluster1):
      c = "C"+l[0]
      c1 = inv_cluster1[id][0]+"|"+c
      clusters[c1][id].value = 1
      inv_cluster[id]=[c1,id1]
  del inv_cluster1
  return(clusters,inv_cluster)

def Compare_CDR3s(cdr3_sub, end_indent, n_mutations_allowed,G):
  cdr3_sub.sort(key=lambda x: x[1],reverse=True)
  for i1 in range(0,len(cdr3_sub)):
    cdr3_1 = cdr3_sub[i1][0]
    words1 = cdr3_1[end_indent:len(cdr3_1)-end_indent]
    l1 = int(math.floor(len(words1)/n_mutations_allowed))
    words1 = [words1[((i)*l1):((i+1)*l1)] for i in range(0,n_mutations_allowed)]
    words1 = [words1[i] for i in range(0,len(words1)) if len(words1[i])>4]
    index1 = [cdr3_1.index(words1[i]) for i in range(0,len(words1))]
    for j in range(i1+1,len(cdr3_sub)):
      cdr3_2 = cdr3_sub[j][0]
      matches = [i for i in range(0,len(words1)) if cdr3_2.count(words1[i])!=0]
      if(len(matches)>0):
        for i in matches:
          index2i = cdr3_2.index(words1[i])
          index1i = index1[i]
          if(index2i>index1i):
            s1 = cdr3_1
            s2 = cdr3_2[index2i-index1i:cdr3_sub[j][1]]
          else:
            s1 = cdr3_1[index1i-index2i:len(cdr3_1)]
            s2 = cdr3_2
        mm= [i for i in range(0,min([len(s1),len(s2)])) if s1[i] != s2[i]]
        if(len(mm)<=n_mutations_allowed):
            G.add_edge(cdr3_1, cdr3_2)
  return(G)

def Get_clusters_merged_SC (cluster_file_1, cluster_file_2, merge_cluster_file, merge_cluster_file_count):
  clusters,inv_cluster = Get_joint_clusters_SC(cluster_file_1, cluster_file_2)
  out,ind = "#cluster_ID\tSequence_ID\tSample_ID\n",0
  files = Tree()
  total_found = 0
  for c in clusters: 
    for id in clusters[c]:
      file_id = id.split("||")[1]
      out=out+"\t".join([c,id,file_id])+"\n"
      files[c][file_id].value = 1
      ind = ind+1
      total_found = total_found+1
      if(ind>500):
        Write_output(out, merge_cluster_file)
        out, ind = '',0
  Write_output(out, merge_cluster_file)
  out, ind = '',0
  out, ind = '#cluster_ID\tnumber_of_sequences\tnumber_of_samples\n',0
  for c in clusters:
    out=out+"\t".join(map(str,[c,len(clusters[c]), len(files[c])]))+"\n"
    ind =ind+1
    if(ind>500):
      Write_output(out, merge_cluster_file_count)
      out, ind = '',0
  Write_output(out, merge_cluster_file_count)
  out, ind = '',0
  return()

def OLD():
  #n_mutations_allowed = 6
  #length_var = 10
  #end_indent = 3
  #clusters,inv_cluster = Get_clusters_SC(tmp_file1+"clust.bak.clstr", tmp_file2+"clust.bak.clstr")
  fh=open(CDR3_info_file,"r")
  cdr3s,cdr3_inv,potential_clusters = Tree(),{},Tree()
  clusters_merged = Tree()
  unique_cluster_ind,total_IDs = 0,0
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id = l[1]+"||"+l[0]
      total_IDs = total_IDs+1
      CDR3 = l[4].upper()
      if(id in inv_cluster):
        CDR3_length, CDR3_features = len(CDR3), l[2]+"|"+l[3]+"|"+inv_cluster[id][0]
        cdr3s[CDR3_features][CDR3_length][CDR3][id].value = 1
        potential_clusters[CDR3_features][id].value = 1
      else:
        unique_cluster_ind = unique_cluster_ind+1
        clusters_merged["C"+str(unique_cluster_ind)][id].value = 1
  fh.close()
  done = {}
  for CDR3_features in cdr3s:
    if(len(potential_clusters[CDR3_features])==1):
      unique_cluster_ind = unique_cluster_ind+1
      for id in potential_clusters[CDR3_features]:
        if(id not in done):
          clusters_merged["C"+str(unique_cluster_ind)][id].value = 1
          done[id] = 1
        else:print ("Error 1",id)
    else:#if(len(potential_clusters[CDR3_features])>1):
      lengths = []
      G=nx.Graph()
      cdr3_IDs_sub,len_cdr3s= Tree(),0
      for l in cdr3s[CDR3_features]:
        lengths.append(l)
        len_cdr3s = len_cdr3s+len(cdr3s[CDR3_features][l])
      lengths.sort()
      cdr3_count_sub = {}
      for l in lengths:
        cdr3_sub = []
        for i in range(l-length_var, l+length_var):
          for cdr3 in cdr3s[CDR3_features][i]:
            cdr3_sub.append([cdr3,i])
            for id in cdr3s[CDR3_features][i][cdr3]:
              cdr3_IDs_sub[cdr3][id].value = 1
              cdr3_count_sub[cdr3] = 1
        G = Compare_CDR3s(cdr3_sub, end_indent, n_mutations_allowed,G)
      con= nx.connected_components(G)
      cdr3_done = {}
      for c in con:
        unique_cluster_ind = unique_cluster_ind+1
        for cdr3 in c:
          cdr3_done[cdr3] = 1
          for id in cdr3_IDs_sub[cdr3]:
            if(id not in done):
              clusters_merged["C"+str(unique_cluster_ind)][id].value = 1
              done[id] = 1
            #else:print "Error 2",id
        else:
          unique_cluster_ind = unique_cluster_ind+1
          for cdr3 in cdr3_sub:
            cdr3_done[cdr3[0]] = 1
            for id in cdr3_IDs_sub[cdr3[0]]:
              if(id not in done):
                clusters_merged["C"+str(unique_cluster_ind)][id].value = 1
                done[id] = 1
              #else:print "Error 3",id
      for cdr3 in cdr3_count_sub:
        if(cdr3 not in cdr3_done):
          unique_cluster_ind = unique_cluster_ind+1
          for id in cdr3_IDs_sub[cdr3]:
            if(id not in done):
              clusters_merged["C"+str(unique_cluster_ind)][id].value = 1
              done[id] = 1
            else:print ("Error 4",id)
  out,ind = "#cluster_ID\tSequence_ID\tSample_ID\n",0
  files = Tree()
  total_found = 0
  for c in clusters_merged:
    for id in clusters_merged[c]:
      total_found = total_found+1
      if(id in inv_cluster):
        file_id=inv_cluster[id][1].split("|")
        file_id = file_id[len(file_id)-1]
      else:file_id = id.split("|")[1]
      files[c][file_id].value = 1
      out=out+"\t".join([c,id,file_id])+"\n"
      ind = ind+1
      if(ind>500):
        Write_output(out, merge_cluster_file)
        out, ind = '',0
  Write_output(out, merge_cluster_file)
  out, ind = '#cluster_ID\tnumber_of_sequences\tnumber_of_samples\n',0
  for c in clusters_merged:
    out=out+"\t".join(map(str,[c,len(clusters_merged[c]), len(files[c])]))+"\n"
    ind =ind+1
    if(ind>500):
      Write_output(out, merge_cluster_file_count)
      out, ind = '',0
  Write_output(out, merge_cluster_file_count)
  out, ind = '',0
  print ("total ids:",total_IDs,"found:", total_found, len(done))
  return()

def Get_clusters_merged (tmp_file, combined_sequences_file,merge_cluster_file,CDR3_info_file,merge_cluster_file_count):
  fh=open(CDR3_info_file,"r")
  code = Tree()
  ids = {}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id,cdr3 = l[1]+"|"+l[0],l[5]
      code[cdr3][id].value = 1
      ids[id] =1
  fh.close()
  print (len(ids))
  n_mutations_allowed = 10
  length_var = 5
  end_indent = 3
  clusters,inv_cluster = Get_clusters(tmp_file+"clust.clstr")
  fh=open(CDR3_info_file,"r")
  cdr3s,cdr3_inv,potential_clusters = Tree(),{},Tree()
  clusters_merged = Tree()
  unique_cluster_ind,total_IDs = 0,0
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id = l[1]+"|"+l[0]
      total_IDs = total_IDs+1
      CDR3 = l[4].upper()
      if(id in inv_cluster):
        CDR3_length, CDR3_features = len(CDR3), l[2]+"|"+l[3]+"|"+inv_cluster[id][0]
        cdr3s[CDR3_features][CDR3_length][CDR3][id].value = 1
        potential_clusters[CDR3_features][id].value = 1
      else:
        unique_cluster_ind = unique_cluster_ind+1
        clusters_merged["C"+str(unique_cluster_ind)][id].value = 1
        print ("Unclustered",id, l)
  fh.close()
  done = {}
  for CDR3_features in cdr3s:
    if(len(potential_clusters[CDR3_features])==1):
      unique_cluster_ind = unique_cluster_ind+1
      for id in potential_clusters[CDR3_features]:
        if(id not in done):
          clusters_merged["C"+str(unique_cluster_ind)][id].value = 1
          done[id] = 1
        else:print ("Error 1",id)
    else:#if(len(potential_clusters[CDR3_features])>1):
      lengths = []
      G=nx.Graph()
      cdr3_IDs_sub,len_cdr3s= Tree(),0
      for l in cdr3s[CDR3_features]:
        lengths.append(l)
        len_cdr3s = len_cdr3s+len(cdr3s[CDR3_features][l])
      lengths.sort()
      cdr3_count_sub = {}
      for l in lengths:
        cdr3_sub = []
        for i in range(l-length_var, l+length_var):
          for cdr3 in cdr3s[CDR3_features][i]:
            cdr3_sub.append([cdr3,i])
            for id in cdr3s[CDR3_features][i][cdr3]:
              cdr3_IDs_sub[cdr3][id].value = 1
              cdr3_count_sub[cdr3] = 1
        G = Compare_CDR3s(cdr3_sub, end_indent, n_mutations_allowed,G)
      con= nx.connected_components(G)
      cdr3_done = {}
      for c in con:
        unique_cluster_ind = unique_cluster_ind+1
        for cdr3 in c:
          cdr3_done[cdr3] = 1
          for id in cdr3_IDs_sub[cdr3]:
            if(id not in done):
              clusters_merged["C"+str(unique_cluster_ind)][id].value = 1
              done[id] = 1
            #else:print "Error 2",id
        else:
          unique_cluster_ind = unique_cluster_ind+1
          for cdr3 in cdr3_sub:
            cdr3_done[cdr3[0]] = 1
            for id in cdr3_IDs_sub[cdr3[0]]:
              if(id not in done):
                clusters_merged["C"+str(unique_cluster_ind)][id].value = 1
                done[id] = 1
              #else:print "Error 3",id
      for cdr3 in cdr3_count_sub:
        if(cdr3 not in cdr3_done):
          unique_cluster_ind = unique_cluster_ind+1
          for id in cdr3_IDs_sub[cdr3]:
            if(id not in done):
              clusters_merged["C"+str(unique_cluster_ind)][id].value = 1
              done[id] = 1
            else:print ("Error 4",id)
  ##overall clustering with identical CDR3s included
  G=nx.Graph()
  for c in clusters_merged: 
    prev = ''
    for id in clusters_merged[c]: 
      if(prev!=''):G.add_edge(prev,id)
      prev = id
  for cdr3 in code: 
    prev = ''
    for id in code[cdr3]:
      if(prev!=''):G.add_edge(prev,id)
      prev = id
  con= nx.connected_components(G)
  clusters_final = Tree()
  clust_ind = 0
  done = {}
  for c in con:
    clust_ind =clust_ind+1
    for id in list(c): 
        clusters_final["A"+str(clust_ind)][id].value = 1
        done[id] = 1
  for id in ids: 
    if(id not in done):
      clust_ind =clust_ind+1
      clusters_final["A"+str(clust_ind)][id].value = 1
      done[id] = 1
  print ("done", len(done))
  clusters_merged=clusters_final
  print (len(clusters_merged),len(clusters_final))
  out,ind = "#cluster_ID\tSequence_ID\tSample_ID\n",0
  files = Tree()
  total_found = 0
  for c in clusters_merged:
    for id in clusters_merged[c]:
      total_found = total_found+1
      if(id in inv_cluster):
        file_id=inv_cluster[id][1].split("|")
        file_id = file_id[len(file_id)-1]
      else:file_id = id.split("|")[1]
      files[c][file_id].value = 1
      out=out+"\t".join([c,id,file_id])+"\n"
      ind = ind+1
      if(ind>500):
        Write_output(out, merge_cluster_file)
        out, ind = '',0
  Write_output(out, merge_cluster_file)
  out, ind = '#cluster_ID\tnumber_of_sequences\tnumber_of_samples\n',0
  for c in clusters_merged:
    out=out+"\t".join(map(str,[c,len(clusters_merged[c]), len(files[c])]))+"\n"
    ind =ind+1
    if(ind>500):
      Write_output(out, merge_cluster_file_count)
      out, ind = '',0
  Write_output(out, merge_cluster_file_count)
  out, ind = '',0
  print ("total ids:",total_IDs,"found:", total_found, len(done))
  return()

def Initialise_clonal_analysis(pat, output_dir, merge_cluster_file, merge_cluster_file_count,clone_alignment_files,fasttree):
  clone_sizes = []
  fh = open(merge_cluster_file_count,"r")
  for l in fh:
    if(l[0]!="#"):
      l = l.strip().split()
      clone, size, n_samples = l[0],int(l[1]),int(l[2])
      if(size>0 or n_samples>1):
        clone_sizes.append((clone, size, n_samples))
        #if(n_samples>=2): print l
  fh.close()
  clone_sizes.sort(key=lambda x: x[1], reverse = True)
  clones_plot= {}
  add_large_clones ="F"
  out = "#pat\tclone\tsize\tn_samples\n"
  if(add_large_clones =="True"):
    n_plot = min([12,len(clone_sizes)])
    for i in range(0,n_plot):
      clones_plot[clone_sizes[i][0]] = str(clone_sizes[i][0])+"_"+str(clone_sizes[i][1])+"_"+str(clone_sizes[i][2])
      out = out+pat+"\t"+clones_plot[clone_sizes[i][0]]+"\t"+str(clone_sizes[i][1])+"\t"+str(clone_sizes[i][2])+"\n"
    clone_sizes.sort(key=lambda x: x[2], reverse = True)
    for i in range(0,n_plot):
      if(clone_sizes[i][0] not in clones_plot):
        if(clone_sizes[i][2]>1):
          print ("\t", clone_sizes[i])
          clones_plot[clone_sizes[i][0]] = str(clone_sizes[i][0])+"_"+str(clone_sizes[i][1])+"_"+str(clone_sizes[i][2])
          out = out+pat+"\t"+clones_plot[clone_sizes[i][0]]+"\t"+str(clone_sizes[i][1])+"\t"+str(clone_sizes[i][2])+"\n"
  else:
    clone_sizes.sort(key=lambda x: x[2], reverse = True)
    for i in range(0,len(clone_sizes)):
      if(clone_sizes[i][2]>=2):
        out = out+pat+"\t"+clone_sizes[i][0]+"\t"+str(clone_sizes[i][1])+"\t"+str(clone_sizes[i][2])+"\n"
        clones_plot[clone_sizes[i][0]] = 1
  clone_summary = clone_alignment_files+"summary.txt"
  fh =open(clone_summary,"w")
  fh.write(out)
  fh.close()
  return()

def Get_clone_characteristics_between_samples(clone_alignment_files,merge_cluster_file_count):
  clone_pass_file = clone_alignment_files+"passed_clones_summary.txt"
  clones_plot_filtered = {}
  fh =open(clone_pass_file,"r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split() 
      clones_plot_filtered[l[1]] = l[1]+"_"+l[2]
  fh.close()
  count_uniq, count_all = {},{}
  for clone in clones_plot_filtered:
    split_fasta_file = clone_alignment_files+"split_filtered_"+clone+".fasta"
    fh=open(split_fasta_file,"r")
    for header,sequence in fasta_iterator(fh):
      if(header.count("__")!=0):
        f,clas = map(int, header.split("__")[1].split("|")[0].split("_")), header.split("|")[1].split("_")
        sample = header.split("|")[2]
        nz = [i for i in range(len(f)) if f[i]!=0]
        classes = {}
        for i in nz:
          c = clas[i].split("*")[0]
          name = sample+"\t"+clone+"\t"+c
          if(name in count_all):count_all[name] = count_all[name]+f[i]
          else:count_all[name]=f[i]
          classes[c] = 1
        for c in classes:
          name = sample+"\t"+clone+"\t"+c
          if(name in count_uniq):count_uniq[name] = count_uniq[name]+1
          else:count_uniq[name]=1
      else:
        sample = header.split("|")[2]
        name = sample+"\t"+clone+"\tAll"
        if(name in count_uniq):count_uniq[name] = count_uniq[name]+1
        else:count_uniq[name]=1
        if(name in count_all):count_all[name] = count_all[name]+1
        else:count_all[name] = 1
  out = "#sample\tclonotype\tisotype\tnumber of BCRs\tnumber of unique BCRs\n"
  for c1 in count_all:
    out = out +"\t".join([c1,str(count_all[c1]), str(count_uniq[c1])])+"\n"
  summary_file = clone_alignment_files+"clone_frequencies.txt"
  fh=open(summary_file,"w")
  fh.write(out)
  fh.close()
  return()

def Get_sequences_clones_single_chain(clones_plot, combined_sequences_file1, merge_cluster_file,clone_alignment_files):
  fh = open(merge_cluster_file,"r")
  seqs,clone_inverse = {},{}
  all_clones = {}
  for l in fh:
    if(l[0]!="#"):
      l = l.strip().split()
      pas = 0
      if(l[0] in clones_plot):
        if(l[1].count("-")==0):
          pas = 1
        else:pas = 1
        if(pas==1):
          name = clones_plot[l[0]]
          clone_inverse[l[1]] = name
          all_clones[name] = 1
  fh.close()
  chains,files = ["A"],[combined_sequences_file1]
  seqs = {}
  retained_clones = {}
  for i1 in range(0,len(files)):
    fh = open(files[i1],"r")
    for header,sequence in fasta_iterator(fh):
      header1 = header.split("__")[0].split("|")[0]+"|"+header.split("|")[len(header.split("|"))-1]
      if(header1 in clone_inverse):
        name = clone_inverse[header1]+"_"+chains[i1]
        if(name in seqs):seqs[name] = seqs[name] +">"+header+"\n"+sequence+"\n"
        else:seqs[name] = ">"+header+"\n"+sequence+"\n"
    fh.close()
  for clone in all_clones:
    seq_align_join = {}
    clone_filtered_file = clone_alignment_files+"_filtered_"+str(clone)+".fasta"
    sequences_aligned_trimmed1, sequences_aligned_trimmed2 ={},{}
    for i1 in range(0,len(files)):
      clone1 = clone+"_"+chains[i1]
      clone_raw_file = clone_alignment_files+chains[i1]+"_"+str(clone1)+".fasta"
      if(seqs[clone1].count(">")>5):
        fh =open(clone_raw_file,"w")
        fh.write(seqs[clone1])
        fh.close()
        clone_align_file = clone_alignment_files+"_"+chains[i1]+"_align_"+str(clone1)+".fasta"
        mafft = "/usr/local/bin/mafft "
        insert = ''
        command1 = mafft+" --retree 2 "+insert+" "+clone_raw_file+" > "+clone_align_file
        os.system(command1)
        if(i1 ==0):sequences_aligned_trimmed1 = Trim_sequences(clone_align_file)
        if(i1 ==1):sequences_aligned_trimmed2 = Trim_sequences(clone_align_file)
    clone_align_file = clone_alignment_files+"_"+chains[i1]+"_"+str(clone)+"_joined_aligned.fasta"
    print ("sequences_aligned_trimmed1",len(sequences_aligned_trimmed1), len(sequences_aligned_trimmed2))
    out = ''
    seq_out = {}
    for id in sequences_aligned_trimmed1:
      if(id in sequences_aligned_trimmed2):
        seq = sequences_aligned_trimmed1[id]+"NNN"+sequences_aligned_trimmed2[id]
        seq_out[seq]=1
        out=out+">"+id+"\n"+sequences_aligned_trimmed1[id]+"NNN"+sequences_aligned_trimmed2[id]+"\n"
      else:print (id,"NOT MATCHED", clone)
    if(len(out)>0):
      fh =open(clone_filtered_file,"w")
      fh.write(out)
      fh.close()
      retained_clones[clone] =str(out.count(">"))+"\t"+str(len(seq_out))
  out = "#pat\tclone\tsize\tn_samples\n"
  for c in retained_clones:
    out = out+"\t".join([pat, c,str(retained_clones[c])])+"\n"
  clone_pass_file = clone_alignment_files+"passed_clones_summary.txt"
  fh=open(clone_pass_file,"w")
  fh.write(out)
  fh.close()
  return()

def Get_and_align_largest_clones_single_chain (pat, output_dir, combined_sequences_file1, merge_cluster_file, merge_cluster_file_count,clone_alignment_files,fasttree):
  Initialise_clonal_analysis(pat, output_dir, merge_cluster_file, merge_cluster_file_count,clone_alignment_files,fasttree)
  clones_plot = {}
  clone_summary = clone_alignment_files+"summary.txt"
  fh =open(clone_summary,"r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split() 
      clones_plot[l[1]] = l[1]+"_"+l[2]
  fh.close()
  Get_sequences_clones_single_chain(clones_plot, combined_sequences_file1,  merge_cluster_file,clone_alignment_files)
  Generate_tree(clones_plot,clone_alignment_files,fasttree)
  return()

def Get_and_align_largest_clones (pat, output_dir, combined_sequences_file1, combined_sequences_file2, merge_cluster_file, merge_cluster_file_count,clone_alignment_files,fasttree):
  Initialise_clonal_analysis(pat, output_dir, merge_cluster_file, merge_cluster_file_count,clone_alignment_files,fasttree)
  clones_plot = {}
  clone_summary = clone_alignment_files+"summary.txt"
  fh =open(clone_summary,"r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      clones_plot[l[1]] = l[1]+"_"+l[2]
  fh.close()
  Get_sequences_clones(clones_plot, combined_sequences_file1, combined_sequences_file2, merge_cluster_file,clone_alignment_files)
  ######Split_by_clones(clone_alignment_files, clones_plot)
  Generate_tree(clones_plot,clone_alignment_files,fasttree)
  #Get_clone_characteristics_between_samples(clone_alignment_files,merge_cluster_file_count)
  return()

def Split_by_clones(clone_alignment_files, clones_plot):
  passes_clones= "#clone\tsplit_clone\tn_sequences\tn_samples\n"
  for clone1 in clones_plot:
    clone = clones_plot[clone1]
    clone_filtered_file = clone_alignment_files+"filtered_"+str(clone)+".fasta"
    clone_split_file = clone_alignment_files+"split_"+str(clone)+".fasta"
    exists = os.path.isfile(clone_filtered_file)
    if exists:
      #clone_tree_file = clone_alignment_files+"tree_"+str(clone)+".txt"
      fh= open(clone_filtered_file,"r")
      ### generate network
      G=nx.Graph()
      seqs,ids,seq_dict = [],[],{}
      for header,sequence in fasta_iterator(fh):
        seqs,ids =seqs+[sequence],ids+[header]
        seq_dict[header] = sequence
      fh.close()
      lens = len(seqs[0])
      for i in range(0,len(ids)):
        for j in range(i, len(ids)):
          if(i<j):
            diff = len([i1 for i1 in range(lens) if seqs[i][i1]!=seqs[j][i1]])
            if(diff<=30):G.add_edge(ids[i],ids[j])
      con= nx.connected_components(G)
      unique_cluster_ind = 0
      for c in con:
        unique_cluster_ind = unique_cluster_ind+1
        ids_cluster,samples_cluster = [], {}
        for id in c:
          ids_cluster.append(id)
          ids1 = id.split("|")
          samples_cluster[ids1[len(ids1)-1]]=1
        if(len(samples_cluster)>1):
          split_clone_name = clone+"_"+str(unique_cluster_ind)
          passes_clones = passes_clones+clone+"\t"+split_clone_name+"\t"+str(len(ids_cluster))+"\t"+str(len(samples_cluster))+"\n"
          clone_split_file = clone_alignment_files+"split_filtered_"+str(split_clone_name)+".fasta"
          out = ''
          for id in ids_cluster:
            out=out+">"+id+"\n"+seq_dict[id]+"\n"
          fh=open(clone_split_file,"w")
          fh.write(out)
          fh.close()
          del out
  clone_pass_file = clone_alignment_files+"passed_clones_summary.txt"
  fh=open(clone_pass_file,"w")
  fh.write(passes_clones)
  fh.close()
  return()

def Get_sequences_clones(clones_plot, combined_sequences_file1, combined_sequences_file2, merge_cluster_file,clone_alignment_files):
  fh = open(merge_cluster_file,"r")
  seqs,clone_inverse = {},{}
  all_clones = {}
  for l in fh:
    if(l[0]!="#"):
      l = l.strip().split()
      pas = 0
      if(l[0] in clones_plot):
        if(l[1].count("-")==0):
          pas = 1
        else:pas = 1
        if(pas==1):
          name = clones_plot[l[0]]
          clone_inverse[l[1]] = name
          all_clones[name] = 1
  fh.close()
  chains,files = ["A","B"],[combined_sequences_file1, combined_sequences_file2]
  chains,files = ["A"],[combined_sequences_file1]
  seqs = {}
  retained_clones = {}
  for i1 in range(0,len(files)):
    fh = open(files[i1],"r")
    for header,sequence in fasta_iterator(fh):
      header1 = header.split("__")[0].split("|")[0]+"||"+header.split("|")[len(header.split("|"))-1]
      if(header1 in clone_inverse):
        name = clone_inverse[header1]+"_"+chains[i1]
        if(name in seqs):seqs[name] = seqs[name] +">"+header+"\n"+sequence+"\n"
        else:seqs[name] = ">"+header+"\n"+sequence+"\n"
    fh.close()
  for clone in all_clones:
    seq_align_join = {}
    clone_filtered_file = clone_alignment_files+"_filtered_"+str(clone)+".fasta"
    sequences_aligned_trimmed1, sequences_aligned_trimmed2 ={},{}
    for i1 in range(0,len(files)):
      clone1 = clone+"_"+chains[i1]
      clone_raw_file = clone_alignment_files+chains[i1]+"_"+str(clone1)+".fasta"
      if(seqs[clone1].count(">")>5):
        fh =open(clone_raw_file,"w")
        fh.write(seqs[clone1])
        fh.close()
        clone_align_file = clone_alignment_files+"_"+chains[i1]+"_align_"+str(clone1)+".fasta"
        mafft = "/apps/well/mafft/7.149/bin/mafft "
        insert = ''
        command1 = mafft+" --retree 2 "+insert+" "+clone_raw_file+" > "+clone_align_file
        os.system(command1)
        if(i1 ==0):sequences_aligned_trimmed1 = Trim_sequences(clone_align_file)
        if(i1 ==1):sequences_aligned_trimmed2 = Trim_sequences(clone_align_file)
    clone_align_file = clone_alignment_files+"_"+chains[i1]+"_"+str(clone)+"_joined_aligned.fasta"
    print ("sequences_aligned_trimmed1",len(sequences_aligned_trimmed1), len(sequences_aligned_trimmed2))
    out = ''
    seq_out = {}
    for id in sequences_aligned_trimmed1:
      if(id in sequences_aligned_trimmed2):
        seq = sequences_aligned_trimmed1[id]+"NNN"+sequences_aligned_trimmed2[id]
        seq_out[seq]=1
        out=out+">"+id+"\n"+sequences_aligned_trimmed1[id]+"NNN"+sequences_aligned_trimmed2[id]+"\n"
      else:print (id,"NOT MATCHED", clone)
    if(len(out)>0):
      fh =open(clone_filtered_file,"w")
      fh.write(out)
      fh.close()
      retained_clones[clone] =str(out.count(">"))+"\t"+str(len(seq_out))
  out = "#pat\tclone\tsize\tn_samples\n"
  for c in retained_clones:
    out = out+"\t".join([pat, c,str(retained_clones[c])])+"\n"
  clone_pass_file = clone_alignment_files+"passed_clones_summary.txt"
  fh=open(clone_pass_file,"w")
  fh.write(out)
  fh.close()
  return()

def Trim_sequences(clone_align_file):
  sequences_aligned_trimmed = {}
  fh =open(clone_align_file,"r")
  start,end = -1,1000000000
  seqs_sub = {} 
  starts,ends = [],[] 
  for header,sequence in fasta_iterator(fh):
    if(sequence.count("-")==0):
      start1,end1 = 0,len(sequence)
      starts,ends =starts+[start1],ends+[end1]
      seqs_sub[header]=[sequence, start1,end1]
    else: 
      for i in range(0,50):
        if(sequence[i]!="-"):
          break
      start1=i
      len_s = len(sequence)
      for i in range(0,50):
        if(sequence[len_s-i-1]!="-"):
          break
      end1 = len_s-i
      starts,ends =starts+[start1],ends+[end1]
      seqs_sub[header]=[sequence, start1,end1]
  fh.close()
  print (starts, ends)
  unique_starts = np.unique(starts)
  counts_starts = [starts.count(unique_starts[x]) for x in range(len(unique_starts))]
  cumulative = 0
  for i in range(0,len(unique_starts)):
    cumulative = cumulative+counts_starts[i]*1.0/sum(counts_starts)
    if(cumulative>=0.95):break
  start = unique_starts[i]
  unique_starts = np.unique(ends)#.tolist
  unique_starts = unique_starts[::-1]
  counts_starts = [ends.count(unique_starts[x]) for x in range(len(unique_starts))]
  cumulative = 0
  for i in range(0,len(unique_starts)):
    cumulative = cumulative+counts_starts[i]*1.0/sum(counts_starts)
    if(cumulative>=0.95):break
  end = unique_starts[i]
  out = ''
  removed,kept = 0,0
  for id in seqs_sub:
    seq = seqs_sub[id][0][start:end].upper()
    sequences_aligned_trimmed[id] = seq
    kept = kept+1
  return(sequences_aligned_trimmed)

def Generate_tree(clones_plot,clone_alignment_files,fasttree):
  clone_pass_file = clone_alignment_files+"passed_clones_summary.txt"
  fh=open(clone_pass_file,"r")
  clones_plot1={}
  for l in fh:
    if(l[0]!="#"):
      l = l.strip().split()
      if(int(l[3])>1):
        clones_plot1[l[1]] = 1
  for clone in clones_plot1:
    #clone = clones_plot1[clone1]
    clone_filtered_file = clone_alignment_files+"_filtered_"+str(clone)+".fasta"
    clone_tree_file = clone_alignment_files+"tree_"+str(clone)+".txt"
    command1 = fasttree+"  -nt -boot 1000 -gtr -log "+clone_filtered_file+"LOG < "+clone_filtered_file+" > "+clone_tree_file
    os.system(command1)
    print (clone)
  return()

class Tree(defaultdict):
  def __init__(self, value=None):
    super(Tree, self).__init__(Tree)
    self.value = value

def Get_info(file):
  fh=open(file,"r")
  groups,group_info,group_order = Tree(),{},[]
  for l in fh:
    if(l[0]!="#"):
      l = l.strip().split()
      group_id,sample,annot_file,seq_file,output_dir,annot_type = l[0],l[1],l[2],l[3],l[4],l[5]
      groups[group_id][sample].value = 1
      group_info[sample] = [group_id,sample,annot_file,seq_file,output_dir,annot_type]
      if(group_id not in group_order):group_order.append(group_id)
  fh.close()
  for g in groups: 
      for sample in groups[g]:
        print (g, sample)
  return(groups,group_info,group_order)

###########################
if(len(sys.argv)!=4):
    print ("\nUsage:\npython3 MRDARCY_2.2.py <input sample id file>  <sample number from input file> <command number>\n")
    print ("<input sample id file> should be a tab delimited file with the following columns:")
    print ("\tColumn 1: Group name\n\tColumn 2: Sample name\n\tColumn 3: Annotation file location\n\tColumn 4: Fasta file location\n\tColumn 5: Output directory\n\tColumn 6: Annotation file type (IMGT or ImmuneReceptorAnnotator or Single_Cell)\n")
    print ("<command number>:\n\t1: Merge files\n\t2: Align sequences from largest clones\n")
    print ("e.g.\npython3 MRDARCY_2.0.py Sample_example.txt 1 1\n\n")
    print ("Edit line 1152 to point to cd-hit location\n\n")
    print ("Note that the single cell version needs to be updated for non-Oxford cluster usage")
    quit()

infile = sys.argv[1]                 ### full location name for input file
index_required = int(sys.argv[2])-1  ### index of sample group to run
cd_hit_directory = "/Users/sammut01/software/cd-hit-v4.8.1-2019-0228/"
command_ind = int(sys.argv[3])
###########################
groups,group_info,group_order =Get_info(infile)

if(index_required>len(group_order)-1):
  print ("\tError: index too high for sample file")
  quit()
###########################
edge_lengths = 0.95
pat = group_order[index_required]
seq_files, annot_files,samples,annot_type,output_dir = [],[],[],'',''
for sample in groups[pat]: 
  annot_file,seq_file,output_dir,annot_type = group_info[sample][2],group_info[sample][3],group_info[sample][4],group_info[sample][5]
  seq_files, annot_files,samples = seq_files+[seq_file], annot_files + [annot_file], samples+[sample]

combined_sequences_file1 =output_dir+"Combined_fully_reduced_A_sequences_"+pat+".fasta"
CDR3_info_file = output_dir+"CDR3_information_"+pat+".txt"
merge_cluster_file = output_dir+"Merge_clustering_"+pat+".txt"
tmp_file = output_dir+"Tmp_A_"+pat+"_"
merge_cluster_file_count = output_dir+"Count_clusters_merged_"+pat+".txt"
summary_merge_clustering_count= output_dir+"Count_per_sample_clusters_merged_"+pat+".txt"
combined_sequences_file2 =output_dir+"Combined_fully_reduced_B_sequences_"+pat+".fasta"
tmp_file2 = output_dir+"Tmp_B_"+pat+"_"


if(command_ind ==1):
  if(annot_type=="IMGT"):
    CDR3_defined_sequences_IMGT(pat, output_dir, seq_files, annot_files,samples,cd_hit_directory,tmp_file,edge_lengths,merge_cluster_file_count,combined_sequences_file1,summary_merge_clustering_count)
  elif(annot_type=="10X"):
    CDR3_defined_sequences_SC(pat, output_dir, seq_files, annot_files,samples,cd_hit_directory,tmp_file1,tmp_file2,edge_lengths,merge_cluster_file_count,combined_sequences_file1,combined_sequences_file2)


if(command_ind ==2):
  clone_alignment_files = output_dir+"CLONES/Clone_alignment_"+pat+"_"
  fasttree = "/Users/sammut01/software/fast-tree/FastTree"
  if(annot_type=="IMGT"):
    Get_and_align_largest_clones_single_chain (pat, output_dir, combined_sequences_file1,  merge_cluster_file, merge_cluster_file_count,clone_alignment_files,fasttree)
  elif(annot_type=="10X"):
    Get_and_align_largest_clones (pat, output_dir, combined_sequences_file1, combined_sequences_file2, merge_cluster_file, merge_cluster_file_count,clone_alignment_files,fasttree)




