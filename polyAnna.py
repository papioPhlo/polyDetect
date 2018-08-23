#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import sys
import re
from sam import *
from essentials import *
from meiDetect import *
from sra import *
import linecache
import numpy as np
import os
from subprocess import Popen, PIPE
from timeit import default_timer as timer

            
def processSam_polyA(path2samP,path2samS,outfq,outif):
    """The Purpose is to: 1.) Cut sam and output fastq"""    
    with open(path2samP,'r') as samfileP:
        with open(path2samS,'r') as samfileS:
            with open(outfq,'w') as fastq:
                with open(outif,'w') as info:
                    x=1
                    for lines in samfileP:
                        if (lines[0]!='@'):
                            if (lines.split()[5]!='*') and ('H' not in lines.split()[5]):
                                s=sam(lines)
                                s_span,s_char=splitCigar(s.cigar)
                                i,right_UM=int(find_last(s_char,'M'))+1,0
                                while i<len(s_char):
                                    right_UM+=int(s_span[i])
                                    i+=1
                                if (s_char[-1]!='M') and (right_UM>=20):
                                    sequence,n=s.seq[(0-right_UM):],0
                                    while sequence[n]=='A':
                                        if (len(sequence)>n+1):
                                            n+=1
                                        else:
                                            break
                                    if len(sequence[n:])>=20:
                                        fastq.write('@'+str(zeroThehero(12,x))+'\n'+sequence[n:]+'\n+\n'+s.qual[(0-right_UM+n):]+'\n')
                                        info.write(s.seq+'\n')
                                        x+=1
                    for lines in samfileS:
                        if (lines[0]!='@'):
                            if (lines.split()[2]!='*') and ('H' not in lines.split()[5]):
                                s=sam(lines)
                                s_span,s_char=splitCigar(s.cigar)
                                i,right_UM=int(find_last(s_char,'M'))+1,0
                                while i<len(s_char):
                                    right_UM+=int(s_span[i])
                                    i+=1
                                if (s_char[-1]!='M') and (right_UM>=20):
                                    sequence,n=s.seq[(0-right_UM):],0
                                    while sequence[n]=='A':
                                        if (len(sequence)>n+1):
                                            n+=1
                                        else:
                                            break
                                    if len(sequence[n:])>=20:
                                        fastq.write('@'+str(zeroThehero(12,x))+'\n'+sequence[n:]+'\n+\n'+s.qual[(0-right_UM+n):]+'\n')
                                        info.write(s.seq+'\n')
                                        x+=1
                                
def predictPolyA(path2sam,outfile,path2info):
    infoList=lister(path2info)
    with open(path2sam,'r') as samfile:
        with open(outfile,'w') as output:
            final,chromosome_list=[],[]
            for i in range(1,21):
                final.append([])
                chromosome_list.append(str(i))
            for lines in samfile:
                if lines[0]!='@':
                    if (lines.split()[2]!='*') and ('H' not in lines.split()[5]):
                        s=sam(lines)
                        if s.rname[3:] in chromosome_list:
                            if s.orientation=='0':
                                final[int(s.rname[3:])-1].append(s.rname[3:]+'\t'+s.pos+'\t300\tA\tN\t000000000000:000000000000\t'+infoList[int(s.qname)-1]+'\t0\n')
                            else:
                                s_span,s_char=splitCigar(s.cigar)
                                sequence,n,p,mapped=infoList[int(s.qname)-1],find_last(s_char,'M'),0,0
                                while p<=n:
                                    mapped+=int(s_span[p])
                                    p+=1
                                final[int(s.rname[3:])-1].append(s.rname[3:]+'\t'+str(int(s.pos)+mapped)+'\t300\tA\tN\t000000000000:000000000000\tAAAAAAAAAAAAA'+sequence+'\t1\n')
            for f in final:
                f_sort=sortMe(f,2)
                for line in f_sort:
                    output.write(line+'\n')
                    
#processSam1(sys.argv[1],sys.argv[2],sys.argv[3])
#bowtie2(sys.argv[1],sys.argv[2],sys.argv[3],'S')
#predictPolyA(sys.argv[1],sys.argv[2],sys.argv[3])
