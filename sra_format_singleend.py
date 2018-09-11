"""
Created on Mon Sep  3 13:11 2018

@author: jessica_mfs
"""


import sys
from J_essentials import *

mylist=open(sys.argv[1]).readlines()
fastq_length=len(mylist)
with open(sys.argv[2],'w') as outfile:
    i,final,number,length = 0,[],1,12
    while i < (fastq_length):
        final.append("@" + zeroThehero(length, number)+'\t'+'OP:i:2'+'\n')
        final.append(mylist[i+1])
        final.append("+"+'\n')
        final.append(mylist[i+3])
        i+=4
        number+=1
            
    for line in final:
        outfile.write(line)
