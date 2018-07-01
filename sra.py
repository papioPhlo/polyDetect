"""
Created on Fri Apr  8 14:31:28 2016

@author: vallmer
"""


import sys
import subprocess
from sam import *
from essentials import *
from meiDetect import *
import os

def numberFastq(fastq_list):
    """The purpose of this function is to take a fastq list and replace the read names with numbers. These numbers will be used
    later to sort the sam file reads after they have been aligned. The options input will be used to generate different output options
    if the user would simply like the numbered fastq file, then the enter the integer 0 as the option. If they would like the fastq,
    a list of numbers, and sequences"""
    i,final,number,length = 0,[],1,len(str(len(fastq_list)))
    while i < (len(fastq_list)-4):
        final.append("@" + zeroThehero(length, number))
        final.append(fastq_list[i+1])
        final.append("+")
        final.append(fastq_list[i+3])
        i+=4
        number+=1
    return(final)
        
