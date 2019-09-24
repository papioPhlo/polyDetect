#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 13:45:57 2019

@author: batzerlab
"""
import re, sys

with open(sys.argv[1],'r') as infile:
    with open(sys.argv[2],'w') as outfile:
        for line in infile:
            spl=line.split()
            seq=spl[6]
            alu_head=re.search(r"(G|C)CCGGG(C|T)(A|G|C)(C|T)(G|T|A)",seq)
            #if int(spl[2]) <= 5:
             #   if spl[4]=='S' or spl[4]=='P':
            if alu_head:
                flanking_len=len(seq[0:alu_head.start()])
                if flanking_len >= 30:
                    outfile.write(line)