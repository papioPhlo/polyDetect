# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 16:48:39 2016

@author: vallmer
"""

import sys

from sam import *
from essentials import *

def refCheck(refList,rCol,qryList,qCol):
    rLen,qLen,r,q,shared,not_shared = len(refList),len(qryList),0,0,[],[]
    while (r < rLen) and (q < qLen):
        rsplit,qsplit = refList[r].strip("\n").split("\t"),qryList[q].strip("\n").split("\t")
        r_locus,q_locus = int(rsplit[int(rCol)-1]),int(qsplit[int(qCol)-1])
        diff,abs_diff = r_locus - q_locus,abs(r_locus - q_locus)
        if abs_diff <= 15:
            shared.append(refList[r])
            r+=1
            q+=1
        elif diff < 0:
            not_shared.append(refList[r])
            r+=1
        elif diff > 0:
            q+=1
    if q >= qLen:
        while r < rLen:
            not_shared.append(refList[r])
            r+=1
    return(shared,not_shared)
