"""
Created on Mon Feb 15 09:55:28 2016
@author: vallmer
"""

def orientConverter(FLAG):
    if FLAG != "0":
        binary = bin(int(FLAG))[2:]
        return(str(binary[-5]))
    else:
        return("0")

class sam:
    "A simple class for reads in SAM format"
    def __init__(self, line):
        self.qname=line.split()[0]
        self.flag=line.split()[1]
        self.rname=line.split()[2]
        self.pos=line.split()[3]
        self.mapq=line.split()[4]
        self.cigar=line.split()[5]
        self.rnext=line.split()[6]
        self.pnext=line.split()[7]
        self.tlen=line.split()[8]
        self.seq=line.split()[9]
        self.qual=line.split()[10]
        self.orientation=orientConverter(line.split()[1])
        
class info:
    "A simple class for info lines"
    def __init__(self, line):
        self.neWname=line.split()[0]
        self.olDname=line.split()[1]
        self.pos=line.split()[2]
        self.cigar=line.split()[3]
        self.seq=line.split()[4]
        self.flag=line.split()[5]
        self.orientation=line.split()[6]
        
class predict:
    "A simple class for predict lines"
    def __init__(self, line):
        self.chromosome=line.split(':')[0][3:]
        self.locus=line.split()[0].split(':')[1]
        self.alu_position=line.split()[1]
        self.pos_tag=line.split()[2]
        self.map_tag=line.split()[3]
        self.new_name=line.split()[4].split(':')[0]
        self.old_name=line.split()[4].split(':')[1]
        self.orientation=line.split()[5]

class predictSort:
    "A simple class for predict lines"
    def __init__(self, line):
        self.chromosome=line.split()[0]
        self.locus=line.split()[1]
        self.alu_position=line.split()[2]
        self.pos_tag=line.split()[3]
        self.map_tag=line.split()[4]
        self.new_name=line.split()[5].split(':')[0]
        self.old_name=line.split()[5].split(':')[1]
        self.orientation=line.split()[6]
        
class rptmsk:
    "A simple class for RepeatMasker output"
    def __init__(self, line):
        self.SWscore=line.split()[0]
        self.perc_div=line.split()[1]
        self.perc_del=line.split()[2]
        self.perc_ins=line.split()[3]
        self.chrom=line.split()[4][3:]
        self.locus_start=line.split()[5]
        self.locus_end=line.split()[6]
        self.orient=line.split()[8]
        self.repeat=line.split()[9]
        self.repeat_end=line.split()[12]
        if self.orient=='+':
            self.repeat_start=line.split()[11]
            self.remainder=line.split()[13]
            self.head=line.split()[5]
            self.tail=line.split()[6]
        else:
            self.repeat_start=line.split()[13]
            self.remainder=line.split()[11]
            self.head=line.split()[6]
            self.tail=line.split()[5]
        self.family=line.split()[10]
        self.id=line.split()[14]
        
def nGenerator(number):
    Ns='NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
    return(Ns[:number])
        
    

import re
        
def sam2fastq(samList):
    fastq = []
    for sam in samList:
        sam_split = re.split('\t|_',sam)
        if sam_split[2] != "0":
            if str(bin(int(sam_split[2]))[-5]) == "1":
                fastq.append("@" + sam_split[3] + "_" + sam_split[4] + "_" + sam_split[6] + "_" + sam_split[2] + "\n" + reverseComplement(sam_split[1]) + "\n+\n" + "HHHHEGHHHHHGHHHHFHHGHHHHHHHHHHDGGGGHHHHFHHHGHHHHHGEGEHGHHFHHHHHHBHHHHHHHHHHFGFGGGHHHHHEHHHHHHHHGFHHEHHHHHHHGHHHHHFGHHHHHHCHHGHFHHHHHHHHGHHHFHHHGHCFHFHGHGHHHHHHHHHHHHEDHHFGFBG?FEGEFFCD@9EDEEE<CECG@7B;A")
        else:
            fastq.append("@" + sam_split[3] + "_" + sam_split[4] + "_" + sam_split[6] + "_" + sam_split[2] + "\n" + sam_split[1] + "\n+\n" + "HHHHEGHHHHHGHHHHFHHGHHHHHHHHHHDGGGGHHHHFHHHGHHHHHGEGEHGHHFHHHHHHBHHHHHHHHHHFGFGGGHHHHHEHHHHHHHHGFHHEHHHHHHHGHHHHHFGHHHHHHCHHGHFHHHHHHHHGHHHFHHHGHCFHFHGHGHHHHHHHHHHHHEDHHFGFBG?FEGEFFCD@9EDEEE<CECG@7B;A")
    return(fastq)

def find_last(lst, sought_elt):
    for r_idx, elt in enumerate(reversed(lst)):
        if elt == sought_elt:
            return(len(lst) - 1 - r_idx)

def splitCigar(cigar):
    list_of_characters,index = ['M','I','D','N','S','H','P','X'],[]   
    for c in list_of_characters:
        ind = [m.start() for m in re.finditer(c,str(cigar))]
        for i in ind:
            index.append(str(i))
    sorted_indices = sorted(index,key=int)
    cigar_span_list, cigar_character_list = [],[]
    left = 0
    for s in sorted_indices:
        right = int(s)
        cigar_span_list.append(cigar[left:right])
        cigar_character_list.append(cigar[int(s)])
        #print(cigar[left:right])
        #print(cigar[int(s)])
        left = right + 1
    return(cigar_span_list,cigar_character_list)
    
def flipCigar(cigar_span_list,cigar_character_list):
    flipped = ""
    cigar_span_rv,cigar_character_rv = cigar_span_list[::-1],cigar_character_list[::-1]
    for i in range(len(cigar_span_rv)):
        flipped+=(cigar_span_rv[i] + cigar_character_rv[i])
    return(flipped)
    
from essentials import *

def tupleCigarSpans(cgr_spn):
    final, length, i, start, max_num = [], len(cgr_spn), 0, 1, int(cgr_spn[0])
    while i < length:
        spanned = []
        while start <= max_num:
            spanned.append(start)
            start +=1
        final.append(spanned)
        i+=1
        if i < length:
            max_num += int(cgr_spn[i])    
    return(final)        
    
def mListcombine(indexlist,spanlist):
    final = []
    for i in indexlist:
        final.extend(spanlist[i])
    return(final)

def shared(sam):
    """"This function checks the portion of the read that maps in the reference against the portion that maps to the Alu. 
                    REMEMBER:::::: The QNAME must be formatted chr#_leftmostref_cigarstr_..."""
    read = samRead()
    split = read.inputRead(sam)
    qname = read.QNAME(split) 
    cgr_spn_Alu,cgr_char_Alu = splitCigar(read.CIGAR(split))
    splitQ = qname.split("_")
    cgr_spn_ref,cgr_char_ref = splitCigar(splitQ[2])
    flag = read.orientation(split[1])
    ref_spanned,alu_spanned = tupleCigarSpans(cgr_spn_ref),tupleCigarSpans(cgr_spn_Alu)
    if flag == "0":
        M_index_ref,M_index_alu = allOccurences("M",cgr_char_ref), allOccurences("M", cgr_char_Alu)
        M_span_ALU, M_span_REF = tuple(mListcombine(M_index_alu,alu_spanned)), tuple(mListcombine(M_index_ref,ref_spanned))
        bpMapped = common_elements(M_span_ALU, M_span_REF)
        return(len(bpMapped) / len(M_span_ALU))
    elif flag == "1":
        alu_spanned = tupleCigarSpans(cgr_spn_Alu[::-1])
        M_index_ref,M_index_alu = allOccurences("M",cgr_char_ref), allOccurences("M", (cgr_char_Alu[::-1]))
        M_span_ALU, M_span_REF = tuple(mListcombine(M_index_alu,alu_spanned)), tuple(mListcombine(M_index_ref,ref_spanned))
        bpMapped = common_elements(M_span_ALU, M_span_REF)
        return(len(bpMapped) / len(M_span_ALU))

def distanceM(cigar):
    cigar_span,cigar_char = splitCigar(cigar)
    index_1st_M = cigar_char.index("M")
    i = index_1st_M -1    
    dist = 0
    while i > -1:
        dist += int(cigar_span[i])
        i -= 1
    return(dist) 

def flagRev(flag):
    if flag == "0":
        return("1")
    elif flag == "1":
        return("0")
import sys 

def aluPredict_old(samList):
    final = []
    for sam in samList:
        sam_split = re.split('\t|_',sam)
        num,ref_chr,ref_start,ref_cig,alu_start,alu_cig = sam_split[0],sam_split[6],sam_split[7],sam_split[9],sam_split[3],sam_split[4]
        dist_ref = distanceM(ref_cig)
        if sam_split[1] == "0":
            ref_tag = "0"
        else:
            ref_tag = str(bin(int(sam_split[1])))[-5]
        if sam_split[5] == "0":
            alu_tag = "0"
        else:
            alu_tag = str(bin(int(sam_split[5])))[-5]
        if ref_tag == "0":
            if alu_tag == "0":
                dist_Alu = distanceM(alu_cig)
                start = str(int(ref_start) - int(dist_ref) + int(dist_Alu) - int(alu_start))
                orientation = "0"
            elif alu_tag == "1":
                alu_span,alu_char = splitCigar(alu_cig)
                index_M_Alu = alu_char.index("M")
                alu_cig_flipped = flipCigar(alu_span,alu_char)
                dist_Alu = distanceM(alu_cig_flipped)
                start = str(int(ref_start) - int(dist_ref) + int(alu_span[index_M_Alu]) + int(dist_Alu) + int(alu_start))
                orientation = "1"
        elif ref_tag  == "1":
            if alu_tag == "0":
                alu_span,alu_char = splitCigar(alu_cig)
                index_M_Alu = alu_char.index("M")
                alu_cig_flipped = flipCigar(alu_span,alu_char)
                dist_Alu = distanceM(alu_cig_flipped)
                start = str(int(ref_start) - int(dist_ref) + int(alu_span[index_M_Alu]) + int(dist_Alu) + int(alu_start))
                orientation = "1"
            elif alu_tag == "1":
                dist_Alu = distanceM(alu_cig)
                start = str(int(ref_start) - int(dist_ref) + int(dist_Alu) - int(alu_start))
                orientation = "0"
        final.append(ref_chr + "\t" + start + "\t" + orientation + "\t" + num)
    return(final)
    
def chopSeq(cigar,seq_len):
    read_len = (int(seq_len)/2)
    span_list,char_list = splitCigar(cigar)
    cig_span = 0
    span_a,span_b,char_a,char_b = [],[],[],[]
    i=0
    while cig_span < read_len:
        cig_span+=int(span_list[i])
        if cig_span <= read_len:
            span_a.append(span_list[i])
            char_a.append(char_list[i])
            i+=1
        else:
            if (cig_span-int(span_list[i])) < read_len:
                span_a.append(str(int(read_len)-(cig_span-int(span_list[i]))))
                char_a.append(char_list[i])
                span_b.append(str(int(span_list[i]) - (int(read_len)-(cig_span-int(span_list[i])))))
                char_b.append(char_list[i])
                i+=1
    while i < len(span_list):
        span_b.append(span_list[i])
        char_b.append(char_list[i])
        i+=1
    if "M" in char_a:
        return(span_a,char_a,int(read_len))
        #print(char_a)
    elif "M" in char_b:
        return(span_b,char_b,int(read_len))
        #print(char_b)
    
def returnMappedBlocks(cigar,locus,seq_len):
    mapped_block_list = []
    span_list,char_list = splitCigar(cigar)
    range_start = char_list.index("M")
    position = 0
    if range_start!= 0:
        i = range_start -1
        while i >= 0:
            position+=int(span_list[i])
            i-=1        
    start = int(locus)
    for i in range(range_start,len(char_list)):
        if char_list[i] == "M":
            end = start + int(span_list[i]) - 1
            mapped_block_list.append(str(start) + "_" + str(end) + "_" + str(seq_len - position - int(span_list[i]) + 1))
            position = position + int(span_list[i]) 
            start = end + 1
        elif char_list[i] == "S":                
            end = start + int(span_list[i]) - 1
            start = end + 1
            position = position + int(span_list[i])
    return(mapped_block_list)
    
def returnMappedBlocksB(span_list,char_list,locus,seq_len):
    mapped_block_list = []
    range_start = char_list.index("M")
    position = 0
    if range_start!= 0:
        i = range_start -1
        while i >= 0:
            position+=int(span_list[i])
            i-=1        
    start = int(locus)
    for i in range(range_start,len(char_list)):
        if char_list[i] == "M":
            end = start + int(span_list[i]) - 1
            mapped_block_list.append(str(start) + "_" + str(end) + "_" + str(seq_len - position - int(span_list[i]) + 1))
            position = position + int(span_list[i]) 
            start = end + 1
        elif char_list[i] == "S":                
            end = start + int(span_list[i]) - 1
            start = end + 1
            position = position + int(span_list[i])
    return(mapped_block_list)
    
def zygoPercent_return(block_list,orientation,locus,window_start,window_end):
    positive_blocks,total = 0,0
    if orientation == "0":
        for block in block_list:
            block_split = block.split("_")
            start,end,remaining = block_split[0],block_split[1],block_split[2]
            if (int(start) <= (int(locus) - int(window_start))) and ((int(locus) + int(window_end)) <= (int(end) + int(remaining))):
                total+=1
                if (int(locus) + int(window_end)) <= int(end):
                    positive_blocks +=1
    elif orientation == "1":
        for block in block_list:
            block_split = block.split("_")
            start,end,remaining = block_split[0],block_split[1],block_split[2]
            if (int(start) <= (int(locus) - int(window_end))) and ((int(locus) + int(window_start)) <= (int(end) + int(remaining))):
                total+=1
                if (int(locus) + int(window_start)) <= int(end):
                    positive_blocks +=1
    if total == 0:
        return("none\tnone\tnone")
    else:
        return(str(positive_blocks/total) + "\t" + str(positive_blocks) + "\t" + str(total))
                    
        
def predictStart(sam):
    read = samRead()
    split = read.inputRead(sam)
    qname = read.QNAME(split)
    refName = read.RNAME(split)
    flag = read.orientation(split[1])
    locus_TAG = qname.split("_")
    ref_flag = read.orientation(locus_TAG[3])
    chr_NUM, leftmost_base_ref, leftmost_base_Alu = locus_TAG[0][3:], locus_TAG[1], read.POS(split)
    cigar_span_list, cigar_character_list = splitCigar(read.CIGAR(split))
    if ref_flag == "1":
        if cigar_character_list.count('M') == 1:
            cigar_span_list, cigar_character_list = splitCigar(read.CIGAR(split))
            flipped = flipCigar(cigar_span_list,cigar_character_list)
            cigar_span_list, cigar_character_list = splitCigar(flipped)
            index_M_Alu = cigar_character_list.index("M")
            mSpan = int(cigar_span_list[index_M_Alu])
        else:            
            cigar_span_list, cigar_character_list = splitCigar(read.CIGAR(split))
            flipped = flipCigar(cigar_span_list,cigar_character_list)
            cigar_span_list, cigar_character_list = splitCigar(flipped)
            first,last = cigar_character_list.index("M"),''.join(cigar_character_list).rfind("M")
            mSpan = 0
            for i in range(first,last+1):
                mSpan+=int(cigar_span_list[i])
        if flag == "0":
            leftmost_base_Alu = int(leftmost_base_Alu)
            leftmost_base_Alu+=mSpan
            leftmost_base_Alu = str(leftmost_base_Alu)
        elif flag == "1":
            leftmost_base_Alu = int(leftmost_base_Alu)
            leftmost_base_Alu-=mSpan
            leftmost_base_Alu = str(leftmost_base_Alu)
        flag = flagRev(flag)            
    if flag == "0":
        dist = distanceM(read.CIGAR(split))
        return((str(int(leftmost_base_ref) + int(dist) - int(leftmost_base_Alu))) + "\t" + refName + "\t" + flag + "\t" + str(shared(sam)))
    elif flag == "1":
        if shared(sam) > .5:
            cigar_character_list,cigar_span_list = (cigar_character_list[::-1]),(cigar_span_list[::-1])
            total_seq_length = len(read.SEQ(split))
            index_1st_M = cigar_character_list.index("M")
            i = index_1st_M - 1
            distance = 0
            mappedAlu = cigar_span_list[index_1st_M]
            while i>-1:
                distance+=int(cigar_span_list[i])
                i-=1
            dist_ref = distanceM(locus_TAG[2])
            return(str((int(leftmost_base_ref) + int(leftmost_base_Alu) + int(distance) + int(mappedAlu))) + "\t" + refName + "\t" + flag + "\t" + str(shared(sam)))
        else:
            ref_span_list, ref_char_list = splitCigar(locus_TAG[2])
            if ref_char_list[0] == "M":
                return(str(leftmost_base_ref) + "\t" + refName + "\t" + flag + "\t" + str(shared(sam)))
            else:
                i = ref_char_list.index("M")
                distance = ref_span_list[i]
                return((str(int(leftmost_base_ref) + int(distance))) + "\t" + refName + "\t" + flag + "\t" + str(shared(sam)) + "\t*")
            
def samPredict(path2samFile):
    predictionList = []
    samList = lister(path2samFile)
    for sam in samList:
        predictionList.append(predictStart(sam))
    return(predictionList)
    
def filterTAG(samList,TAGlist):
    """filters samfile based on cigar string. TAGlist is a list of unwanted characters"""
    final = []
    for sam in samList:
        read = samRead()
        split = read.inputRead(sam)
        cigar_span_list, cigar_character_list = splitCigar(read.CIGAR(split))
        present = []
        for i in TAGlist:
            if i in cigar_character_list:
                present.append("yes")
        if "yes" not in present:
            final.append(sam)
    return(final)       
