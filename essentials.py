"""
Created on Thu Feb 18 12:35:31 2016

@author: vallmer
"""

import os
import sys
import subprocess
#from twilio.rest import TwilioRestClient
from sam import *
from subprocess import Popen, PIPE
from timeit import default_timer as timer

def lister(path2file):
    with open(path2file, "r") as file:
        query = file.read().splitlines()
        return(query)
        
def lister_add(path2file,string2add): 
    final = []       
    with open(path2file, "r") as file:
        query = file.read().splitlines()
        for line in query:
            query_split = line.split()
            final.append(query_split[0] + "\t" + string2add + "\t" + query_split[2] + "\t" + query_split[3])
    return(final)
        
def kinkos(Qlist):
    for item in Qlist:
        print(item)
        
def allOccurences(char,string):
    length,i,final = len(string),0,[]
    while i < length:
        if char == string[i]:
            final.append(i)
        i+=1
    return(final)
    
def common_elements(list1, list2):
    return(list(set(list1) & set(list2)))

def tupler(qlist):
    qtuples = []
    for q in qlist:
        q_cols = []
        q_split = q.strip("\n").split("\t")
        for s in q_split:
            q_cols.append(s)
        tuple(q_cols)
        qtuples.append(q_cols)
    return(qtuples)
    
def tuplerB(qlist):
    qtuples = []
    for q in qlist:
        q_cols = []
        q_split = q.strip("\n").split("\t")
        col1,col2 = q_split[0],q_split[1]
        q_cols.append(col1)
        q_cols.append(col2)
        tuple(q_cols)
        qtuples.append(q_cols)
    return(qtuples)

def untupler(qtuples):
    finaList=[]
    for q in qtuples:
        q_list=list(q)
        finaList.append('\t'.join(q_list))
    return(finaList)
    
def untuplerB(qtuples):
    finaList=[]
    for q in qtuples:
        col1,col2 = q[0],q[1]
        newline = str(col1) + "\t" + str(col2)  
        finaList.append(newline)
    return(finaList)

def sortMe(qlist,col):
    qtuples = tupler(qlist)
    qlist_sort = untupler(sorted(qtuples,key=lambda support:(int(support[(int(col)-1)]))))
    return(qlist_sort)
    
def sortMeB(qlist):
    qtuples = tuplerB(qlist)
    qlist_sort = untuplerB(sorted(qtuples,key=lambda support:(int(support[1]))))
    return(qlist_sort)
    
def joinLists(listA,listB):
    final = []
    for item in listA:
        final.append(item)
    for item in listB:
        final.append(item)
    return(final)
    
def splitBYorient(column,fyList):
    """Simple function to split a file based on orientation. Input orientation column and listed file"""
    forward,reverse=[],[]    
    for line in fyList:
        split = line.strip("\n").split("\t")
        if split[int(column)-1] == "0":
            forward.append(line)
        elif split[int(column)-1] == "1":
            reverse.append(line)
    return(forward,reverse)    

def filterFile(string,column,fyList,presence):
    """The purpose of this function is to remove lines from a file that either contain or do not contain a
    particular string in a particular column. Yes or no."""
    final = []
    if presence == "no":
        for line in fyList:
            line_split = line.strip("\n").split("\t")
            if line_split[int(column) - 1] != string:
                final.append(line)
        return(final)
    elif presence == "yes":
        for line in fyList:
            line_split = line.strip("\n").split("\t")
            if line_split[int(column) - 1] == string:
                final.append(line)
        return(final)
        
def repeatOffenders_only(fyList):
    i,length,final = 0,len(fyList),[]
    while i < (length-1):
        if i == 0:
            present_split, future_split = fyList[i].strip("\n").split("\t"),fyList[i+1].strip("\n").split("\t")
            present_predict, future_predict = int(present_split[1]),int(future_split[1])
            if present_predict == future_predict:
                final.append(fyList[i])
        else:
            past_split, present_split, future_split = fyList[i-1].strip("\n").split("\t"),fyList[i].strip("\n").split("\t"),fyList[i+1].strip("\n").split("\t")
            past_predict, present_predict, future_predict = int(past_split[1]),int(present_split[1]),int(future_split[1])
            if (past_predict == present_predict) or (present_predict == future_predict):
                final.append(fyList[i])
        i+=1
    if i == (length-1):
         past_split, present_split = fyList[i-1].strip("\n").split("\t"),fyList[i].strip("\n").split("\t")
         past_predict, present_predict = int(past_split[1]),int(present_split[1])
         if present_predict == past_predict:
             final.append(fyList[i])
    return(final)
    
def repeatCluster(List):
    final,i = [],0
    while (i+1) < len(List):
        past,present,n = List[i].split()[1],List[i+1].split()[1],1
        while past == present:
            n+=1
            i+=1
            if (i+1) < len(List):
                past,present = List[i].split()[1],List[i+1].split()[1]
            else:
                break
        line_split = List[i].split()
        chrom,locus,orient = line_split[0],line_split[1],line_split[2]
        final.append(locus + "\t" + locus + "\t" + str(n) + "\t" + orient)
        i+=1
    return(final)

def addOrientation(fyList,orientation):
    final = []
    for line in fyList:
        strip = line.strip("\n")
        final.append(strip + "\t" + orientation)
    return(final)
    
def filterPredict_file(fyList):
    final = []
    for line in fyList:
        line_split = line.strip("\n").split("\t")
        if int(line_split[2]) >= 0:
            final.append(line)
    return(final)
    
def filter2(fyList):
    poly,fixed = [],[]
    for line in fyList:
        line_split = line.strip("\n").split("\t")
        if int(line_split[1].split(":")[0]) > 0:
            if int(line_split[1].split(":")[0]) > 25:
                poly.append(line)
            else:
                fixed.append(line)
    return(poly,fixed)

from local_stats import *

def reject_outliers(data, m = 2):
    final = []
    if len(data) > 1:
        std, average = stdev(data), mean(data)
        for d in data:
            if abs(d - average) <= std * m:
                final.append(d)
        return(final)
    else:
        return(data)
        
def cluster(sorted_list, column, distance):
    """The purpose of this function is to take a listed file containing tab delimited lines and cluster the lines into sub-lists based
    on values present in a particular column"""
    first = sorted_list[0].strip("\n").split("\t")
    final,past,present,length,i,clust = [],first[int(column)],first[int(column)],len(sorted_list),0,[]
    while i < length:
        present_line = sorted_list[i].strip("\n").split("\t")
        present = present_line[int(column)]
        diff = int(present) - int(past)
        if diff <= distance:
            clust.append(int(present))
            past = present
            i+=1
        else:
            if clust == []:
                clust.append(int(present))
            clust = reject_outliers(clust,m=2)
            avg,mde = mean(clust),mode(clust)
            if mde == "N/A":
                mde = int(avg)
            final.append(str(mde) + "\t" + str(avg) + "\t" + str(len(clust)))
            clust = []
            clust.append(int(present))
            past = present
            i+=1
    return(final)    

def clustered(sorted_list, column, distance):
    """The purpose of this function is to take a listed file containing tab delimited lines and cluster the lines into sub-lists based
    on values present in a particular column"""
    first = sorted_list[0].strip("\n").split("\t")
    final,past,length,i,clust,support_ID,support = [],first[int(column)],len(sorted_list),1,[first[1]],[str(first[1]) + ":" + str(first[2]) + ";0"],int(first[2])
    while i < length:
        present_line = sorted_list[i].strip("\n").split("\t")
        present = present_line[int(column)]
        diff = int(present) - int(past)
        if diff <= distance:
            clust.append(present_line[1])
            support_ID.append(str(present_line[1]) + ":" + str(present_line[2]) + ";" + str(diff))
            support+=int(present_line[2])
        else:
            if len(clust) > 0:
                final.append(str(past) + "\t" + str(len(clust)) + ":" + str(support) + "\t" + "_".join(support_ID))
            clust = [present_line[1]]
            support_ID = [str(present_line[1]) + ":" + str(present_line[2]) + ";0"]
            support = int(present_line[2])
            past = present
        i+=1
    return(final)    

def joinList(listA,listB):
    final = []
    for a in listA:
        final.append(a)
    for b in listB:
        final.append(b)
    return(final)
    
def predictRefCheck(predict_list,ref_list):
    predict_0,predict_1 = splitBYorient(3,predict_list)
    predict_roo_0, predict_roo_1 = repeatOffenders_only(predict_0),repeatOffenders_only(predict_1)
    ref_0,ref_1 = splitBYorient(3,ref_list)
    cluster_0,cluster_1 = cluster(predict_roo_0,1,15),cluster(predict_roo_1,1,15)   
    add0,add1 = addOrientation(cluster_0, "0"),addOrientation(cluster_1, "1")
    joined = joinList(add0,add1)
    return(sortMe(joined))

def editREPEATMASTER_output(path2rpmskfle):
    rpmsklst = lister(path2rpmskfle)
    final = []
    for line in rpmsklst:
        line_split = line.split()
        chr_num, prime5, prime3, orient, repeat, family, one, three = line_split[4],line_split[5],line_split[6],line_split[8],line_split[9],line_split[10],line_split[11],line_split[13]
        if repeat[0:3] == "Alu":
            if orient == "+":
                final.append(chr_num + "\t" + str(int(prime5) - int(one)) + "\t" + "0" + "\t" + repeat)
            elif orient == "C":
                final.append(chr_num + "\t" + str(int(prime3) + int(three)) + "\t" + "1" + "\t" + repeat)
    return(final)

def zeroGen(numb):
    string = ""
    for i in range(numb):
        string = string + "0"
    return(string)

def zeroThehero(ideal,actual):
    actual = str(actual)
    actual_Len = len(actual)
    diff = int(ideal) - actual_Len
    zeros = zeroGen(diff)
    new_string = str(zeros) + actual
    return(new_string)
    
def fastq_add_seq(fastq_list):
    i, total_length = 0, len(fastq_list)
    while i < total_length:
        line = fastq_list[i]
        if line[0:4] == "@SRR":
            split = line.split(" ")
            seq = fastq_list[i+1]
            print(split[0] + "_" + seq + " " + split[1] + " " + split[2])
        else:
            print(line)
        i+=1      

def reverseComplement(seq):
    seq_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return("".join([seq_dict[base] for base in reversed(seq)]))
    
def splitFastq(fastqList):
    final,l,total_len,name,zeroIdeal = [],0,len(fastqList),1,len(str(int(len(fastqList)/4)))
    while l < (total_len - 3):
        seq,qual,length = fastqList[l+1],fastqList[l+3],len(fastqList[l+3])
        name_a, name_b, seq_a, seq_b, qual_a, qual_b = zeroThehero(zeroIdeal,name),zeroThehero(zeroIdeal,(name+1)), seq[:(int(length/2))], seq[(int(length/2)):], qual[:(int(length/2))], qual[(int(length/2)):]
        final.append("@" + name_a + "/1")
        final.append(seq_a)
        final.append("+")
        final.append(qual_a)  
        final.append("@" + name_a + "/2")
        final.append(seq_b)
        final.append("+")
        final.append(qual_b) 
        name+=1
        l+=4
    return(final)
    
def splitFastqW_ID(fastqList):
    final,l,total_len = [],0,len(fastqList)
    while l < (total_len - 3):
        name,seq,qual,length = fastqList[l].split()[0],fastqList[l+1],fastqList[l+3],len(fastqList[l+3])
        seq_a, seq_b, qual_a, qual_b = seq[:(int(length/2))], reverseComplement(seq[(int(length/2)):]), qual[:(int(length/2))], ((qual[(int(length/2)):])[::-1])
        final.append(name + "/1")
        final.append(seq_a)
        final.append("+")
        final.append(qual_a)  
        final.append(name + "/2")
        final.append(seq_b)
        final.append("+")
        final.append(qual_b) 
        l+=4
    return(final)
    
def splitFastqFromBAM(fastqList):
    final,l,total_len,name,zeroIdeal = [],0,len(fastqList),1,len(str(int(len(fastqList)/4)))
    while l < (total_len - 3):
        seq,qual,length = fastqList[l+1],fastqList[l+3],len(fastqList[l+3])
        name_a, name_b, seq_a, seq_b, qual_a, qual_b = zeroThehero(zeroIdeal,name),zeroThehero(zeroIdeal,(name+1)), seq[:(int(length/2))], reverseComplement(seq[(int(length/2)):]), qual[:(int(length/2))], ((qual[(int(length/2)):])[::-1])
        final.append("@" + name_a)
        final.append(seq)
        final.append("+")
        final.append(qual) 
        name+=1
        l+=4
    return(final)
    
def appendSeq2mapped(mappedList,fastqList):
    mappedID_list,newID = [],[]
    for m in mappedList:
        m_split = m.strip("\n").split("\t")
        mappedID_list.append(m_split[0])
        newID.append(m_split[0] + "_" + m_split[1] + "_" + m_split[2] + "_" + m_split[3] + "_" + m_split[4])
    final,n = [],0
    for m in mappedID_list:
        i = (int(m) - 1) * 4
        final.append("@" + newID[n])
        final.append(fastqList[i + 1])
        final.append("+")
        final.append(fastqList[i+3])
        n+=1
    return(final)

def disk_usage(path):
    """Return disk usage statistics about the given path.
    
    Returned valus is a named tuple with attributes 'total', 'used' and
    'free', which are the amount of total, used and free space, in bytes.
    """
    st = os.statvfs(path)
    free = st.f_bavail * st.f_frsize
    total = st.f_blocks * st.f_frsize
    used = (st.f_blocks - st.f_bfree) * st.f_frsize
    return(total, used, free)
    
def checkAvailable(path2file):
    path2fastq = path2file
    total,used,free = disk_usage('/') #Check the total, used, and free disk space
    total,used,free = total/1000000000,used/1000000000,free/1000000000 #Convert units from B to GB
    fastq_size = os.path.getsize(path2fastq)/1000000000
    
    if free <= fastq_size + 10.0:
        print("Error:\n Available disk space in current directory is insufficient. Total required space is " + str(fastq_size + 10.0) + " GB. \nTotal available space is " + str(free) + " GB.", file=sys.stderr)
        raise SystemExit

def find_nth(haystack, needle, n):
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+len(needle))
        n -= 1
    return(start)
    
def sendSMS(message):
    client = TwilioRestClient(account= 'AC96970e768653b33b68267efa24150b54',token='335c3747266f29e53b9c25e98feef9a7')
    client.messages.create(from_='(872) 395-6034',to='(773) 517-3665',body=message)
    
def fastqSplit(path2fastq,fastqDir,orgID,srr_number):
    subprocess.call("mkdir " + fastqDir + "tmp",shell=True)
    subprocess.call("mkdir " + "tmp_" + str(orgID),shell=True)
    subprocess.call("split -l 800000 " + path2fastq + " tmp/split" + srr_number,shell=True)
    sendSMS(orgID + ":\n splitA COMPLETE")
    #subprocess.call("rm -f " + path2fastq, shell=True)
    splitfiles = [f for f in os.listdir(fastqDir + "tmp/") if os.path.isfile(os.path.join(fastqDir + "tmp/", f))]
    for file in splitfiles:
        filepath = "tmp/" + file
        reformatted = splitFastq(lister(filepath))
        with open(fastqDir + "tmp_" + str(orgID) + "/" + file,"w") as out:
            for line in reformatted:
                out.write(line + "\n")
        subprocess.call("rm -f " + filepath, shell=True)
    subprocess.call("rmdir " + fastqDir + "tmp", shell=True)
    sendSMS(orgID + ":\nsplitB COMPLETE")
    
def mappedHTsplit(aligned2ALU,fastq_copy):
    final = []
    for line in aligned2ALU:
        line_split = line.split()
        if len(line_split) > 3:
            if (line_split[2][0] == "A") and (int(line_split[4]) >= 30):
                QNAME,FLAG,RNAME,POS,MAPQ,CIGAR = line_split[0],line_split[1],line_split[2],line_split[3],line_split[4],line_split[5]
                temp = ["@" + QNAME + "_" + RNAME[-1] + "_" + FLAG + "_" + POS + "_" + MAPQ + "_" + CIGAR]
                if int(QNAME) != 1:
                    for i in range(1,4):
                        temp.append(fastq_copy[int(((int(QNAME))-1)*4+i)])
                else:
                    for i in range(1,4):
                        temp.append(fastq_copy[int(i)])
                if RNAME[-1] == "H":
                    final.extend(temp)
                elif RNAME[-1] == "T":
                    final.extend(temp)
    return(final)
    
def extractMapped(sam,fastq_copy):
    final,i = [],0
    while i < len(sam):
        line_split = sam[i].split()
        if (len(line_split) > 5):
            next_split = sam[i+1].split()
            if ((str(line_split[2]) + str(line_split[5])) != "**"):
                QNAME,FLAG,POS,MAPQ,CIGAR,PNEXT = line_split[0],line_split[1],line_split[3],line_split[4],line_split[5],line_split[7]
                QNEXT,FNXT,NPS,NMPQ,NXTCR,PLAST = next_split[0],next_split[1],next_split[3],next_split[4],next_split[5],next_split[7]           
                temp = []
                t_string = "@" + QNAME  + ";" + FLAG + ";" + POS + ";" + MAPQ + ";" + CIGAR + ";" + PNEXT + "{" + QNEXT  + ";" + FNXT + ";" + NPS + ";" + NMPQ + ";" + NXTCR + ";" + PLAST
                temp.append(t_string + ";1")
                if int(QNAME) != 1:
                    for n in range(1,4):
                        temp.append(fastq_copy[int(((int(QNAME))-1)*8+n)])
                    temp.append(t_string + ";2")
                    for n in range(5,8):
                        temp.append(fastq_copy[int(((int(QNAME))-1)*8+n)])
                else:
                    for n in range(1,4):
                        temp.append(fastq_copy[int(n)])
                    temp.append(t_string + ";2")
                    for n in range(5,8):
                        temp.append(fastq_copy[int(n)])
                for t in temp:
                    final.append(t)
            i+=2
        else:
            i+=1
            
    return(final)
    
    
def masterPIPE(org_ID,path2split,path2ALU,path2REF):
    start = timer()
    splitfiles = [f for f in os.listdir(path2split) if os.path.isfile(os.path.join(path2split, f))]
    split_len,i = len(splitfiles),1
    print("AluPredict\n" + str(split_len) + " runs to process", file=sys.stderr)
    for fastq in splitfiles:
        start_iter = timer()
        print("Processing " + str(i) + " of " + str(split_len), file=sys.stderr)
        fastq_copy = lister(path2split + fastq)
        print("Aligning reads to Alu reference...", file=sys.stderr)
        aligned2ALU = Popen(['bwa', 'mem', '-t', '8', '-k', '18', '-p', '-B', '2', '-E', '1', '-O', '2', path2ALU, path2split + fastq], stdout=PIPE)
        #kinkos(mappedHTsplit(aligned2ALU.communicate()[0].decode("utf-8").split("\n"),fastq_copy))
        kinkos(extractMapped(aligned2ALU.communicate()[0].decode("utf-8").split("\n"),fastq_copy))
        i+=1
    
        
def orientConverter(FLAG):
    if FLAG != "0":
        binary = bin(int(FLAG))[2:]
        return(str(binary[-5]))
    else:
        return("0")

def predictCandidate(samSplit):
    aluSplit = samSplit[0].split("_")
    ANAM,ATAG,ACIG,APOS,AORN = aluSplit[1],aluSplit[2],aluSplit[5],aluSplit[3],orientConverter(aluSplit[2])
    RNAM,RTAG,RCIG,RPOS,RORN = samSplit[2][3:],samSplit[1],samSplit[5],samSplit[3],orientConverter(samSplit[1])
    combOrient = AORN+RORN
    #if AluHT == 'H':
    
def sam2fasta(sam):
    if len(sam.split('\n')) == 1:
        return('>'+sam.split()[0]+'\n'+sam.split()[9])
    else:
        fasta = []
        for line in sam.split('\n'):
            if len(line.split()) > 3:
                fasta.append('>'+line.split()[0]+'\n'+line.split()[9])
        return(fasta)
        
def extractRefREPEAT(name,path2repeatmasker_library):
    final = []
    rptLIST,i = lister(path2repeatmasker_library),0
    while i < (len(rptLIST)):
        if name in rptLIST[i]:
            final.append(rptLIST[i])
            i+=1
            while (i < len(rptLIST)) and (rptLIST[i][0] != '>') :
                final.append(rptLIST[i])
                i+=1
        else:
            i+=1
    return(final)
    

def filterPredict(samList):       #00:000000000:0:0
    final = []                            #CR:APPROXPOS:H:O
    for line in samList:
        if (len(line.split()) > 5) and ((str(line.split()[2]) + str(line.split()[5])) != "**"):
            final.append(line)
    return(final)
    
def filterPredict2(samList):       #00:000000000:0:0
    i,final = 0,[]                       #CR:APPROXPOS:H:O
    while i < len(samList):
        temp,support,n = [[],[]],[],1
        temp[int(samList[i].split()[0][-1:])-1].append(samList[i])
        support.append(int(samList[i].split()[4])) 
        tag = True
        while (((i+n) < len(samList))) and (tag == True):
            if ((samList[i].split()[0].split(";")[0]) == (samList[i+n].split()[0].split(";")[0])):
                temp[int(samList[i+n].split()[0][-1:])-1].append(samList[i+n])
                support.append(int(samList[i+n].split()[4]))
                i+=1
            else:
                x = 0
                for l in support:
                    if l >= 30:
                        x+=1
                if x != 0:
                    final.extend(temp[0])
                    final.extend(temp[1])
                tag = False 
        i+=1
    return(final)

def set_range(start,stop):
    final=set()
    for i in range(start,stop+1):
        final.add(i)
    return(final)

def generate_PairedEnd(sequence,gap,read_length,depth):
    y,n=0,0
    tick=int(1/(depth/(2*read_length)))
    with open('forward.fq','w') as forward_fq:
        with open('reverse.fq','w') as reverse_fq:
            with open('forward.infoTemp','w') as forward_info:
                with open('reverse.infoTemp','w') as reverse_info:
                    while (y+(2*read_length)+gap) < len(sequence):
                        locusF_start,locusF_stop=y,y+read_length
                        locusR_start,locusR_stop=gap+y+read_length,gap+y+read_length+read_length
                        seqF,seqR=sequence[locusF_start:locusF_stop],sequence[locusR_start:locusR_stop]
                        seq=Seq(seqR, IUPAC.ambiguous_dna)
                        name=str(zeroThehero(12,n+1))
                        forward_fq.write('@'+name+'\tCO:Z:F'+str(y+1)+'\n'+seqF+'\n+\n'+quality+'\n')
                        reverse_fq.write('@'+name+'\tCO:Z:R'+str(gap+y+read_length+1)+'\n'+str(seq.reverse_complement())+'\n+\n'+quality+'\n')
                        forward_info.write(name+'F\t'+str(y+1)+'-'+str(y+read_length)+'\n')
                        reverse_info.write(name+'R\t'+str(gap+y+read_length+1)+'-'+str(gap+y+read_length+read_length)+'\n')
                        if tick>=1:
                            y+=tick
                        n+=1

def rptStEnd(rptLine):# Position inside of repeat str(start,end)
    orient=rptLine.split()[8]
    if orient=='+':
        start,end=rptLine.split()[11],rptLine.split()[12]
    elif orient=='C':
        start,end=rptLine.split()[13],rptLine.split()[12]
    return(start,end)

def rrSpanner(orient,rptS,rptE,seqS,seqE,posRMstart,posRMend,que):
    if que==1:
        rdSpan=str(rptS)+'-'+str(seqE)
        if orient=='+':
            rpSpan=posRMstart+'-'+str(seqE-rptS+int(posRMstart))
        elif orient=='C':
            rpSpan=str(int(posRMend)-(seqE-rptS))+'-'+posRMend
    elif que==2:
        rdSpan=str(rptS)+'-'+str(rptE)
        rpSpan=posRMstart+'-'+posRMend
    elif que==3:
        rdSpan=str(seqS)+'-'+str(seqE)
        if orient=='+':
            rpSpan=str(seqS-rptS+int(posRMstart))+'-'+str(int(posRMend)-(rptE-seqE))
        elif orient=='C':
            rpSpan=str(rptE-seqE+int(posRMstart))+'-'+str(int(posRMend)-(seqS-rptS))
    elif que==4:
        rdSpan=str(seqS)+'-'+str(rptE)
        if orient=='+':
            rpSpan=str(seqS-rptS+int(posRMstart))+'-'+posRMend
        elif orient=='C':
            rpSpan=posRMstart+'-'+str(int(posRMend)-(seqS-rptS))
    if abs((rptE-rptS)-((int(posRMend)-int(posRMstart))))>20:
        rpSpan='N/A'
    return(rdSpan,rpSpan)
                            
def addMEinfo(rmList,infoList,tag,rptCount):
    with open(tag+'.info','w') as outfile:
        i,r,iLen,rLen=0,0,len(infoList),len(rmList)
        while i<iLen:
            if r<rLen:
                rptS,rptE,seqS,seqE=int(rmList[r].split()[5]),int(rmList[r].split()[6]),int(infoList[i].split()[1].split('-')[0]),int(infoList[i].split()[1].split('-')[1])
                locus,skip=str(seqS)+'-'+str(seqE),False
                readID,meiInfo,readSpan,rptSpan,rptN=infoList[i].split()[0],'N/A\tN/A','N/A','N/A','N\A'                   
                if (seqS<=rptS):
                    if (seqE>rptS):
                        meiInfo,orient=rmList[r].split()[10]+'\t'+rmList[r].split()[9],rmList[r].split()[8]
                        posRMstart,posRMend=rptStEnd(rmList[r])
                        rptN=str(r)
                        if (seqE<=rptE):
                            readSpan,rptSpan=rrSpanner(orient,rptS,rptE,seqS,seqE,posRMstart,posRMend,1)
                        elif (seqE>rptE):
                            readSpan,rptSpan=rrSpanner(orient,rptS,rptE,seqS,seqE,posRMstart,posRMend,2)
                            if seqE>int(rmList[r+1].split()[5]):
                                classR,subR,readSpanL,rptSpanL,r2,rptNL=[rmList[r].split()[10]],[rmList[r].split()[9]],[readSpan],[rptSpan],r+1,[str(r)]
                                rptS2,rptE2=int(rmList[r2].split()[5]),int(rmList[r2].split()[6])
                                posRMstart2,posRMend2=rptStEnd(rmList[r2])
                                while (seqE>rptS2):
                                    classR.append(rmList[r2].split()[10])
                                    subR.append(rmList[r2].split()[9])
                                    if (seqE>rptE2):
                                        rdSpan2,rpSpan2=rrSpanner(orient,rptS2,rptE2,seqS,seqE,posRMstart2,posRMend2,2)
                                    else:
                                        rdSpan2,rpSpan2=rrSpanner(orient,rptS2,rptE2,seqS,seqE,posRMstart2,posRMend2,1)
                                    readSpanL.append(rdSpan2)
                                    rptSpanL.append(rpSpan2)
                                    rptNL.append(str(r2))
                                    r2+=1
                                    if r2<rLen:
                                        rptS2,rptE2=int(rmList[r2].split()[5]),int(rmList[r2].split()[6])
                                    else:
                                        break
                                readSpan,rptSpan,meiInfo,rptN=','.join(readSpanL),','.join(rptSpanL),','.join(classR)+'\t'+','.join(subR),','.join(rptNL)
                elif (seqS>rptS):
                    if (seqS<rptE):
                        meiInfo,orient=rmList[r].split()[10]+'\t'+rmList[r].split()[9],rmList[r].split()[8]
                        posRMstart,posRMend=rptStEnd(rmList[r])
                        rptN=str(r)
                        if seqE<=rptE:
                            readSpan,rptSpan=rrSpanner(orient,rptS,rptE,seqS,seqE,posRMstart,posRMend,3)
                        elif seqE>rptE:
                            readSpan,rptSpan=rrSpanner(orient,rptS,rptE,seqS,seqE,posRMstart,posRMend,4)
                            if seqE>int(rmList[r+1].split()[5]):
                                classR,subR,readSpanL,rptSpanL,r2,rptNL=[rmList[r].split()[10]],[rmList[r].split()[9]],[readSpan],[rptSpan],r+1,[str(r)]                                                       
                                rptS2,rptE2=int(rmList[r2].split()[5]),int(rmList[r2].split()[6])
                                posRMstart2,posRMend2=rptStEnd(rmList[r2])
                                while (seqE>rptS2):
                                    classR.append(rmList[r2].split()[10])
                                    subR.append(rmList[r2].split()[9])
                                    if (seqE>rptE2):
                                        rdSpan2,rpSpan2=rrSpanner(orient,rptS2,rptE2,seqS,seqE,posRMstart2,posRMend2,4)
                                    else:
                                        rdSpan2,rpSpan2=rrSpanner(orient,rptS2,rptE2,seqS,seqE,posRMstart2,posRMend2,3)
                                    readSpanL.append(rdSpan2)
                                    rptSpanL.append(rpSpan2)
                                    rptNL.append(str(r2))
                                    r2+=1
                                    if r2<rLen:
                                            rptS2,rptE2=int(rmList[r2].split()[5]),int(rmList[r2].split()[6])
                                    else:
                                        break
                                readSpan,rptSpan,meiInfo,rptN=','.join(readSpanL),','.join(rptSpanL),','.join(classR)+'\t'+','.join(subR),','.join(rptNL)
                if (seqS>=rptE):
                    skip=True
                    r+=1
                if skip==False:
                    outfile.write(readID+'\t'+locus+'\t'+meiInfo+'\t'+readSpan+'\t'+rptSpan+'\t'+rptN+'\n')
                    i+=1
                    if readSpan!='N/A':
                        nList=rptN.split(',')
                        if len(nList)==1:
                            if len(rptCount)>=r+1:
                                rptCount[r]+=1
                            else:
                                rptCount.append(1)
                        else:
                            for item in nList:
                                if len(rptCount)>=int(item)+1:
                                    rptCount[int(item)]+=1
                                else:
                                    rptCount.append(1)
        print(str(r)+'\t'+str(i),file=sys.stderr)
    return(rptCount)

def aluSupport(rptCount,repeatMaskerfile):
    with open('aluSupportFL.txt','w') as fullLength:
        with open('aluSupportTR.txt','w') as truncated:
            with open('noSupport.txt','w') as noSupport:
                r,fl,tr,na=0,[],[],0
                while r<len(rptCount):
                    rptSplit=repeatMaskerfile[r].split()
                    if rptSplit[9][:3]=='Alu':
                        if rptCount[r]>0:
                            start,end=rptStEnd(repeatMaskerfile[r]) 
                            if (int(start)<4) and (int(end)>266):
                                fullLength.write(rptSplit[9]+'\t'+start+'\t'+end+'\t'+str(rptCount[r])+'\n')
                                fl.append(rptCount[r])
                            else:
                                truncated.write(rptSplit[9]+'\t'+start+'\t'+end+'\t'+str(rptCount[r])+'\n')
                                tr.append(rptCount[r])
                        else:
                            noSupport.write(rptSplit[9]+'\t'+start+'\t'+end+'\t'+str(rptCount[r])+'\n')
                            na+=1
                    r+=1
    print('\nFull Length:\ntotal = '+str(len(fl))+'\naverage = '+str(np.mean(fl))+'\n',file=sys.stderr)
    print('Truncated:\ntotal = '+str(len(tr))+'\naverage = '+str(np.mean(tr))+'\n',file=sys.stderr)
    print('No Support:\ntotal = '+str((na)),file=sys.stderr)
    
def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False
    
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
