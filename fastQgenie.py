"""
Created on Fri Apr 15 15:52:36 2016

@author: vallmer
"""

import sys
import re
from sam import *
from essentials import *
from sra import *
import linecache

def fastqDump(accession,outfile): #NCBI SRA accession number
    subprocess.call('fastq-dump --split-spot -Z ' + accession + '> ' + outfile , shell=True)  
    
def nesoni(sra):
    subprocess.call('nesoni clip: clipped --gzip no interleaved: ' + 'out' , shell=True)
    subprocess.call('rm -f out', shell=True)
    subprocess.call('rm -f clipped_log.txt', shell=True)

def bwaRef(path2ref,paired,outfile):
    if paired == 'yes':
        subprocess.call('bwa mem -p -t 8 -C -T 20 ' + path2ref + ' out2 > ' + outfile , shell=True)
    elif paired == 'no':
        subprocess.call('bwa mem -t 8 -C -T 20 ' + path2ref + ' alu_s_mapped.fastq > ' + outfile , shell=True)
    
def fixmate(path2sam,output):
    subprocess.call('samtools fixmate -O sam ' + path2sam + ' ' + output,shell=True)
    #subprocess.call('rm -f '+path2sam, shell=True)

def samSort(output,path2bam):
    subprocess.call('samtools sort -O sam -o ' + output + ' -T /tmp/lane_temp ' + path2bam,shell=True)
    #subprocess.call('rm -f '+path2bam, shell=True)
    
def mappedORmappedMate(inSamp,inSams):
    with open('alu_p_mapped.fastq','w') as out:
        with open(inSamp,'r') as inpute:
            for line in inpute:
                splyt = line.strip('\n').split()
                if len(splyt) > 10:
                    flag = splyt[1]
                    if flag != "0":
                        binary = bin(int(flag))[2:]
                        if (str(binary[-3]) == "0") or (str(binary[-4]) == "0"):
                            out.write('@'+splyt[0]+"\n"+splyt[9]+"\n+\n"+splyt[10]+"\n")
                            #out.write(line+'\n')
    with open(inSams,'r') as inputes:
        with open('alu_s_mapped.fastq','w') as outs:
            for line in inputes:
                splyt=line.strip('\n').split()
                if len(splyt) > 10:
                    outs.write('@'+splyt[0]+"\tOP:i:1\n"+splyt[9]+"\n+\n"+splyt[10]+"\n")     

def mappedORmappedMate_SAM(inSamp):
    with open('alu_p_mapped.sam','w') as out:
        with open(inSamp,'r') as inpute:
            for line in inpute:
                splyt = line.strip('\n').split()
                if len(splyt) > 10:
                    flag = splyt[1]
                    if flag != "0":
                        binary = bin(int(flag))[2:]
                        if (str(binary[-3]) == "0") or (str(binary[-4]) == "0"):
                            out.write(line)

def bwaAlu(path2aluConsensus):
    subprocess.call('bwa mem -C -p -t 8 -k 10 -T 15 ' + path2aluConsensus + ' clipped_paired.fq > alu_p.sam', shell=True)#############################
    subprocess.call('bwa mem -t 8 -k 10 -T 15 ' + path2aluConsensus + ' clipped_single.fq > alu_s.sam', shell=True)#############################
    #subprocess.call('samtools view -F 4 alu_p.sam > alu_p_mapped.sam', shell=True)
    subprocess.call('samtools view -F 4 alu_s.sam > alu_s_mapped.sam', shell=True)

def processMapped(inSam,alu_or_ref,PorS):
    n = 0
    with open(alu_or_ref+'.'+PorS,'w') as out:
        with open(inSam,'r') as inpute:
            for line in inpute:
                splyt = line.strip('\n').split()
                if len(splyt) > 10:
                    if (splyt[0][:3] == 'SRR') and (int(splyt[4]) >= 25) and('OP:i:' in line):
                        out.write(splyt[0].split('.')[1] + '.' + str(line.split('OP:i:')[1][0]) + str(99-int(splyt[4])) + '\t' + splyt[1] + '\t' + splyt[2] + '\t' + splyt[3] + '\t' + splyt[5] + '\t' + splyt[4] + '\n')
                        n+=1
            subprocess.call('sort -nk1 --parallel=8 '+alu_or_ref+'.'+PorS + ' >' + alu_or_ref+'.sort.'+PorS,shell=True)
            #subprocess.call('rm -f '+alu_or_ref+'.'+PorS, shell=True)
    with open(alu_or_ref+'.uniq.sort.'+PorS,'w') as uniq:
        fp,u = open(alu_or_ref+'.sort.'+PorS),0
        for x, line in enumerate(fp):
            if x==0:
                past=line
            elif x==1:
                present=line
            elif x < (n-1):
                future=line
                past_split,present_split,future_split = past.strip('\n').split('\t'),present.strip('\n').split('\t'),future.strip('\n').split('\t')
                if (past_split[0]!=present_split[0])and(present_split[0]!=future_split[0]):
                    uniq.write(present)
                    past=present
                    present=future
                    u+=1
            else:
                future=line
                past_split,present_split,future_split = past.strip('\n').split('\t'),present.strip('\n').split('\t'),future.strip('\n').split('\t')
                if present_split[0]!=future_split[0]:
                    uniq.write(future)
                    u+=1
                    if present_split[0]!=past_split[0]:
                        uniq.write(present)
                        u+=1
    #subprocess.call('rm -f '+alu_or_ref+'.sort.'+PorS, shell=True)
    return(u)

def processMapped_ALU_S(inSam):
    n = 0
    with open('alu.s','w') as out:
        with open(inSam,'r') as inpute:
            for line in inpute:
                splyt = line.strip('\n').split()
                if len(splyt) > 10:
                    if (splyt[0][:3] == 'SRR') and (int(splyt[4]) >= 25):
                        out.write(splyt[0].split('.')[1] + '.1' + str(99-int(splyt[4])) + '\t' + splyt[1] + '\t' + splyt[2] + '\t' + splyt[3] + '\t' + splyt[5] + '\t' + splyt[4] + '\n')
                        n+=1
            subprocess.call('sort -nk1 --parallel=8 alu.s > alu.sort.s',shell=True)
            #subprocess.call('rm -f '+alu_or_ref+'.'+PorS, shell=True)
    with open('alu.uniq.sort.s','w') as uniq:
        fp,u = open('alu.sort.s'),0
        for x, line in enumerate(fp):
            if x==0:
                past=line
            elif x==1:
                present=line
            elif x < (n-1):
                future=line
                past_split,present_split,future_split = past.strip('\n').split('\t'),present.strip('\n').split('\t'),future.strip('\n').split('\t')
                if (past_split[0]!=present_split[0])and(present_split[0]!=future_split[0]):
                    uniq.write(present)
                    past=present
                    present=future
                    u+=1
            else:
                future=line
                past_split,present_split,future_split = past.strip('\n').split('\t'),present.strip('\n').split('\t'),future.strip('\n').split('\t')
                if present_split[0]!=future_split[0]:
                    uniq.write(future)
                    u+=1
                    if present_split[0]!=past_split[0]:
                        uniq.write(present)
                        u+=1
    #subprocess.call('rm -f '+alu_or_ref+'.sort.'+PorS, shell=True)
    return(u)
    
def processFastq(inFastq,processedFq):
    with open(inFastq,'r') as inFile:
        with open(processedFq,'w') as fq:
            for line in inFile:
                if len(line.split()) > 2:
                    if (line.split()[2][0:6]) == 'length':
                        fq.write(line.split()[0] + '\n')
                else:
                    fq.write(line)
    subprocess.call('rm -f '+ inFastq, shell=True)
    with open('out','a') as out:
        fp = open(processedFq)
        for x, line in enumerate(fp):
            if ((x + 8)/8) == (int((x + 8)/8)):
                out.write(line.strip('\n') +  " OP:i:1\n")
            elif ((x + 7)/8) == (int((x + 7)/8)):
                out.write(line)
            elif ((x + 6)/8) == (int((x + 6)/8)):
                out.write(line.strip('\n') + "\n")
            elif ((x + 5)/8) == (int((x + 5)/8)):
                out.write(line)
            elif ((x + 4)/8) == (int((x + 4)/8)):
                out.write(line.strip('\n') + " OP:i:2\n")
            elif ((x + 3)/8) == (int((x + 3)/8)):
                out.write(line)
            elif ((x + 2)/8) == (int((x + 2)/8)):
                out.write(line.strip('\n') + "\n")
            elif ((x + 1)/8) == (int((x + 1)/8)):
                out.write(line)
        fp.close()
    subprocess.call('rm -f '+ processedFq, shell=True)

def processFastq2(processedFq):
    with open('out2','a') as out:
        fp = open(processedFq)
        for x, line in enumerate(fp):
            if ((x + 8)/8) == (int((x + 8)/8)):
                out.write(line.strip('\n') +  " OP:i:1\n")
            elif ((x + 7)/8) == (int((x + 7)/8)):
                out.write(line)
            elif ((x + 6)/8) == (int((x + 6)/8)):
                out.write(line.strip('\n') + "\n")
            elif ((x + 5)/8) == (int((x + 5)/8)):
                out.write(line)
            elif ((x + 4)/8) == (int((x + 4)/8)):
                out.write(line.strip('\n') + " OP:i:2\n")
            elif ((x + 3)/8) == (int((x + 3)/8)):
                out.write(line)
            elif ((x + 2)/8) == (int((x + 2)/8)):
                out.write(line.strip('\n') + "\n")
            elif ((x + 1)/8) == (int((x + 1)/8)):
                out.write(line)
        fp.close()
    #subprocess.call('rm -f '+ processedFq, shell=True)
        
def spanHead(alu_cig,aluStart):
    previous,charList = alu_cig.split('M')[-1],['I','D','N','S','H','P','X']
    index = [0]
    for i in range(len(previous)):
        if previous[i] in charList:
            index.append(i)
    index.append(len(previous)-1)
    if len(index) == 2:
        return('1')
    else:
        distance,d1,d2 = 0,0,1
        for i in range(len(index)-2):
            distance+=int(previous[index[d1+1]:index[d2]])
            d1+=1
            d2+=1
        if (int(aluStart)-distance) > -20:
            return('0')
        else:
            return('1')
            
def predictGenie(alu_split,ref_split):
    alu_Orient,ref_Orient,splitHead = orientConverter(alu_split[1]),orientConverter(ref_split[1]),spanHead(alu_split[4],alu_split[3])
    if int(alu_split[3]) < 60:
        alu_pos = 'H'
    elif int(alu_split[3]) > 200:
        alu_pos = 'T'
    else:
        alu_pos = 'M'
    if ref_split[4] != '*':
        mapped = '0'
    else:
        mapped = '1'
    if (alu_Orient == '0') and (ref_Orient == '0'):
        if ref_split[4] != '*':
            chrom,position,orient = ref_split[2][3:],zeroThehero(9,(int(ref_split[3]) - int(alu_split[3]))),'0'
        else:
            chrom,position = ref_split[2][3:],zeroThehero(9,(int(ref_split[3]) - int(alu_split[3])))

def aluPredicter(org,sra_acc_num,path2Alu,path2Ref,PorS):
    subprocess.call('cat '+ path2Alu + ' ' + path2Ref + ' > uniq.cat.' + PorS,shell=True)
    subprocess.call('sort -nk1 --parallel=8 uniq.cat.' + PorS + ' > uniq.sort.' + PorS,shell=True)
    fp = open('uniq.sort.'+PorS)
    with open(org+'_'+sra_acc_num+'_Alupredictions_'+PorS+'.txt','a') as final:
        for x, line in enumerate(fp):
            if x==0:
                past=line
            else:
                present=line
                past_split,present_split = past.strip('\n').split('\t'),present.strip('\n').split('\t')
                if (past_split[0][:-2]==present_split[0][:-2]):
                    final.write(past.strip('\n')+';'+present.strip('\n')+'\n')
                past=present
    
def removeDuplicates(inpath,outpath):
    subprocess.call('java -Xmx8G -jar /home/vallmer/programs/picard.jar MarkDuplicates INPUT=' + inpath + ' OUTPUT=' + outpath + ' METRICS_FILE=metrics.txt REMOVE_DUPLICATES=true READ_NAME_REGEX=null', shell=True)
    #subprocess.call('rm -f '+inpath, shell=True)
    

def pipeline(org,sra_acc_num,path2ref,path2aluConsensus):
    sra_list = lister(sra_acc_num)
    for sra in sra_list:
        fastqDump(sra,sra+'.fq')
        processFastq(sra+'.fq',sra+'.processed.fq')
        subprocess.call('rm -r /home/vallmer/ncbi/public/sra/*', shell=True)
    nesoni('out')
    
pipeline('30388','sra','/home/vallmer/programs/genomes/rheMac8.fa','/home/vallmer/programs/genomes/aluConsensus.fa')

#processFastq(sys.argv[1],sys.argv[2])
