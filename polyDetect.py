import sys
import re
from sam import *
from essentials import *
from sra import *
import linecache
import os
from subprocess import Popen, PIPE
from timeit import default_timer as timer

path2aluConsensus,path2referenceGenome,path2_PE_fq,path2_SE_fq,path2outdir,orgID,chromLen,path2polyA=sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8]

def editFQ_title(inpath,outpath):
    with open(inpath,'r') as infastq:
        with open(outpath,'w') as outfastq:
            for line in infastq:
                if line[:4]=='@SRR':
                    outfastq.write('@'+str(zeroThehero(12,(line.split()[0].split('.')[1])))+'\t'+line.split()[1]+'\n')
                else:
                    outfastq.write(line)

def bwaMeM(path2ref,path2fsq,P_S,S_L,outpath):
    """The Purpose of This is to Excute Bwa Mem"""
    with open(outpath,'w') as out:
        if S_L=='S': #This means that we are using the BWA standard parameters
            if P_S=='P':
                subprocess.call(['bwa', 'mem', '-t', '8', '-p','-C', path2ref, path2fsq],stdout=out)
            elif P_S=='S':
                subprocess.call(['bwa', 'mem', '-t', '8', '-C', path2ref, path2fsq],stdout=out)
        elif S_L=='L': #This means that we are using the BWA more liberal parameters
            if P_S=='P':
                subprocess.call(['bwa', 'mem', '-t', '8', '-p','-C', '-k', '10', '-T', '15', path2ref, path2fsq],stdout=out)
            elif P_S=='S':
                subprocess.call(['bwa', 'mem', '-t', '8', '-C', '-k', '10', '-T', '15', path2ref, path2fsq],stdout=out)

def bowtie2(path2ref,path2fsq,P_S,outpath):
    """The Purpose of This is to Excute bowtie2"""
    with open(outpath,'w') as out:
        if P_S=='S': #This means that we are using the BWA standard parameters
            subprocess.call(['bowtie2', '-p', '8', '--very-sensitive', '-x', path2ref, '-U', path2fsq, '-S', outpath],stdout=out)
        elif P_S=='P': #This means that we are using the BWA more liberal parameters
            subprocess.call(['bowtie2', '-p', '8', '--very-sensitive', '-x', path2ref, '--interleaved', path2fsq, '-S', outpath],stdout=out)
            
def mappOnly(in_sam,out_sam):
    with open(in_sam,'r') as inSam:
        with open(out_sam,'w') as outSam:
            for line in inSam:
                if line[0]!='@':
                    if line.split()[2]!='*':
                        outSam.write(line)

def aluClipped(samline):
    s,tooShort=sam(samline),False
    if (s.cigar!='*') and (int(s.pos)<=270):
        cigar_span,cigar_char=splitCigar(s.cigar)
        position='M'
        if (len(cigar_span)>1) and ('H' not in s.cigar):
            mFirst,mLast=cigar_char.index('M'),find_last(cigar_char,'M')
            i,m=mFirst,0
            while i<=mLast:
                m+=int(cigar_span[i])
                i+=1
            if int(s.pos)<=80:
                m+=(int(s.pos)-1)
                position='H'
            elif (int(s.pos)+m)>=267:
                if (int(s.pos)+m)>282:
                    position='P'
                else:
                    position='T'
            else:
                m=len(s.seq)
                position='M'                
            if mFirst==0:
                if (len(s.seq)-m)>=20:
                    seq,qual=s.seq[m:],s.qual[m:]
                else:
                    tooShort=True
            elif mLast==(len(cigar_char)-1):
                if (len(s.seq)-m)>=20:
                    seq,qual=s.seq[:(len(s.seq)-m)],s.qual[:(len(s.seq)-m)]
                else:
                    tooShort=True
            else:
                begin,end,b,e=0,0,0,(mLast+1)
                while b<mFirst:
                    begin+=int(cigar_span[b])
                    b+=1
                while e<len(cigar_char):
                    end+=int(cigar_span[e])
                    e+=1
                if (len(s.seq)-m)>=20:
                    if b>e:
                        if b>=20:
                            seq,qual=s.seq[:b],s.qual[:b]
                        else:
                            tooShort=True
                    elif e>b:
                        if e>=20:
                            seq,qual=s.seq[:e],s.qual[:e]
                        else:
                            tooShort=True
                    else:
                        tooShort=True
                            
                else:
                    tooShort=True
        elif (len(cigar_span)==1):
            mFirst,mLast=cigar_char.index('M'),find_last(cigar_char,'M')
            i,m=mFirst,0
            while i<=mLast:
                m+=int(cigar_span[i])
                i+=1
            if int(s.pos)<=80:
                m+=(int(s.pos)-1)
                position='H'
            elif (int(s.pos)+m)>=267:
                if (int(s.pos)+m)>282:
                    position='P'
                else:
                    position='T'
            else:
                m=len(s.seq)
                position='M' 
            tooShort=True
        else:
            tooShort=True
        if tooShort==True:
            return('','','fully_clipped',position)
        else:
            if s.orientation=='0':
                return(seq,qual,'partially_clipped',position)
            elif s.orientation=='1':
                qual=qual[::-1]
                return(reverseComplement(seq),qual,'partially_clipped',position)
    else:
        if s.cigar!='*':
            if s.orientation=='0':
                return(s.seq,s.qual,'unclipped','N')
            elif s.orientation=='1':
                qual=s.qual[::-1]
                return(reverseComplement(s.seq),qual,'unclipped','N')
        else:
            return(s.seq,s.qual,'unclipped','N')
        
def clipOutAluPFU(cigar,pos,seq,qual):
    if cigar=='*':
        return(seq,qual,'unclipped','N')
    elif (int(pos)>270):
        return(seq,qual,'unclipped','N')
    else:
        cigar_span_list,cigar_character_list=splitCigar(cigar)
        mcount,n=0,0
        while n<len(cigar_character_list):
            if cigar_character_list[n]=='M':
                mcount+=int(cigar_span_list[n])
            n+=1
        if int(pos)<=50:
            position='H'
        elif (int(pos)+mcount)>=250:
            if (int(pos)+mcount)>=282:
                position='P'
            else:
                position='T'
        else:
            position='M'
    if 'S' not in cigar:
        return('','','fully_clipped',position)
    else:
        cigar_span_list,cigar_character_list=splitCigar(cigar)
        frontEnd,backEnd,f,b,part=0,0,0,-1,False
        while cigar_character_list[f]!='M':
            frontEnd+=int(cigar_span_list[f])
            f+=1
        while cigar_character_list[b]!='M':
            backEnd+=int(cigar_span_list[b])
            b-=1
        if ((frontEnd-int(pos))>=15) and (backEnd==0):
            part=True
            return(seq[:frontEnd],qual[:frontEnd],'partially_clipped',position)
        elif (backEnd>=15) and ((int(pos)+(len(seq))-(backEnd+frontEnd))>=267) and (int(pos)<=270):
            part=True
            return(seq[backEnd:],qual[backEnd:],'partially_clipped',position)
        if part==False:
            return('','','fully_clipped',position)

def processSamA(path2Psam,path2Ssam,reference):
    with open(reference+'.DS.info','w') as DS_info:
        with open(reference+'.SP.info','w') as SP_info:
            with open(reference+'.DS.fq','w') as fastqDS:
                with open(reference+'.SP.fq','w') as fastqSP:
                    with open(path2Psam,'r') as samfileP:
                        with open(path2Ssam,'r') as samfileS:
                            pe,ds,sp=1,1,1
                            temp,n=[],0
                            for line in samfileP:
                                if line[0]!='@':
                                    if n!=0:
                                        if line.split()[0]!=n:
                                            if len(temp)==2:
                                                for t in temp:
                                                    if 'OP:i:1' in t:
                                                        fwd_seq,fwd_qual,fwd_tag,positionF=clipOutAluPFU(t.split()[5],t.split()[3],t.split()[9],t.split()[10])
                                                        forward_read=t
                                                    elif 'OP:i:2' in t:
                                                        rvs_seq,rvs_qual,rvs_tag,positionR=clipOutAluPFU(t.split()[5],t.split()[3],t.split()[9],t.split()[10])
                                                        reverse_read=t
                                                posPlus=positionF+positionR
                                                if (posPlus!='NN'):
                                                    if (fwd_tag!='fully_clipped') and (rvs_tag!='fully_clipped') and (posPlus!='NN'):
                                                        nothin='nothin'
                                                    elif (fwd_tag=='fully_clipped'):
                                                        if rvs_tag=='partially_clipped':
                                                            rSam=sam(reverse_read)
                                                            fastqSP.write('@'+zeroThehero(12,sp)+'\n'+rvs_seq+'\n+\n'+rvs_qual+'\n')
                                                            SP_info.write(zeroThehero(12,sp)+'\t'+rSam.qname+'\t'+rSam.pos+'\t'+rSam.cigar+'\t'+rSam.seq+'\t'+positionR+'\t'+rSam.orientation+'\n')
                                                            sp+=1
                                                        elif rvs_tag=='unclipped':
                                                            fSam=sam(forward_read)
                                                            fastqDS.write('@'+zeroThehero(12,ds)+'\n'+rvs_seq+'\n+\n'+rvs_qual+'\n')
                                                            DS_info.write(zeroThehero(12,ds)+'\t'+fSam.qname+'\t'+fSam.pos+'\t'+fSam.cigar+'\t'+fSam.seq+'\t'+positionF+'\t'+fSam.orientation+'\n')
                                                            ds+=1
                                                    if rvs_tag=='fully_clipped':
                                                        if fwd_tag=='partially_clipped':
                                                            fSam=sam(forward_read)
                                                            fastqSP.write('@'+zeroThehero(12,sp)+'\n'+fwd_seq+'\n+\n'+fwd_qual+'\n')
                                                            SP_info.write(zeroThehero(12,sp)+'\t'+fSam.qname+'\t'+fSam.pos+'\t'+fSam.cigar+'\t'+fSam.seq+'\t'+positionF+'\t'+fSam.orientation+'\n')
                                                            sp+=1
                                                        elif fwd_tag=='unclipped':
                                                            rSam=sam(reverse_read)
                                                            fastqDS.write('@'+zeroThehero(12,ds)+'\n'+fwd_seq+'\n+\n'+fwd_qual+'\n')
                                                            DS_info.write(zeroThehero(12,ds)+'\t'+rSam.qname+'\t'+rSam.pos+'\t'+rSam.cigar+'\t'+rSam.seq+'\t'+positionR+'\t'+rSam.orientation+'\n')
                                                            ds+=1
                                            temp=[line.strip('\n')]
                                            n=line.split()[0]
                                        else:
                                            temp.append(line.strip('\n'))
                                            n=line.split()[0]
                                    else:
                                        temp.append(line.strip('\n'))
                                        n=line.split()[0]
                            for sline in samfileS:
                                if sline[0]!='@':
                                    fwd_seq,fwd_qual,fwd_tag,positionF=clipOutAluPFU(sline.split()[5],sline.split()[3],sline.split()[9],sline.split()[10])
                                    if fwd_tag=='partially_clipped':
                                        fSam=sam(sline)
                                        fastqSP.write('@'+zeroThehero(12,sp)+'\n'+fwd_seq+'\n+\n'+fwd_qual+'\n')
                                        SP_info.write(zeroThehero(12,sp)+'\t'+fSam.qname+'\t'+fSam.pos+'\t'+fSam.cigar+'\t'+fSam.seq+'\t'+positionF+'\t'+fSam.orientation+'\n')
                                        sp+=1

    

def processSam1(path2Psam,path2Ssam,reference):
    """The Purpose is to: 1.) Condense Each Line by Removing the SEQ and QUAL 2.) Separate the Forward and Reverse 3.) Parse out only the mapped reads""" 
    with open(reference+'.PE.fq','w') as fastqPE:
        with open(reference+'.PE.info','w') as PE_info:
            with open(path2Psam,'r') as samfileP:
                with open(path2Ssam,'r') as samfileS:
                    with open(reference+'_fwd.dup.sam','w') as fwd_dup:
                        with open(reference+'_rvs.dup.sam','w') as rvs_dup:
                            tempF,tempR,n,mapped,f,r=[],[],1,False,0,0
                            pe,ds,sp=1,1,1
                            for line in samfileP:
                                if line[0]!='@':
                                    s=sam(line)
                                    if int(s.qname)==n:
                                        if 'OP:i:1' in line:
                                            tempF.append([line])
                                            f+=1
                                            if (s.cigar!='*') and (int(s.pos)<=270):
                                                mapped=True
                                        elif 'OP:i:2' in line:
                                            tempR.append([line])
                                            r+=1
                                            if (s.cigar!='*') and (int(s.pos)<=270):
                                                mapped=True
                                    else:
                                        if mapped==True:
                                            if (f==1) and (r==1):
                                                fSam,rSam=sam(tempF[0][0]),sam(tempR[0][0])
                                                #forward.write(tempF[0][0])
                                                #reverse.write(tempR[0][0])
                                                fwd_seq,fwd_qual,fwd_tag,positionF=aluClipped(tempF[0][0])
                                                rvs_seq,rvs_qual,rvs_tag,positionR=aluClipped(tempR[0][0])
                                                posPlus=positionF+positionR
                                                if (fwd_tag!='fully_clipped') and (rvs_tag!='fully_clipped') and (posPlus!='NN'):
                                                    fastqPE.write('@'+zeroThehero(12,pe)+'\n'+fwd_seq+'\n+\n'+fwd_qual+'\n@'+zeroThehero(12,pe)+'\n'+rvs_seq+'\n+\n'+rvs_qual+'\n')
                                                    PE_info.write(zeroThehero(12,pe)+'\t'+fSam.qname+'\t'+fSam.pos+'\t'+fSam.cigar+'\t'+fSam.seq+'\t'+positionF+'\t'+fSam.orientation+'\n'+zeroThehero(12,pe)+'\t'+rSam.qname+'\t'+rSam.pos+'\t'+rSam.cigar+'\t'+rSam.seq+'\t'+positionR+'\t'+rSam.orientation+'\n')
                                                    pe+=1
                                        n,tempF,tempR,mapped,f,r=int(s.qname),[],[],False,0,0
                                        if 'OP:i:1' in line:
                                            tempF.append([line])
                                            f+=1
                                            if s.cigar!='*':
                                                mapped=True
                                        elif 'OP:i:2' in line:
                                            tempR.append([line])
                                            r+=1
                                            if s.cigar!='*':
                                                mapped=True
                            
def matchedReads(info1,info2,read):
    r,i1,i2,match1,match2=sam(read),info(info1),info(info2),False,False
    if (r.seq in i1.seq) or (r.seq in reverseComplement(i1.seq)):
        match1=True
        final=info1
    if (r.seq in i2.seq) or (r.seq in reverseComplement(i2.seq)):
        match2=True
        final=info2
    if (match1==True) and (match2==True):
        double=True
    else:
        double=False
    return(final,double)

def mapSpan(cigar):
    span,char=splitCigar(cigar)
    i,mapp=char.index('M'),0
    while i<len(span):
        mapp+=int(span[i])
        i+=1
    return(str(mapp))
    

def predictLocus(read,in_fo,read_type):
    r,i=sam(read),info(in_fo)
    dist_Alu,dist_ref,aluMap=distanceM(i.cigar),distanceM(r.cigar),0
    alu_span,alu_char=splitCigar(i.cigar)
    alu_index=alu_char.index('M')
    ref_span,ref_char=splitCigar(r.cigar)
    refMap=0
    for char in range(len(ref_char)):
        if ref_char[char]=='M':
            refMap+=int(ref_span[char])    
    while alu_index<len(alu_char):
        aluMap+=int(alu_span[alu_index])
        alu_index+=1
    if read_type!='D':
        if r.orientation==i.orientation:
            orientation='0'
            if int(i.pos)<4:
                start=str(int(r.pos) - int(dist_ref) + int(dist_Alu) - int(i.pos))
                return(r.rname+':'+start+'\t'+i.pos+'\t'+i.flag+'\t'+read_type+'\t'+i.neWname+':'+i.olDname+'\t'+i.seq+'\t'+orientation)
            elif ((int(i.pos)+aluMap)>267) and (alu_char[-1]!='M'):
                start=r.pos
                return(r.rname+':'+start+'\t'+i.pos+'\t'+i.flag+'\t'+read_type+'\t'+i.neWname+':'+i.olDname+'\t'+i.seq+'\t'+orientation)
            else:
                return('Null')
        else:
            if int(i.pos)<4:
                start = str(r.pos)
                orientation='1'
                return(r.rname+':'+start+'\t'+i.pos+'\t'+i.flag+'\t'+read_type+'\t'+i.neWname+':'+i.olDname+'\t'+i.seq+'\t'+orientation)
            elif ((int(i.pos)+aluMap)>267) and (alu_char[-1]!='M'):
                start=str(int(r.pos)+refMap)
                orientation='1'
                return(r.rname+':'+start+'\t'+i.pos+'\t'+i.flag+'\t'+read_type+'\t'+i.neWname+':'+i.olDname+'\t'+i.seq+'\t'+orientation)
            else:
                return('Null')
    else:
        if r.orientation=='0':
            start=str(int(r.pos)+aluMap)
            if i.orientation=='0':
                orientation='1'
            elif i.orientation=='1':
                orientation='0'
        elif r.orientation=='1':
            start=r.pos
            if i.orientation=='0':
                orientation='0'
            elif i.orientation=='1':
                orientation='1'
        return(r.rname+':'+start+'\t'+i.pos+'\t'+i.flag+'\t'+read_type+'\t'+i.neWname+':'+i.olDname+'\t'+i.seq+'\t'+orientation)
                                                                    
def pairedEndPredict(read1,read2,info1,info2):
    r1,r2,i1,i2=sam(read1),sam(read2),info(info1),info(info2)
    predictions=[]
    info_r1,double1=matchedReads(info1,info2,read1)
    info_r2,double2=matchedReads(info1,info2,read2)
    inf1,inf2=info(info_r1),info(info_r2)
    if (int(r1.mapq)>=1) and (int(r2.mapq)>=1):
        if ('M' in inf1.cigar):
            if double1==False:
                predictions.append(predictLocus(read1,info_r1,'P'))
        if ('M' in inf2.cigar):
            if double2==False:
                predictions.append(predictLocus(read2,info_r2,'P'))
    else:
        if (int(r1.mapq)>=1):
            if ('M' in inf1.cigar):
                if double1==False:
                    predictions.append(predictLocus(read1,info_r1,'S'))
            if ('M' in inf2.cigar):
                if double2==False:
                    predictions.append(predictLocus(read1,info_r2,'D'))
        elif (int(r2.mapq)>=1):
            if ('M' in inf2.cigar):
                if double2==False:
                    predictions.append(predictLocus(read2,info_r2,'S'))
            if ('M' in inf1.cigar):
                if double1==False:
                    predictions.append(predictLocus(read2,info_r1,'D'))
    if predictions!=[]:
        return(predictions)
    else:
        return(['Null'])

def predictInsertionOG(head_tag,outfile):
    with open(head_tag+'.PE.sam','r') as pairedEnd:
        with open(head_tag+'.SP.sam','r') as splitRead:
            with open(head_tag+'.DS.sam','r') as discordant:
                with open(outfile,'w') as final:
                    temp,finaList,ds_info,sp_info,pe_info=[],[],lister(head_tag+'.DS.info'),lister(head_tag+'.SP.info'),lister(head_tag+'.PE.info')
                    for line in pairedEnd:
                        if line[0]!='@':
                            if temp==[]:
                                temp.append(line)
                            else:
                                r1,present=sam(temp[0]),sam(line)
                                if r1.qname==present.qname:
                                    temp.append(line)
                                elif len(temp)==2:
                                    r2=sam(temp[1])
                                    predictions=pairedEndPredict(temp[0],temp[1],pe_info[(int(r2.qname)*2)-2],pe_info[(int(r2.qname)*2)-1])
                                    for p in predictions:
                                        if p!='Null':
                                            finaList.append(p)
                                    temp=[line]
                                else:
                                    temp=[line]
                    for line in splitRead:
                        if line[0]!='@':
                            if line.split()[2]!='*':
                                sprd=sam(line)
                                if int(sprd.mapq)>=1:
                                    prediction=predictLocus(line,sp_info[(int(sprd.qname)-1)],'S')
                                    if prediction!='Null':
                                        finaList.append(prediction)
                    for line in discordant:
                        if line[0]!='@':
                            if line.split()[2]!='*':
                                dscr=sam(line)
                                if int(sprd.mapq)>=1:
                                    prediction=predictLocus(line,ds_info[(int(dscr.qname)-1)],'D')
                                    if prediction!='Null':
                                        finaList.append(prediction)
                    for f in finaList:
                        try:
                            final.write(f+'\n')
                        except TypeError:
                            print(f,file=sys.stderr)
                                                                    
def predictInsertion(head,tag,outfile):
    with open(head+'.'+tag+'.sam','r') as infile:
        with open(outfile,'a') as final:
            temp,finaList,info=[],[],lister(head+'.'+tag+'.info')
            if tag=='PE':                
                for line in infile:
                    if line[0]!='@':
                        if temp==[]:
                            temp.append(line)
                        else:
                            r1,present=sam(temp[0]),sam(line)
                            if r1.qname==present.qname:
                                temp.append(line)
                            elif len(temp)==2:
                                r2=sam(temp[1])
                                predictions=pairedEndPredict(temp[0],temp[1],info[(int(r2.qname)*2)-2],info[(int(r2.qname)*2)-1])
                                for p in predictions:
                                    if p!='Null':
                                        final.write(p+'\n')
                                temp=[line]
                            else:
                                temp=[line]
            elif tag=='SP':
                for line in infile:
                    if line[0]!='@':
                        if line.split()[2]!='*':
                            sprd=sam(line)
                            if int(sprd.mapq)>=1:
                                prediction=predictLocus(line,info[(int(sprd.qname)-1)],'S')
                                if prediction!='Null':
                                    final.write(prediction+'\n')
            elif tag=='DS':
                for line in infile:
                    if line[0]!='@':
                        if line.split()[2]!='*':
                            dscr=sam(line)
                            if int(dscr.mapq)>=1:
                                prediction=predictLocus(line,info[(int(dscr.qname)-1)],'D')
                                if prediction!='Null':
                                    final.write(prediction+'\n')
                                
                                        
                                    
                    

def removeDupsSAM(path2sam,outpath):
    with open(path2sam[:-3]+'rpt','w') as sam:
        with open(path2sam[:-3]+'uniq','w') as uniq:
            subprocess.call(['uniq', '-w', '12', '-D', path2sam],stdout=sam)
            subprocess.call(['uniq', '-w', '12', '-u', path2sam],stdout=uniq)
    previous='nothing'
    with open(path2sam[:-3]+'rpt','r') as enput:
        with open(path2sam[:-3]+'norpt','w') as output:
            for line in enput:
                if line[:12]!=previous:
                    output.write(line)
                previous=line[:12]
    with open(path2sam[:-3]+'cat','w') as cat:
        with open(outpath,'w') as out:
            subprocess.call(['cat', path2sam[:-3]+'uniq', path2sam[:-3]+'norpt'],stdout=cat)
            subprocess.call(['sort','-nk1',path2sam[:-3]+'cat'],stdout=out)
    subprocess.call(['rm -f '+path2sam[:-3]+'rpt'], shell=True)
    subprocess.call(['rm -f '+path2sam[:-3]+'norpt'], shell=True)
    subprocess.call(['rm -f '+path2sam[:-3]+'cat'], shell=True)
    subprocess.call(['rm -f '+path2sam[:-3]+'uniq'], shell=True)
    
def aluLine_process(line,s_or_l):
    aluSplit=line.split()
    tag,cigar=aluSplit[1],aluSplit[5]
    orient=orientConverter(tag)
    aluPos,mapping_score=aluSplit[3],aluSplit[4]
    tags=aluSplit[-1].split('_')
    for t in tags:
        if t[:4]=='MD:Z':
            tagMD=t
    if orient=='0':
        readPos=distanceM(cigar)
        cigar_span_list,cigar_character_list=splitCigar(cigar)
    elif orient=='1':
        cigar_span_list,cigar_character_list=splitCigar(cigar)
        readPos=distanceM(flipCigar(cigar_span_list,cigar_character_list))
        cigar_span_list,cigar_character_list=splitCigar(flipCigar(cigar_span_list,cigar_character_list))
    last_M_index,first_M_index=''.join(cigar_character_list).rfind('M'),cigar_character_list.index('M')
    if last_M_index==first_M_index:
        indexRange_MpBs=(str(readPos)+':'+str(int(readPos)+int(cigar_span_list[first_M_index])))
    else:
        total,i=0,first_M_index
        while i<=last_M_index:
            total+=int(cigar_span_list[i])
            i+=1
        indexRange_MpBs=str(readPos)+':'+str(int(readPos)+total)
    return(str(aluPos)+'_'+str(int(readPos)+1)+'_'+s_or_l+'_'+orient+'_'+mapping_score+'_'+cigar+'_'+indexRange_MpBs+'_'+tagMD)


def clippers(temp,seq):
    seq=seq.strip('\n')
    temp_split=temp[0].split()[1].split('DEW')
    id_N,coZ=temp[0].split()[0],temp_split[0]
    if temp_split[2]!='N/A':
        info_split=temp_split[2].split('_')
    else:
        info_split=temp_split[1].split('_')
    alu_pos,read_pos,orient,index_range,cigar=int(info_split[0]),int(info_split[1]),info_split[3],info_split[6].split(':'),info_split[5]
    mappedBases=int(index_range[1])-int(index_range[0])
    tagMD=info_split[-1]
    tag='N/A'
    if alu_pos<=70:        
        if orient=='0':
            if (read_pos-alu_pos)>=30:
                seq=seq[:(read_pos-alu_pos+1)]
                tag='Head'
            else:
                seq='NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
        elif orient=='1':
            if mappedBases<=70:
                seq=seq[mappedBases-1:]
                tag='Head'
            else:
                seq='NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
    elif alu_pos>=180:        
        if orient=='0':
            if int(index_range[1])<=70:
                seq=seq[int(index_range[1]):]
                tag='Tail0'                 
                if (alu_pos+mappedBases)>=285:
                    tag='Tail1'
            else:
                seq='NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
        elif orient=='1':
            if int(index_range[0])>=30:
                seq=seq[:int(index_range[0])]
                tag='Tail0'
                if (alu_pos+mappedBases)>=285:
                    tag='Tail1'
            else:
                seq='NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
    else:
        seq='NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
    additionalInfo='\tCO:Z:aluMD='+tagMD+'\tCO:Z:orientation='+orient+'\tCO:Z:aluHT='+tag+'\tCO:Z:cigar='+cigar
    return(id_N+'\t'+coZ+'\tCO:Z:aluPos='+str(alu_pos)+'\tCO:Z:readPos='+str(read_pos)+additionalInfo+'\n'+seq+'\n')
    
def reformatFQ_wAluData(path2fq,path2aluS,fwd_or_rv):
    with open(path2fq,'r') as alu_fq:
        alu_S=lister(path2aluS)
        with open('/home/vallmer/papio_pipe/30388/aligned/clippedAlu.'+fwd_or_rv+'.fq','w') as outFastQ:
            temp,clipped,seq,qual=[],False,False,False
            for line in alu_fq:
                if (line[0]=='@') and (len(line.split()[0])==13):
                    readN=int(line.split()[0][1:])-1
                    if alu_S[readN].split()[5] != '*':
                        aluS_tag=aluLine_process(alu_S[readN],'S')
                        clipped=True
                    else:
                        aluS_tag='N/A'
                    if alu_L[readN].split()[5] != '*':
                        aluL_tag=aluLine_process(alu_L[readN],'L')
                        clipped=True
                    else:
                        aluL_tag='N/A'    
                    if clipped==True:
                        temp.append(line.split()[0]+'\t'+line.split()[1]+'DEW'+aluS_tag+'DEW'+aluL_tag)
                        seq,clipped=True,False
                    else:
                        outFastQ.write(line)
                elif (seq==False) and (qual==False):
                    outFastQ.write(line)
                elif seq==True:
                    clipString=clippers(temp,line)
                    seqLen=len(clipString.split('\n')[1])
                    outFastQ.write(clippers(temp,line))
                    seq,temp,=False,[]
                    qual='Plus'
                elif qual=='Plus':
                    qual=True
                    outFastQ.write(line)
                elif qual==True:
                    outFastQ.write(line[:seqLen]+'\n')
                    qual=False
                    
def toSplit_or_not2split(path2fwd,path2rvs):
    clipped_fwd,masked_fwd,unaltered_fwd=set(),set(),set()
    clipped_rvs,masked_rvs,unaltered_rvs=set(),set(),set()
    with open(path2fwd,'r') as fqF:
        for line in fqF:
            if 'CO:Z:aluHT=' in line:
                if 'N/A' not in line:
                    clipped_fwd.add(int(line.split()[0][1:]))
                else:
                    masked_fwd.add(int(line.split()[0][1:]))
            elif 'CO:Z:F' in line:
                unaltered_fwd.add(int(line.split()[0][1:]))
    with open(path2rvs,'r') as fqR:
        for line in fqR:
            if 'CO:Z:aluHT=' in line:
                if 'N/A' not in line:
                    clipped_rvs.add(int(line.split()[0][1:]))
                else:
                    masked_rvs.add(int(line.split()[0][1:]))
            elif 'CO:Z:R' in line:
                unaltered_rvs.add(int(line.split()[0][1:]))
    final=set()
    for readNum in clipped_fwd:
        if (readNum in unaltered_rvs):
            final.add(readNum)
    for readNum in clipped_rvs:
        if (readNum in unaltered_fwd):
            final.add(readNum)
    return(final)

def parseSplitFQ(path2fq,rnSet,path2outFQ):
    head,seq,plus=False,False,False
    with open(path2fq,'r') as inFq:
        with open(path2outFQ,'w') as outFq:
            for line in inFq:
                if head==True:
                    outFq.write(line)
                    if seq==True:
                        if plus==True:
                            head,seq,plus=False,False,False
                        else:
                            plus=True
                    else:
                        seq=True
                elif 'CO:Z:' in line:
                    if int(line.split()[0][1:]) in rnSet:
                        head=True
                        outFq.write(line)
                        
def sortPredict(path2predict,outfile,chromLen):
    with open(outfile,'w') as out:
        sorting_piles,chrom_list=[],[]
        for i in range(1,(int(chromLen)+1):
            sorting_piles.append([])
            chrom_list.append(str(i))
        for p in lister(path2predict):
            candidate=predict(p)
            if candidate.chromosome in chrom_list:
                cols3plus=p.split()[1:]
                sorting_piles[int(candidate.chromosome)-1].append(candidate.chromosome+'\t'+candidate.locus+'\t'+'\t'.join(cols3plus))
        for s in sorting_piles:
            s_sort=(sortMe(s,2))
            for item in s_sort:
                out.write(item+'\n')

def cleanup(label,orgID):
    if label!='A':
        subprocess.call(['rm -f '+orgID+'.'+label+'.DS.fq'], shell=True)
        subprocess.call(['rm -f '+orgID+'.'+label+'.DS.info'], shell=True)
        subprocess.call(['rm -f '+orgID+'.'+label+'.DS.sam'], shell=True)
        subprocess.call(['rm -f '+orgID+'.'+label+'.PE.sam'], shell=True)
        subprocess.call(['rm -f '+orgID+'.'+label+'.PE.fq'], shell=True)
        subprocess.call(['rm -f '+orgID+'.'+label+'.PE.info'], shell=True)
        subprocess.call(['rm -f '+orgID+'.'+label+'.SP.sam'], shell=True)
        subprocess.call(['rm -f '+orgID+'.'+label+'.SP.fq'], shell=True)
        subprocess.call(['rm -f '+orgID+'.'+label+'.SP.info'], shell=True)
        subprocess.call(['rm -f '+orgID+'.'+label+'_fwd.dup.sam'], shell=True)
        subprocess.call(['rm -f '+orgID+'.'+label+'_rvs.dup.sam'], shell=True)
        subprocess.call(['rm -f '+orgID+'alu.P.'+label+'.mapped.sam'], shell=True)        
        subprocess.call(['rm -f '+orgID+'alu.P.'+label+'.sam'], shell=True)
        subprocess.call(['rm -f '+orgID+'alu.S.'+label+'.mapped.sam'], shell=True)        
        subprocess.call(['rm -f '+orgID+'alu.S.'+label+'.sam'], shell=True)
        subprocess.call(['rm -f '+orgID+'.'+label+'.predict'], shell=True)
    else:
        subprocess.call(['rm -f '+orgID+'.P.polyA.sam'], shell=True)
        subprocess.call(['rm -f '+orgID+'.S.polyA.sam'], shell=True)
        subprocess.call(['rm -f '+orgID+'.polyA.fq'], shell=True)
        subprocess.call(['rm -f '+orgID+'.polyA.info'], shell=True)        
        subprocess.call(['rm -f '+orgID+'.polyA.sam'], shell=True)

editFQ_title(path2_PE_fq,path2outdir+orgID+'.PE.fq')  
editFQ_title(path2_SE_fq,path2outdir+orgID+'.SE.fq')
subprocess.call(['rm -f clipped_paired.fq'], shell=True)
subprocess.call(['rm -f clipped_single.fq'], shell=True)

bwaMeM(path2aluConsensus,path2outdir+orgID+'.PE.fq','P','S',path2outdir+orgID+'alu.P.S.sam')
bwaMeM(path2aluConsensus,path2outdir+orgID+'.SE.fq','S','S',path2outdir+orgID+'alu.S.S.sam')
mappOnly(path2outdir+orgID+'alu.P.S.sam',path2outdir+orgID+'alu.P.S.mapped.sam')
mappOnly(path2outdir+orgID+'alu.S.S.sam',path2outdir+orgID+'alu.S.S.mapped.sam')
processSam1((path2outdir+orgID+'alu.P.S.sam'),(path2outdir+orgID+'alu.S.S.sam'),(orgID+'.S'))
processSamA((path2outdir+orgID+'alu.P.S.mapped.sam'),(path2outdir+orgID+'alu.S.S.mapped.sam'),(orgID+'.S'))
bowtie2(path2referenceGenome,(orgID+'.S.PE.fq'),'P',(orgID+'.S.PE.sam'))
bowtie2(path2referenceGenome,(orgID+'.S.DS.fq'),'S',(orgID+'.S.DS.sam'))
bowtie2(path2referenceGenome,(orgID+'.S.SP.fq'),'S',(orgID+'.S.SP.sam'))
predictInsertionOG(orgID+'.S',orgID+'.S.predict')
sortPredict(orgID+'.S.predict',orgID+'.S.predict.sort',chromLen)
cleanup('S',orgID)

bwaMeM(path2aluConsensus,path2outdir+orgID+'.PE.fq','P','L',path2outdir+orgID+'alu.P.L.sam')
bwaMeM(path2aluConsensus,path2outdir+orgID+'.SE.fq','S','L',path2outdir+orgID+'alu.S.L.sam')
mappOnly(path2outdir+orgID+'alu.P.L.sam',path2outdir+orgID+'alu.P.L.mapped.sam')
mappOnly(path2outdir+orgID+'alu.S.L.sam',path2outdir+orgID+'alu.S.L.mapped.sam')
processSam1((path2outdir+orgID+'alu.P.L.sam'),(path2outdir+orgID+'alu.S.L.sam'),(orgID+'.L'))
processSamA((path2outdir+orgID+'alu.P.L.mapped.sam'),(path2outdir+orgID+'alu.S.L.mapped.sam'),(orgID+'.L'))
bowtie2(path2referenceGenome,(orgID+'.L.PE.fq'),'P',(orgID+'.L.PE.sam'))
bowtie2(path2referenceGenome,(orgID+'.L.DS.fq'),'S',(orgID+'.L.DS.sam'))
bowtie2(path2referenceGenome,(orgID+'.L.SP.fq'),'S',(orgID+'.L.SP.sam'))
predictInsertionOG(orgID+'.L',orgID+'.L.predict')
sortPredict(orgID+'.L.predict',orgID+'.L.predict.sort',chromLen)
cleanup('L',orgID)

bwaMeM(path2polyA,path2outdir+orgID+'.PE.fq','P','L',(orgID+'.P.polyA.sam'))
bwaMeM(path2polyA,path2outdir+orgID+'.SE.fq','S','L',(orgID+'.S.polyA.sam'))
processSam_polyA((orgID+'.P.polyA.sam'),(orgID+'.S.polyA.sam'),(orgID+'.polyA.fq'),(orgID+'.polyA.info'))
bowtie2(path2referenceGenome,(orgID+'.polyA.fq'),'S',(orgID+'.polyA.sam'))
predictPolyA((orgID+'.polyA.sam'),(orgID+'.polyA.predict'),(orgID+'.polyA.info'))
cleanup('A',orgID)
