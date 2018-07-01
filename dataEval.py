import sys
from essentials import *

def clusterS_predict(indirectory,path2orgList,outPath):
    with open(outPath,'w') as outPredict:
        for org in lister(path2orgList):
            print('Processing '+org,file=sys.stderr)
            wc_len= (subprocess.check_output('wc -l '+indirectory+org+'/'+org+'.S.predict.sort',shell=True)).decode('ascii').split()[0]
            n=0
            with open(indirectory+org+'/'+org+'.S.predict.sort','r') as inPredict:
                with open(indirectory+org+'/'+org+'.S.predict.sort.SnP','w') as outPredictSnP:
                    previous,total='',1
                    for line in inPredict:
                        if line.split()[4]!='D':
                            if previous=='':
                                previous=line
                            else:
                                prevPredict=prediction(previous)
                                currPredict=prediction(line)
                                if (prevPredict.chromosome==currPredict.chromosome) and (prevPredict.locus==currPredict.locus) and (prevPredict.orientation==currPredict.orientation):
                                    total+=1
                                    previous=line
                                else:
                                    outPredict.write(prevPredict.chromosome+'\t'+prevPredict.locus+'\t'+org+':'+str(total)+'\t'+prevPredict.aluLocus+'\t'+prevPredict.HTP+'\t'+prevPredict.SPD+'\t'+prevPredict.readname+'\t'+prevPredict.seq+'\t'+prevPredict.orientation+'\n')
                                    outPredictSnP.write(prevPredict.chromosome+'\t'+prevPredict.locus+'\t'+org+':'+str(total)+'\t'+prevPredict.aluLocus+'\t'+prevPredict.HTP+'\t'+prevPredict.SPD+'\t'+prevPredict.readname+'\t'+prevPredict.seq+'\t'+prevPredict.orientation+'\n')
                                    previous=line
                                    total=1
                        n+=1
                        if len(str(n))>=7:
                            if (str(n)[-6:]=='000000'):
                                print(str(n)+' of '+wc_len+' processed',file=sys.stderr)
            print(org+ ': Processing Complete\n',file=sys.stderr)
            
def clusterL_predict(indirectory,path2orgList,outPath):
    with open(outPath,'w') as outPredict:
        for org in lister(path2orgList):
            print('Processing '+org,file=sys.stderr)
            wc_len= (subprocess.check_output('wc -l '+indirectory+org+'/'+org+'.L.predict.sort',shell=True)).decode('ascii').split()[0]
            n=0
            with open(indirectory+org+'/'+org+'.L.predict.sort','r') as inPredict:
                with open(indirectory+org+'/'+org+'.L.predict.sort.SnP','w') as outPredictSnP:
                    previous,total='',1
                    for line in inPredict:
                        if line.split()[4]!='D':
                            if previous=='':
                                previous=line
                            else:
                                prevPredict=prediction(previous)
                                currPredict=prediction(line)
                                if (prevPredict.chromosome==currPredict.chromosome) and (prevPredict.locus==currPredict.locus) and (prevPredict.orientation==currPredict.orientation):
                                    total+=1
                                    previous=line
                                else:
                                    outPredict.write(prevPredict.chromosome+'\t'+prevPredict.locus+'\t'+org+':'+str(total)+'\t'+prevPredict.aluLocus+'\t'+prevPredict.HTP+'\t'+prevPredict.SPD+'\t'+prevPredict.readname+'\t'+prevPredict.seq+'\t'+prevPredict.orientation+'\n')
                                    outPredictSnP.write(prevPredict.chromosome+'\t'+prevPredict.locus+'\t'+org+':'+str(total)+'\t'+prevPredict.aluLocus+'\t'+prevPredict.HTP+'\t'+prevPredict.SPD+'\t'+prevPredict.readname+'\t'+prevPredict.seq+'\t'+prevPredict.orientation+'\n')
                                    previous=line
                                    total=1
                        n+=1
                        if len(str(n))>=7:
                            if (str(n)[-6:]=='000000'):
                                print(str(n)+' of '+wc_len+' processed',file=sys.stderr)
            print(org+ ': Processing Complete\n',file=sys.stderr)
            
def collectPredictions(path2predictDir,orgID,chromN):
    predictList=[]
    for i in range(int(chromN)):
        predictList.append([])
    orgIDlist=lister(orgID)
    for org in orgIDlist:
        with open(path2predictDir+org+'/'+org+'.S.predict.sort.SnP','r') as predict:
            with open(path2predictDir+org+'/'+org+'.L.predict.sort.SnP','r') as predictL:
                for line in predict:
                    split=line.split()
                    if (int(split[3])<=5) and (split[2].split(':')[1] != '1'):
                        predictList[int(split[0])-1].append(split[1]+split[-1]+'\t'+split[2])
                for line in predictL:
                    split=line.split()
                    predictList[int(split[0])-1].append(split[1]+split[-1]+'\t'+split[2].split(':')[0]+'L:'+split[2].split(':')[1])
    chrom=1
    for p in predictList:
        predictListFinal=[]
        p_sort=sortMe(p,1)
        prev,p,spec,tag='0',1,[],[]
        for pred in p_sort:
            ps=pred.split()[0]
            if (ps!=prev):
                if (prev!='0') and (p>1):
                    S=False
                    noL=[]
                    for species in spec:
                        if species[-1]!='L':
                            S=True
                            noL.append(species)
                        else:
                            noL.append(species[:-1])
                    if S==True:
                        predictListFinal.append(str(chrom)+'\t'+prev[:-1]+'\t'+'\t'.join(noL)+'\t'+prev[-1])
                prev=ps
                if pred.split()[1].split(':')[0][-1]=='L':
                    spec=[pred.split()[1].split(':')[0]]
                    tag=[pred.split()[1]]
                else:
                    spec=[pred.split()[1].split(':')[0]]
                    tag=[pred.split()[1]]
                p=1
            else:
                if pred.split()[1].split(':')[0][-1]!='L':
                    p+=1
                    spec.append(pred.split()[1].split(':')[0])
                    tag.append(pred.split()[1])
                elif (pred.split()[1].split(':')[0][:-1] not in spec):
                    p+=1
                    spec.append(pred.split()[1].split(':')[0])
                    tag.append(pred.split()[1])
        chrom+=1
        kinkos(predictListFinal)
            

clusterS_predict(sys.argv[1],sys.argv[2],sys.argv[3])
clusterL_predict(sys.argv[1],sys.argv[2],sys.argv[3])
collectPredictions(sys.argv[1],sys.argv[2],sys.argv[4])
