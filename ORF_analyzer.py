## ORF_analyzer: identify and summarize nucleotide and amino acid mutations present in a specified ORF across all sequences in a provided FASTA file
##
## Kevin Kuchinski
## Prystajecky Lab, BCCDC Public Health Laboratory
## Pathology and Laboratory Medicine, University of British Columbia
## kevin.kuchinski@bccdc.ca

import argparse
import subprocess
import itertools as it
import numpy as np
import pandas as pd

parser=argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, required=True, help="input FASTA file")
parser.add_argument("-r", "--reference", type=str, required=True, help="referene sequence in FASTA format")
parser.add_argument("-c", "--coord", type=str, required=True, help="nucleotide positions (1-indexed) of ORF (start-end)")
parser.add_argument("-o", "--output", type=str, required=True, help="short description to append to output files")
args=parser.parse_args()

print()
print('ORF_analyzer')

refSeqs={}
with open(args.reference) as referenceFile:
    referenceLines=referenceFile.readlines()
for i in range(len(referenceLines)):
    if referenceLines[i][0]=='>':
        header=referenceLines[i].rstrip()
        refSeqs[header]=''
    else:
        refSeqs[header]+=referenceLines[i].rstrip()

refHeader,refSeq=list(refSeqs.items())[0]
print()
print('Reference:',refHeader.lstrip('>'))

ORF_start=int(args.coord.split('-')[0])
ORF_end=int(args.coord.split('-')[1])
print('ORF coordinates:',ORF_start,'-',ORF_end)

start=max([1,ORF_start-1-200])
end=min([ORF_end-1+200,len(refSeq)])

ntSeqs={}
with open(args.input) as inputFile:
    inputLines=inputFile.readlines()
for i in range(len(inputLines)):
    if inputLines[i][0]=='>':
        header=inputLines[i].rstrip()
        ntSeqs[header]=''
    else:
        ntSeqs[header]+=inputLines[i].rstrip()

print()
print('Query sequences:',len(ntSeqs))

outputFile=open(args.output+'_regions.fa','w')
outputFile.write(refHeader+'\n')
outputFile.write(refSeq[start:end]+'\n')
for header,sequence in ntSeqs.items():
    totalCanonicalBases=sequence[start:end].count('A')+sequence[start:end].count('T')+sequence[start:end].count('G')+sequence[start:end].count('C')
    if len(sequence[start:end])==totalCanonicalBases:
        outputFile.write(header+'\n')
        outputFile.write(sequence[start:end]+'\n')
outputFile.close()

print()
print('Aligning nucleotide sequences...')
command='muscle -in '+args.output+'_regions.fa'+' -out '+args.output+'_regions_aligned.fa'+' -diags -maxiters 2'
subprocess.call(command,shell=True)

ntAlignments={}
with open(args.output+'_regions_aligned.fa') as inputFile:
    inputLines=inputFile.readlines()
for i in range(len(inputLines)):
    if inputLines[i][0]=='>':
        header=inputLines[i].rstrip()
        ntAlignments[header]=''
    else:
        ntAlignments[header]+=inputLines[i].rstrip()

print()
print('Nucleotide sequences aligned:',len(ntAlignments))

refNtPos={}
ri=start+1
ai=0
for s in ntAlignments[refHeader]:
    if s=='-':
        ai+=1
    else:
        refNtPos[ri]=ai
        ri+=1
        ai+=1

ntORFs={}
for header,alignment in ntAlignments.items():
    ntORFs[header]=alignment[refNtPos[ORF_start]:refNtPos[ORF_end]+1].replace('-','')

codons = {
    # 'M' - START, '*' - STOP                                                                                                                                                                             
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATA": "I", "ATT": "I", "ATC": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "*", "TAG": "*", "TGA": "*"
}

aaORFs={}
for header,ORF in ntORFs.items():
    aaORFs[header]=''
    for codon in [ORF[i:i+3] for i in list(range(len(ORF)))[::3]]:
        if codon in codons.keys():
            aaORFs[header]+=codons[codon]
        else:
            aaORFs[header]+='X'

outputFile=open(args.output+'_ORFs.fa','w')
for header,sequence in aaORFs.items():
    if 'X' not in sequence:
        outputFile.write(header+'\n')
        outputFile.write(sequence.split('*')[0]+'\n')
outputFile.close()

print()
print('Aligning amino acid sequences...')
command='muscle -in '+args.output+'_ORFs.fa'+' -out '+args.output+'_ORFs_aligned.fa'+' -maxiters 2 -seqtype protein'
subprocess.call(command,shell=True)                                                                                                                                                                     

aaAlignments={}
with open(args.output+'_ORFs_aligned.fa') as inputFile:
    inputLines=inputFile.readlines()
for i in range(len(inputLines)):
    if inputLines[i][0]=='>':
        header=inputLines[i].rstrip()
        aaAlignments[header]=''
    else:
        aaAlignments[header]+=inputLines[i].rstrip()

print()
print('Amino acid sequences aligned:',len(aaAlignments))

refAaPos={}
ri=1
ai=0
for s in aaAlignments[refHeader]:
    if s=='-':
        ai+=1
    else:
        refAaPos[ri]=ai
        ri+=1
        ai+=1

ntMutations={}
for header,alignment in ntAlignments.items():
    ntMutations[header]=[]
    for pos in range(ORF_start,ORF_end+1):
        refNt=ntAlignments[refHeader][refNtPos[pos]]
        subNt=alignment[refNtPos[pos]]
        if subNt!=refNt and 'N' not in [refNt,subNt]:
            ntMutations[header].append(refNt+str(pos)+subNt)

aaMutations={}
for header,alignment in aaAlignments.items():
    aaMutations[header]=[]
    for pos in range(1,len(refAaPos)+1):
        refAa=aaAlignments[refHeader][refAaPos[pos]]
        subAa=alignment[refAaPos[pos]]
        if subAa!=refAa and 'X' not in [refAa,subAa]:
            aaMutations[header].append(refAa+str(pos)+subAa)

## REPORT NT MUTATIONS
alphabet='ABCDEFGHIJKLMNOPQRSTUVWXYZ-'
outData=pd.Series([len(ntMutations[header]) for header in ntMutations.keys()]).value_counts().reset_index()
outData.columns=['SNP_num','count']
outData['freq']=outData.apply(lambda row: round(row['count']/len(ntAlignments),3), axis=1)
outData.to_csv(args.output+'_SNP_num.tsv',sep='\t',index=False)
combinationsThreshold=int(pd.Series([len(ntMutations[header]) for header in ntMutations.keys()]).quantile(0.9))

outData=pd.Series([item.lstrip(alphabet).rstrip(alphabet) for sublist in ntMutations.values() for item in sublist]).value_counts().reset_index()
outData.columns=['SNP_pos','count']
outData['freq']=outData.apply(lambda row: round(row['count']/len(ntAlignments),3), axis=1)
outData.to_csv(args.output+'_SNP_pos.tsv',sep='\t',index=False)

outData=pd.DataFrame()
for k in range(1,combinationsThreshold+1):
    data=pd.Series([item for sublist in ntMutations.values() for item in it.combinations(sublist,k)]).value_counts().reset_index()
    data.columns=['SNP_combo','count']
    data['freq']=data.apply(lambda row: round(row['count']/len(ntAlignments),3), axis=1)
    data['num_SNPs']=[k]*len(data)
    data.SNP_combo=[', '.join(combo) for combo in data.SNP_combo]
    outData=pd.concat([outData,data],sort=True)
outData=outData[['SNP_combo','num_SNPs','count','freq']]
outData.sort_values(by='freq', ascending=False).to_csv(args.output+'_SNP_combo.tsv',sep='\t',index=False)

### REPORT AA MUTATIONS
outData=pd.Series([len(aaMutations[header]) for header in aaMutations.keys()]).value_counts().reset_index()
outData.columns=['sub_num','count']
outData['freq']=outData.apply(lambda row: round(row['count']/len(aaAlignments),3), axis=1)
outData.to_csv(args.output+'_sub_num.tsv',sep='\t',index=False)
combinationsThreshold=int(pd.Series([len(aaMutations[header]) for header in aaMutations.keys()]).quantile(0.9))

outData=pd.Series([item.lstrip(alphabet).rstrip(alphabet) for sublist in aaMutations.values() for item in sublist]).value_counts().reset_index()
outData.columns=['sub_pos','count']
outData['freq']=outData.apply(lambda row: round(row['count']/len(aaAlignments),3), axis=1)
outData.to_csv(args.output+'_sub_pos.tsv',sep='\t',index=False)

outData=pd.DataFrame()
for k in range(1,combinationsThreshold+1):
    data=pd.Series([item for sublist in aaMutations.values() for item in it.combinations(sublist,k)]).value_counts().reset_index()
    data.columns=['sub_combo','count']
    data['freq']=data.apply(lambda row: round(row['count']/len(aaAlignments),3), axis=1)
    data['num_subs']=[k]*len(data)
    data.sub_combo=[', '.join(combo) for combo in data.sub_combo]
    outData=pd.concat([outData,data],sort=True)
outData=outData[['sub_combo','num_subs','count','freq']]
outData.sort_values(by='freq', ascending=False).to_csv(args.output+'_sub_combo.tsv',sep='\t',index=False)
