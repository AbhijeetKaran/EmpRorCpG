import os
import sys
import py2bit as tbit
import numpy as np
import pandas as pd 


class SequenceFetcher:
    def __init__(self,filename,sformat):
        self.filename = filename
        self.sformat = sformat
        self.filename_index = f"{filename.split('.')[0]}.fai"
        
        # Validate the input provided 
        self.validate()
        self.index = None

        # read index of the file
        if sformat == 'f':
            self.index = self.readFastaIndex()
        elif sformat == 'b':
            self.index = self.read2bitIndex()

    # Validates provided parameters
    def validate(self):
        if not os.path.exists(self.filename):
            print(f"Sequence file do not exist")
            exit()

        if self.sformat == None:
            print(f"File format is not provided, options are f or b [f for fasta and b for 2bit]")
            exit()

    # Read index from fasta file
    def readFastaIndex(self):
        if not os.path.exists(self.filename_index):
            print("Index file do not exist")
            self.createFastaIndex()

        sequenceIndex = dict()
        fh = open(self.filename_index,"r")
        for line in fh:
            line = line.strip()
            line_con = line.split('\t')
            params = []
            for i in range(1,len(line_con)):
                params.append(int(line_con[i]))
            
            sequenceIndex[line_con[0]] = params
        fh.close()
        return sequenceIndex

    # Create index from fasta file
    def createFastaIndex(self):
        print(f"Creating index file for {os.path.basename(self.filename)}")
        fh = open(self.filename,'r')
        fhw = open(self.filename_index,'w')
        data = fh.read().strip()
        seqs = data.split(">")
        seqs.pop(0)
        totalSeqs = len(seqs)
        nsqp = 0
        dataw = ""
        for block in seqs:
            blockSize = len(block)
            blockc = block.strip().split("\n")
            seqId = blockc.pop(0)
            seqIdLen = len(seqId)
            seqLineLength = len(blockc[0])
            seqLineLengthH = seqLineLength+1  
            nsqp += (1+seqIdLen+2)
            seq = ''.join(blockc)
            dataw += f"{seqId}\t{len(seq)}\t{nsqp}\t{seqLineLength}\t{seqLineLengthH}\n"
            nsqp -= (seqIdLen+2)
            nsqp += blockSize
        dataw = dataw.strip()
        fhw.write(dataw)
        fhw.close()
        fh.close()

    # Read index from 2bit file
    def read2bitIndex(self):
        fh = tbit.open(self.filename)
        sequenceIndex = fh.chroms()
        fh.close()
        return sequenceIndex

    # Read sequences
    def readSequences(self,seqid):
        sequence = None
        if self.sformat == 'b':
            fh = tbit.open(self.filename)
            sequence = fh.sequence(seqid)
            fh.close()
        else:
            fh = open(self.filename,'r')
            sequenceParams = self.index[seqid]
            seqLength = sequenceParams[0]
            startpointer = sequenceParams[1]-1
            endpointer = startpointer+seqLength-1
            fh.seek(startpointer)
            total_lines = int(sequenceParams[0]/sequenceParams[2])
            total_words = total_lines*sequenceParams[3]
            sequence = fh.read(total_words)
            sequence = ''.join(sequence).replace('\n','')
            fh.close()

        return sequence.upper()
            
    # Convert sequence into concensus matrix
    def makeConcensus(self,sequence):
        seqarray = np.array(list(sequence))
        nmatrix = np.array([np.where(seqarray=='C',1,0),np.where(seqarray=='G',1,0),np.zeros(len(sequence),dtype=int)])
        
        ## new section for faster marking of the CpG presence
        # contains extra steps of copying 2 arrays and then deletion,appending and comparing
        copy_c = nmatrix[0].copy()
        copy_g = nmatrix[1].copy()
        copy_g = np.delete(copy_g,0)  
        
        # deletion step can be removed, and operation for array copy_g should start from index 1 
        copy_g = np.append(copy_g,0)
        cg = ((np.where(copy_c != 0,1,0) == 1) & (np.where(copy_g != 0,1,0) == 1)).astype(int)

        ## old method for marking CG dinucleotide: slower than the previous method, effect is visible on larger chromosomes
        # contains only 0ne step of comapring by through loop
        # cg = [1 if nmatrix[0][i] == 1 and nmatrix[1][i+1] == 1 else 0 for i in range(0,len(nmatrix[0])-1)]
        # cg.append(0)
        nmatrix[2] = cg
        consensus = np.vstack([list(np.cumsum(row)) for row in nmatrix])
        return consensus


class CountFeatures:
    def __init__(self,bedfile):
        self.bedFilename = bedfile

        # Validate the existence of the file and content
        self.validate()
        
        self.bedTable = self.readBEDFile()

    def validate(self):
        if not os.path.exists(self.bedFilename):
            print(f"Bed file do not exists")
            exit()
    
    def readBEDFile(self,skiplines=False):
        skip_rows = 0
        header = None

        # code for calculation the number of lines to be skipped
        if skiplines:
            fh = open(self.bedFilename,'r')
            for l in fh:
                if l[0] == '#':
                    skip_rows+=1
                    header = l.replace('#','')
                else:
                    break
            fh.close()

        if header != None:
            header=header.strip().split('\t')

        # Reading content into dataFrame
        bedTable = pd.read_csv(self.bedFilename,sep='\t',header=header,skiprows=skip_rows)
        return bedTable            


class CpGCalculator:
    def validate(self,chrid,index):
        if chrid not in index.keys():
            print(f"{chrid} id do not exixt in the index of provided sequence file")
            return False
        return True 

    # Calculate cpg ratio of the complete sequence filtered by feature
    def calByFeature(self,data,cons,chrid,feature):
        data = data.loc[(data[0] == chrid) & (data[2] == feature) & (data[6] == '+')]
        data = data.sort_values(by=[3])
        start = data[3]
        end = data[4]
        length = end-start+1
        cp = ((cons[0][end-1] - cons[0][start-1])+1)/length
        gp = ((cons[1][end-1] - cons[1][start-1])+1)/length
        cpgp = ((cons[2][end-1] - cons[2][start-1])+1)/length
        cpgRatio = round((cp*gp)/cpgp,2)
        outputdata = pd.DataFrame({'chrid': data[0],'feature':data[2],'start':start,'stop':end,'CpG_ratio':cpgRatio})
        return outputdata
        

class Visuals:
    def plotDensityByFeature(self):
        pass

    def plotHeatMapByFeature(self):
        pass



##<< Main body >> 
# main body for testing 2bit sequence file functions 
# ssf = SequenceFetcher("/home/grids3/Desktop/emprorcpg/dm6.2bit",'b')

ssf = SequenceFetcher("/home/grids3/Desktop/emprorcpg/dm6.fa",'f')
ccf = CountFeatures("/home/grids3/Desktop/emprorcpg/dm6.ncbiRefSeq.gtf")
seq = ssf.readSequences('chrM')
cons = ssf.makeConcensus(seq)

# Calculating CpG ratio for feature 'Transcript'
cpgcal = CpGCalculator()
output_data = cpgcal.calByFeature(ccf.bedTable,cons,'chrM','transcript')

# Calculating CpG ratio for all the segments in the genome
# for chrs in ssf.index.keys():  
#     print(chrs)
#     seq = ssf.readSequences(chrs)
#     cons = ssf.makeConcensus(seq)
#     #calculating CpG ratio for feature 'Transcript' for all chromosomes one by one
#     cpgcal = CpGCalculator()
#     output_data = cpgcal.calByFeature(ccf.bedTable,cons,chrs,'transcript')
#     print(output_data)

## logs
# Reading sequences from the file do not take much time : for the example file takes only ~2 seconds 
# Making Concensus matrix takes time: "no other solution"

## speed gain by logic
# 1. Many contig sequence do not have annotation or asked feature in the gtf file 
#    hence sequence reding and convertion for such contigs are of no use. A check 
#    is required which will lower the time comsumption

# 2. A class is required which can filter and then upload the content of the GTF file,
#    Comment lines in the middle of the file will result in reading GTF file into dataFrame. 
#    All the rows of the GTF file is no required hence avoiding that will give advance in memory usage.
