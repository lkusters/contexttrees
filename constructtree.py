# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 14:48:39 2016

@author: Lieneke Kusters
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 10:03:50 2016

@author: Lieneke Kusters

MAIN file:
Called from console
Take input
Call functions
"""

import argparse
from Bio import SeqIO
from contexttree import fulltree
import datetime
import gzip 
#import json
import pickle

parser = argparse.ArgumentParser(description='Load all the data files (zipped fasta/fastq format), construct context tree model, and store it (pickle format)')
parser.add_argument('-i', nargs='+', required=True,
                    help='input filenames (.fna.gz or .fastq.gz)')
parser.add_argument('-o', nargs=1, required=True,
                    help='output filename (storing the context tree model)')
parser.add_argument('-t', type=str ,default = ['no','description','given'],nargs='+', required=False,
                    help='description of the performed experiment',dest='description')
parser.add_argument('-d',nargs=1,type=int,required=True,dest='depth',
                    help='depth of the tree')
parser.add_argument('-r',nargs=1,type=str,required=True,dest='revcomp',
                    help='automatically also include reverse complement of each sequence? y/n')

args = parser.parse_args()
filenamesin = args.i
filenameout = args.o[0]
depth = args.depth[0]
if args.revcomp[0] == 'y':
    revcom = True
elif args.revcomp[0] == 'n':
    revcom = False
else:
    raise ValueError("not clear wether reverse complement should be included, choose -r y/n",args.revcomp[0])

startedtaskat = "{0}".format(str(datetime.datetime.now()))
print('-Performing construction task: '+' '.join(args.description)+'-')
print(' started at: '+startedtaskat)

print(' load data and count symbols')
tree = fulltree(depth)
for filename in filenamesin:
    print(' loading: {0}'.format(filename))
    handle = gzip.open(filename,'rt')
    #handle = open(filename, "rU")
    checkextension = filename.split('.')
    if checkextension[-2] == 'fna':
        for record in SeqIO.parse(handle, "fasta"):
            sequences = record.seq.split('N')
            for sequence in sequences:
                tree.updatesymbolcounts(str(sequence))
                if revcom:
                    tree.updatesymbolcounts(str(sequence.reverse_complement()))
        handle.close()
    elif checkextension[-2] == 'fastq':
        for record in SeqIO.parse(handle, "fastq"):
            sequences = record.seq.split('N')
            for sequence in sequences:
                tree.updatesymbolcounts(str(sequence))
                if revcom:
                    tree.updatesymbolcounts(str(sequence.reverse_complement()))
        handle.close()
    else:
        print("filename extension {0} not recognised".format(checkextension[-2]))

completedtaskat = "{0}".format(str(datetime.datetime.now()))
print(' completed at: '+completedtaskat)
print(' storing results in file: '+filenameout)

result = dict()
result['tree'] = tree
result['args'] = args
result['description'] = ' '.join(args.description)
result['starttask'] = startedtaskat
result['endtask'] = completedtaskat
with open(filenameout, 'wb') as output:
    pickle.dump(result, output)
    
print('-end-')


   
#with open(filenameout, 'w') as f:
#    json.dump(result,f)
    


