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
#import json
import pickle

parser = argparse.ArgumentParser(description='Load all the data files (fasta format), construct context tree model, and store it (pickle format)')
parser.add_argument('-i', nargs='+', required=True,
                    help='input filenames (fasta)')
parser.add_argument('-o', nargs=1, required=True,
                    help='output filename (storing the context tree model)')
parser.add_argument('-t', type=str ,default = ['no','description','given'],nargs='+', required=False,
                    help='description of the performed experiment',dest='description')
parser.add_argument('-d',nargs=1,type=int,required=True,dest='depth',
                    help='depth of the tree')

args = parser.parse_args()
filenamesin = args.i
filenameout = args.o[0]
depth = args.depth[0]

startedtaskat = "{0}".format(str(datetime.datetime.now()))
print('-Performing construction task: '+' '.join(args.description)+'-')
print(' started at: '+startedtaskat)

print(' load data and count symbols')
tree = fulltree(depth)
for filename in filenamesin:
    print(' loading: {0}'.format(filename))
    handle = open(filename, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        sequences = record.seq.split('N')
        for sequence in sequences:
            tree.updatesymbolcounts(sequence)
    handle.close()

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
    


