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
import datetime
#import json
import pickle

parser = argparse.ArgumentParser(description='Load the tree models (stored in pickle format) and compare them')
parser.add_argument('-i', nargs='+', required=True,
                    help='input filenames (pickle)')
parser.add_argument('-o', nargs=1, required=True,
                    help='output filename (storing a table of divergences)')
parser.add_argument('-t', type=str ,default = ['no','description','given'],nargs='+', required=False,
                    help='description of the performed experiment',dest='description')

args = parser.parse_args()
filenamesin = args.i
filenameout = args.o[0]

startedtaskat = "{0}".format(str(datetime.datetime.now()))
print('-Performing comparison task: '+' '.join(args.description)+'-')
print(' started at: '+startedtaskat)

print(' load the models and compare them one by one')
divergence = dict()
for filename1 in filenamesin:
    print(' loading: {0}'.format(filename1))
    with open(filename1, 'rb') as treefile:
        input1 = pickle.load(treefile)
    tree1 = input1['tree']
    
    del(input1)
    divs = dict()
    
    for filename2 in filenamesin:
        with open(filename2, 'rb') as treefile:
            input1 = pickle.load(treefile)
        
            tree2 = input1['tree']
        
        divs[filename2] = tree1.getdivergence(tree2)
        
    divergence[filename1] = divs

completedtaskat = "{0}".format(str(datetime.datetime.now()))
print(' completed at: '+completedtaskat)
print(' storing results in file: '+filenameout)

result = dict()
result['divergence'] = divergence
result['args'] = args
result['description'] = ' '.join(args.description)
result['starttask'] = startedtaskat
result['endtask'] = completedtaskat
with open(filenameout, 'wb') as output:
    pickle.dump(result, output)
    
print('-end-')


   
#with open(filenameout, 'w') as f:
#    json.dump(result,f)
    


