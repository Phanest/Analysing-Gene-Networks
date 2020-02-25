from os import listdir
from os.path import isfile, join
import csv

def findAttributes(files):
    d = {}

    for f in files:
        with open(f, 'r') as doc:
            doc = csv.reader(doc, delimiter = '\t')
            for row in doc:
                for column in row:
                    if column not in d:
                        d[column] = 0
                    d[column] += 1
                break
    return d

path = './Data/Transposed Gene Expression'

files = [f for f in listdir(path) if isfile( join(path, f) ) ]

normal = []
cancer = []

for f in files:
    print (join(path, f) )
    if 'normal' in f.lower():
        normal.append( join(path, f) )
    else:
        cancer.append( join(path, f) )

normalattributes = findAttributes( normal )
cancerattributes = findAttributes( cancer )

print ('Normal Tissue')
print ('min', min( normalattributes.values() ) )
print ('max', max( normalattributes.values() ) )
print ('length', len( normalattributes.values() ) )
print ('files', len( normal ) )
s = set( normalattributes.values() )
print (s)
print ('Cancer')
print ('min', min( cancerattributes.values() ) )
print ('max', max( cancerattributes.values() ) )
print ('length', len( cancerattributes.values() ) )
print ('files', len( cancer ) )
s = set( cancerattributes.values() )
print (s)