# Read all Tab-delimited rows from stdin.
from os import listdir
from os.path import isfile, join
import csv


def transpose( name, fro, to):
    with open(fro, 'r') as doc, open( join(to, name), 'w' ) as new:
        data = []
        transposed = []

        doc = csv.reader(doc, delimiter = '\t')
        new = csv.writer( new, delimiter = '\t' )

        for row in doc:
            data.append( row )

        for row in range(len(data[0])):
            new_row = []
            for col in range(len(data)):
                new_row.append(data[col][row])
            transposed.append(new_row)

        new.writerows( transposed )


path = './Data/Gene Expression'
copyTo = './Data/Transposed Gene Expression'

files = [f for f in listdir(path) if isfile( join(path, f) ) ]

for f in files:
    transpose( f, join(path, f), copyTo)