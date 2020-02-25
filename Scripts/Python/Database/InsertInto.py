import sqlite3
import csv
from os import listdir
from os.path import isfile, join

def findAttributes(file):
    d = {}

    with open(file, 'r') as doc:
        doc = csv.reader(doc, delimiter = ' ')
        for row in doc:
            d[row[0]] = row[1]
    return d

def coerce(val, d, key):
    if str(val).lower() in ['na', '!value', '#value!']:
        return None
    if d[key] == 'INT':
        return int(val)
    elif d[key] == 'REAL':
        return float(val)
    else:
        return val

def insert(file, d):
    path = './cancer'

    connection = sqlite3.connect(path)

    c = connection.cursor()

    command = 'INSERT INTO patient( '

    nrow = dict.fromkeys(d)

    with open(file, 'r') as doc:
        doc = csv.reader(doc, delimiter='\t')
        m = 0
        for row in doc:
            if m == 0:
                m += 1
                continue
            for key, column in zip(nrow.keys(), row):
                nrow[key] = coerce(column, d, key)

            s = ''
            v = ''
            i = 0

            for key, value in nrow.items():
                i += 1
                if value is None:
                    continue
                if i == len(nrow)-1:
                    s += key
                    v += ':' + key
                    break
                s += key + ', '
                v += ':' + key + ', '

            ncommand = command + s + ') VALUES(' + v + ');'

            c.execute(ncommand, nrow)

    connection.commit()
    connection.close()

path = './Data/Attributes.txt'

#files = [f for f in listdir(path) if isfile( join(path, f) ) ]

#normal = []

#for f in files:
#    if 'normal' in f.lower():
#        normal.append( join(path, f) )

d = findAttributes( path )

path = './Data/Clinical Data/All.tsv'

insert(path, d)


