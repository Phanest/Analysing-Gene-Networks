import sqlite3
import csv
from os import listdir
from os.path import isfile, join

#filename indicates whether we want to insert cancer or normal tissue information into database
#path folder in which files are present

def insert(filename, disease, table, expression):

    path = './cancer'

    connection = sqlite3.connect(path)

    c = connection.cursor()

    d = {}

    command = 'INSERT INTO ' + table + '(barcode, disease) VALUES (:barcode, :disease);'
    acommand = 'INSERT INTO ' + expression + '(barcode, disease, gene, value) VALUES (:barcode, :disease, :gene, :value);'

    with open(filename, 'r') as doc, open('./Data/' + table + 'Er.txt', 'a')\
            as writer:

        doc = csv.reader(doc, delimiter = '\t')

        m = 0

        for row in doc:
            if len(row) == 0:
                continue
            if m == 0:
                for column in row:
                    if column == 'gene_id':
                        d['barcode'] = None
                    else:
                        d[ column ] = None
            if m != 0:
                n = 0
                for key, column in zip( d.keys(), row):
                    if n == 0:
                        d['barcode'] = column
                        print(d['barcode'])
                        n += 1
                    else:
                        d[key] = coerce(column)
                        try:
                            c.execute(acommand, {'barcode': d['barcode'], 'disease': disease, 'gene': key, 'value': d[key]})
                        except:
                            writer.write(expression + ' barcode: ' + d['barcode'] + ' disease ' + disease + ' gene ' +
                                         key + ' value ' + str( d[key] ) )
                try:
                    c.execute(command, {'barcode': d['barcode'], 'disease': disease})
                except:
                    writer.write(table + ' barcode ' + d['barcode'] + ' disease ' + disease)

            m += 1

    connection.commit()
    connection.close()

def coerce(val):
    if str(val).lower() in ['na', '!value', '#value!']:
        return None
    return float(val)

filename = 'normal'
expression = 'n_expression'
path = './Data/Transposed Gene Expression'

files = [f for f in listdir(path) if isfile( join(path, f) ) and filename in f.lower() ]

filename = 'normal'

for f in files:
    parts = f.split('_')
    insert( join(path, f), parts[0], filename, expression)