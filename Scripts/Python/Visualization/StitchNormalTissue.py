"""
Stitches all normal tissue into one using only the given immune genes
"""
import csv
import sqlite3
from os import listdir
from os.path import isfile, join

#Variables
Src = "./Data/Transposed Gene Expression"
Destination = "./Data/Normal Immune Genes/NormalTissue_Combined.tsv"
Genes = "./Data/Trajanoski_CellRep_immune_signatures.txt"
Label = "Normal"
ImmuneGenes = {}

#Adds Immune Gene names in ImmuneGene
with open(Genes, 'r') as doc:
    doc = csv.reader(doc, delimiter = "\t")
    ImmuneGenes['barcode'] = None
    for row in doc:
        geneName = row[0].strip()
        geneNumber = row[len(row)-1].strip()
        ImmuneGenes[geneName + "|" + geneNumber] = None
    ImmuneGenes.pop("Metagene|EntrezGene.ID")

def coerce(val):
    if str(val).lower() in ['na', '!value', '#value!']:
        return None
    return float(val)

files = [join(Src, f) for f in listdir(Src) if isfile(join(Src, f)) and Label in f]

#Open connection
path = './cancer'
connection = sqlite3.connect(path)
c = connection.cursor()
#Get all distinct barcode's
command = "Select distinct barcode from n_expression"
c.execute(command)
barcodes = c.fetchall()
#Print ImmuneGenes
print(ImmuneGenes)
with open('./Data/Normal Immune Genes/NormalTissue_Combined.tsv', 'a') as new:
    new = csv.writer(new, delimiter = '\t')
    #Add header
    new.writerow(list(ImmuneGenes.keys()))
    #Add rest of rows
    for barcode in barcodes:
        #Get details for barcode
        command = 'select * from n_expression where barcode=?'
        c.execute(command, barcode)
        barcode = barcode[0].strip()
        ImmuneGenes['barcode'] = barcode
        for row in c.fetchall():
            gene_id = row[2]
            val = coerce(row[3])
            ImmuneGenes[gene_id] = val
        new.writerow(list(ImmuneGenes.values()))
        #Clean ImmuneGenes
        for key in ImmuneGenes.keys():
            ImmuneGenes[key] = None