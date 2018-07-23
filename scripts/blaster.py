# Python
# $ blaster.py
# BLASTs genes against various fungal databases and returns the best result(s)

from subprocess import call
import glob

query = '../test/A6ZVQ9.fasta'
searched_genomes = []
for genome in glob.iglob('../db/*genes.fasta'):
    searched_genomes.append(genome)
    print (genome.split('/')[-1][0:4] + query.split('/')[-1][:-3])
    # tblastx search for the best hits among the genomes present in /db
    call('tblastn -db ../db/' + genome + ' -query ' + query + ' -out output/fat/' + genome.split('/')[-1][0:4] + query.split('/')[-1][:-3] + '.blast' + ' -max_target_seqs 1 -outfmt "7 sseqid evalue"', shell=True)
    #-evalue 0.0001 Shouldn't matter because we're only going to look at the best hit anyway

matches = []

for output in glob.iglob('output/fat/*.blast'):
    infile = open(output, 'r')
    i = 0
    for line in infile.readlines():
        i += 1
        if i==1 or i==2 or i==4 or i==5:
            continue
        if line[0] == '#':
            genome = line.split('/')[-1].rstrip()
            continue
        seqid, e_val = line.split()
        if float(e_val) < float(1e-04):
            matches.append(['../db/' + genome, seqid])
            break
    infile.close()

print (matches)
