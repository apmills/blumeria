# Python
# $ pathway.py <Pathway Name>
# Pathway Name arg should be a folder in ../test containing fasta files of AA sequences for relevant proteins
# Finds the best matches for a collection of S. cerevisiae genes in the various databases

from subprocess import call
import sys, glob, requests, time, re
# import glob
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
# import requests
# import time
from Bio import Phylo, AlignIO, SeqIO
from Bio.Phylo.Applications import PhymlCommandline
# from Bio import AlignIO
# from Bio import SeqIO
# import re

def blasterMaster(query, blu=False):
    # This function handles the BLAST search (through the CLI) and the results
    searched_genomes = []
    if blu:
        call('tblastn -db ../db/blumeria/latest/bluGenes.fa -query ' + query + ' -out output/' + pathway + '/blu/' + query.split('/')[-1][:-3] + '.blast' + ' -num_threads 4 -max_target_seqs 1 -outfmt "7 sseqid evalue"', shell=True)
        match = '' # Only one match for Blumeria
        infile = open('output/' + pathway + '/blu/' + query.split('/')[-1][:-3] + '.blast', 'r')
        i = 0
        for line in infile.readlines():
            i += 1
            if i==1 or i==2 or i==4 or i==5:
                # Skip over some irrelevant lines output by BLAST
                continue
            if line[0] == '#':
                genome = line.split('/')[-1].rstrip()
                continue
            seqid, e_val = line.split()
            if float(e_val) < float(1e-04):
                match = seqid
                break
            else:
                # input(query + ' not found in ' + genome + '!')    # Warn the user that the gene was missing in a reference genome
                break
        infile.close()
        return match
    else:
        for genome in glob.iglob('../db/*genes.fasta'):
            searched_genomes.append(genome)
            print (genome.split('/')[-1][0:4] + query.split('/')[-1][:-3])
            # tblastx search for the best hits among the genomes present in /db
            call('tblastn -db ' + genome + ' -query ' + query + ' -out output/' + pathway + '/' + genome.split('/')[-1][0:4] + query.split('/')[-1][:-3] + '.blast' + ' -num_threads 4 -max_target_seqs 1 -outfmt "7 sseqid evalue"', shell=True)
            #-evalue 0.0001 Shouldn't matter because we're only going to look at the best hit anyway

        matches = []
        #print ('Searched ' + ','.join(searched_genomes[:-1]) + ', and ' + searched_genomes[-1] + ' for ' + query)

        # Read through BLAST output and return the seqid of any hits with an e-value less than 1e-4
        for output in glob.iglob('output/' + pathway + '/*' + query.split('/')[-1][:-3] + '.blast'):
            infile = open(output, 'r')
            i = 0
            for line in infile.readlines():
                i += 1
                if i==1 or i==2 or i==4 or i==5:
                    # Skip over some irrelevant lines output by BLAST
                    continue
                if line[0] == '#':
                    genome = line.split('/')[-1].rstrip()
                    continue
                seqid, e_val = line.split()
                if float(e_val) < float(1e-04):
                    matches.append(['../db/' + genome, seqid])
                    break
                else:
                    # input(query + ' not found in ' + genome + '!')    # Warn the user that the gene was missing in a reference genome
                    break
            infile.close()

        #print (matches)
        return (matches)

def fastaFinder(lines, seqid):
    # Collects all of the sequence after a desired FASTA header
    i = 0
    seqs = []
    collect = False
    while i < len(lines):
        if lines[i][0] == '>':
            if re.search(seqid + r'\s', lines[i]):
                collect = True
                print('Found ' + seqid + ' in ' + lines[i])
            elif collect:
                break
        elif collect:
            seqs.append(lines[i].rstrip().upper())
        i += 1

    return ''.join(seqs)

def seqRetriever(seqinfo):
    # Not really a necessary function, but it calls fastaFinder() on all the items in a list of lists
    output = []
    for genome, seqid in seqinfo:
        print ('Retrieving ' + seqid + ' from ' + genome + '...')
        infile = open(genome, 'r')
        seq = fastaFinder(infile.readlines(), seqid)
        output.append([genome, seqid, seq])
    return output

pathway = sys.argv[1]

#####################################################
# Identify S. cerevisiae genes in reference genomes #
#####################################################
enzymes = [] # [[Enzyme ID X, Fungus Y genome, Enzyme X gene ID in Y],...]
for enzymeId in glob.iglob('../test/' + pathway + '/*'):
    hits = blasterMaster(enzymeId)
    for genome, seqid in hits:
         enzymes.append([enzymeId.split('/')[-1].split('.')[0], genome, seqid])
    # print (enzymes)

#########################################
# Collect 1 kbp upstream of those genes #
#########################################

# ids = [z for x, y, z in enzymes]
seqFile = open('output/' + pathway + '/upstream.fa', 'w')
controlFile = open('output/' + pathway + '/control.fa', 'w')
for x, genome, ide in enzymes:
    upper = genome[:-11] + 'upstream.fasta'
    gnome = open(upper, 'r')
    seq = fastaFinder(gnome.readlines(), ide)
    if genome == '../db/a_n_genes.fasta':
        controlFile.write('>' + ide + '\n')
        controlFile.write(seq + '\n')
    else:
        seqFile.write('>' + ide + '\n')
        seqFile.write(seq + '\n')
seqFile.close()
controlFile.close()

#######################################
# Do the above two steps for Blumeria #
#######################################

bluEnzymes = []
for enzymeId in glob.iglob('../test/' + pathway + '/*'):
    hit = blasterMaster(enzymeId, blu=True)
    bluEnzymes.append([enzymeId.split('/')[-1].split('.')[0], hit])
print (bluEnzymes)

seqFile = open('output/' + pathway + '/bluUpstream.fa', 'w')
for x, ide in bluEnzymes:
    upper = genome[:-11] + 'upstream.fasta'
    gnome = open('../db/blumeria/latest/bluUpstream0.fa', 'r')
    lines = gnome.readlines()
    i = 0
    seqs = []
    collect = False
    while i < len(lines):
        if lines[i][0] == '>':
            if lines[i].rstrip() == '>' + ide:
                collect = True
                print('Found ' + ide + ' in ' + lines[i])
            elif collect:
                break
        elif collect:
            seqs.append(lines[i].rstrip().upper())
        i += 1

    seq = ''.join(seqs)
    seqFile.write('>' + ide + '\n')
    seqFile.write(seq + '\n')
    gnome.close()
seqFile.close()

################################################
# Identify motifs in those upstreams with MEME #
################################################

print ('Running MEME...')
# Currently looks for 5 motifs to make sure that the relevant one is found
call('~/meme/bin/meme -oc output/' + pathway + '/meme output/' + pathway + '/upstream.fa -dna -mod oops -revcomp -nmotifs 3 -objfun se', shell=True)

memeOut = open('output/' + pathway + '/meme/meme.txt', 'r')
motifs = []
for line in memeOut.readlines():
    if len(line) < 6:
        continue
    if line[:5] == 'MOTIF':
        sline = line.split()
        motifs.append(sline[1])
print (len(motifs))
useMotif = ''
while useMotif == '':
    for motif in motifs:
        answer = input('Use this motif: ' + motif + '? (y/n)')
        if answer == 'y':
            useMotif = motif
            break

#######################################################
# Look for those motifs in the same genes of Blumeria #
#######################################################

print ('Running FIMO...')
call('~/meme/bin/fimo -oc output/' + pathway + '/fimo/control/ -thresh 0.001 -motif ' + motif + ' output/' + pathway + '/meme/meme.xml output/' + pathway + '/control.fa', shell=True)

call('~/meme/bin/fimo -oc output/' + pathway + '/fimo/blumeria/ -thresh 0.001 -motif ' + motif + ' output/' + pathway + '/meme/meme.xml output/' + pathway + '/bluUpstream.fa', shell=True)

quit()

# motif = Seq('CCTCGG')
# print (ids)
# allMotifs = [] # [[[Gene ID, # Motifs in upstream 1kb], [Another gene, #]], [[Different Organism]], ...]
# for upper in glob.iglob('../db/*upstream.fasta'):
#     motifs = []
#     infile = open(upper, 'r')
#     for record in SeqIO.parse(infile, 'fasta'):
#         if record.id[1:-1] not in ids:
#             continue
#         count = len(re.findall(str(motif), str(record.seq).upper()))
#         # print('>' + record.id[1:-1])
#         # print(str(record.seq))
#         count += len(re.findall(str(motif.reverse_complement()), str(record.seq).upper()))
#         if count > 0:
#             motifs.append([record.id[1:-1], count])
#         # if count > 1:
#         #     print(record.id + " " + str(count))
#     tot = 0
#     for ide, count in motifs:
#         tot += count
#     #print ('Found ' + str(tot) + ' occurences of "' + motif + '" in ' + upper)
#     #print (len(motifs))
#
#     allMotifs.append(motifs)
# print (allMotifs)



"""
- List of enzymes on pathway
- For each enzyme, list of 6 equivalents in other fungi
- For each equivalent, the upstream sequence
- For each upstream seq, how many motifs are in it
- Return # Motifs and enzyme name
"""
