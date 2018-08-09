# Python
# $ pathway.py <Pathway Name>
# Pathway Name arg should be a folder in ../test containing fasta files of AA sequences for relevant proteins
# Finds the best matches for a collection of S. cerevisiae genes in the various databases

from subprocess import call
import sys, glob, requests, time, re
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import Phylo, AlignIO, SeqIO
from Bio.Phylo.Applications import PhymlCommandline

pathway = sys.argv[1]


def matrixLine(base):
    base = base.upper()
    if base == 'A':
        return (['12', '0', '0', '0'])
    elif base == 'C':
        return (['0', '12', '0', '0'])
    elif base == 'G':
        return (['0', '0', '12', '0'])
    elif base == 'T':
        return (['0', '0', '0', '12'])
    elif base == 'R':
        return (['6', '0', '6', '0'])
    elif base == 'Y':
        return (['0', '6', '0', '6'])
    elif base == 'S':
        return (['0', '6', '6', '0'])
    elif base == 'W':
        return (['6', '0', '0', '6'])
    elif base == 'K':
        return (['0', '0', '6', '6'])
    elif base == 'M':
        return (['6', '6', '0', '0'])
    elif base == 'B':
        return (['0', '4', '4', '4'])
    elif base == 'D':
        return (['4', '0', '4', '4'])
    elif base == 'H':
        return (['4', '4', '0', '4'])
    elif base == 'V':
        return (['4', '4', '4', '0'])
    elif base == 'N':
        return (['3', '3', '3', '3'])
    else:
        print('Invalid Character: "' + base + '"!')
        print('Please only use parentheticals thoughtfully')
        print('You may get erroneous results')
        quit()

entry = ''
while entry == '':
    print('Please enter a consensus sequence in IUPAC nucleotide notation')
    print('Or enter a single "H" for a description of this notation.')
    entry = input('Consensus: ')
    if entry == 'H' or entry == 'h':
        print('---')
        print('IUPAC Nucleotide Codes:')
        print('A,C,G,T - specify an unambiguous nucleotide')
        print('R - A or G')
        print('Y - C or T')
        print('S - G or C')
        print('W - A or T')
        print('K - G or T')
        print('M - A or C')
        print('B - C, G, or T')
        print('D - A, G, or T')
        print('H - A, C, or T')
        print('V - A, C, or G')
        print('N - Any base')
        print('Do not include gap characters (- or .)')
        print('To indicate an ambiguous region use e.g. N(2,3)')
        print('This corresponds to NN or NNN, the program will output multiple files if needed')
        print('Note that only one such variable region is supported')
        print('---')
        entry = ''
    else:
        entry = entry.strip('-')
        entry = entry.upper()
        if entry.strip('ACGTRYSWKMBHVN(,)1234567890') != '':
            print('---')
            print('Invalid characters: ' + entry.strip('ACGTRYSWKMBHVN(,)1234567890'))
            print('Please only use IUPAC notation without gaps (- or .)')
            print('---')
            entry = ''
            continue
        if '(' in entry and ',' in entry and ')' in entry:
            match = re.search(r'\w\(.+\)', entry)
            repeat = match[0]
            base = repeat[0]
            start = re.search(r'\((\d+),', match[0])
            start = start[1]
            stop = re.search(r',(\d+)\)', match[0])
            stop = stop[1]
            # reps = repeat[-2] - repeat[2] + 1
            for i in range(int(start), int(stop) + 1):
                newEntry = entry[:match.start()] + base * i + entry[match.end():]
                outfile = open('output/' + pathway + '/seqMatrix_'+ base + str(i) + '.txt', 'w')
                for i in range(len(newEntry)):
                    outfile.write(' '.join(matrixLine(newEntry[i])) + '\n')
                outfile.close()
        else:
            outfile = open('output/' + pathway + '/seqMatrix.txt', 'w')
            for i in range(len(entry)):
                outfile.write(' '.join(matrixLine(entry[i])) + '\n')
            outfile.close()

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
            if seqid in lines[i]:
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


#####################################################
# Identify S. cerevisiae genes in reference genomes #
#####################################################
enzymes = [] # [[Enzyme ID X, Fungus Y genome, Enzyme X gene ID in Y, Enzyme name],...]
indexfile = open('../test/' + pathway + '/index.csv', 'r')
accIndex = {}
enzymeIndex = {}
for line in indexfile.readlines():
    sline = line.split(',')
    accIndex[sline[0]] = sline[1].rstrip()
for enzymeId in glob.iglob('../test/' + pathway + '/*.fasta'):
    '../test/fat/A81283.fasta'
    acc = enzymeId.split('/')[-1].split('.')[0]
    hits = blasterMaster(enzymeId)
    for genome, seqid in hits:
         enzymes.append([enzymeId.split('/')[-1].split('.')[0], genome, seqid, accIndex[acc]])
         enzymeIndex[seqid] = accIndex[acc]
    # print (enzymes)

#########################################
# Collect 1 kbp upstream of those genes #
#########################################

# ids = [z for x, y, z in enzymes]
seqFile = open('output/' + pathway + '/upstream.fa', 'w')
controlFile = open('output/' + pathway + '/control.fa', 'w')
for x, genome, ide, fullName in enzymes:
    upper = genome[:-11] + 'upstream.fasta'
    gnome = open(upper, 'r')
    seq = fastaFinder(gnome.readlines(), ide)
    if genome == '../db/a_n_genes.fasta':
        controlFile.write('>' + ide + '\n')
        controlFile.write(seq + '\n')
    else:
        seqFile.write('>' + ide + ' | ' + fullName + '\n')
        seqFile.write(seq + '\n')
seqFile.close()
controlFile.close()

#######################################
# Do the above two steps for Blumeria #
#######################################

bluEnzymes = []
for enzymeId in glob.iglob('../test/' + pathway + '/*.fasta'):
    acc = enzymeId.split('/')[-1].split('.')[0]
    hit = blasterMaster(enzymeId, blu=True)
    enzymeIndex[hit] = accIndex[acc]
    bluEnzymes.append([enzymeId.split('/')[-1].split('.')[0], hit, accIndex[acc]])
print (bluEnzymes)

seqFile = open('output/' + pathway + '/bluUpstream.fa', 'w')
for x, ide, fullName in bluEnzymes:
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
    seqFile.write('>' + ide + ' | ' + fullName + '\n')
    seqFile.write(seq + '\n')
    gnome.close()
seqFile.close()


print (enzymeIndex)
################################################
# Identify motifs in those upstreams with MEME #
################################################

print ('Running MEME...')
matrices = [matrix for matrix in glob.glob('output/' + pathway + '/seqMatrix*')]
locs = ['output/' + pathway + '/fimo/blumeria/fimo.tsv','output/' + pathway + '/fimo/control/fimo.tsv','output/' + pathway + '/fimo/training/fimo.tsv']
organisms = {'Af':'A. fumigatus','AN':'A. nidulans','SS':'S. sclerotiorum','NC':'N. crassa','MG':'M. oryzae','BC':'B. cinerea','Bg':'B. graminis'}
results = [] # [[Seq ID, Enzyme name, strand, p-value, matched seq],[Another motif result], ...]
print (matrices)
for matrix in matrices:
    call('cat ' + matrix + ' | ~/meme/libexec/meme-5.0.1/matrix2meme | ~/meme/bin/fimo -oc output/' + pathway + '/fimo/control/ -thresh 0.001 - output/' + pathway + '/control.fa', shell=True)
    call('cat ' + matrix + ' | ~/meme/libexec/meme-5.0.1/matrix2meme | ~/meme/bin/fimo -oc output/' + pathway + '/fimo/training/ -thresh 0.001 - output/' + pathway + '/upstream.fa', shell=True)
    call('cat ' + matrix + ' | ~/meme/libexec/meme-5.0.1/matrix2meme | ~/meme/bin/fimo -oc output/' + pathway + '/fimo/blumeria/ -thresh 0.001 - output/' + pathway + '/bluUpstream.fa', shell=True)
    for loc in locs:
        fimoFile = open(loc, 'r')
        for line in fimoFile.readlines():
            sline = line.split('\t')
            if sline[0] == 'motif_id':
                continue # Skip first line
            if len(line) < 4:
                break # Stop before the end
            ide = sline[2]
            strand = sline[5]
            pval = sline[7]
            match = sline[-1].rstrip()
            results.append([ide, enzymeIndex[ide], strand, pval, match])
        fimoFile.close()

resultFile = open('output/' + pathway + '/results.tsv', 'w')
resultFile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format('Organism','Seq ID','Enzyme','Strand','p-value','Match Seq'))
for ide, name, strand, pval, match in results:
    resultFile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(organisms[ide[:2]],ide,name,strand,pval,match))
resultFile.close()

# results = [] # [[Seq ID, Enzyme name, strand, p-value, matched seq],[Another motif result], ...]
# for loc in locs:
#     fimoFile = open(loc, 'r')
#     for line in fimoFile.readlines():
#         sline = line.split('\t')
#         if sline[0] == 'motif_id':
#             continue # Skip first line
#         if len(line) < 4:
#             break # Stop before the end
#         ide = sline[2]
#         strand = sline[5]
#         pval = sline[7]
#         match = sline[-1].rstrip()
#         results.append([ide, enzymeIndex[ide], strand, pval, match])
#     fimoFile.close()
#
# resultFile = open('output/' + pathway + '/results.tsv', 'w')
# resultFile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format('Organism','Seq ID','Enzyme','Strand','p-value','Match Seq'))
# for ide, name, strand, pval, match in results:
#     resultFile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(organisms[ide[:2]],ide,name,strand,pval,match))
# resultFile.close()

quit()
# Currently looks for 5 motifs to make sure that the relevant one is found
# call('~/meme/bin/meme -oc output/' + pathway + '/meme output/' + pathway + '/upstream.fa -dna -mod oops -revcomp -objfun se -cons TCSNNNNNNNNSGA -cons GCMNNNNNNNNKGC', shell=True)

# memeOut = open('output/' + pathway + '/meme/meme.txt', 'r')
# motifs = []
# for line in memeOut.readlines():
#     if len(line) < 6:
#         continue
#     if line[:5] == 'MOTIF':
#         sline = line.split()
#         motifs.append(sline[1])
# print (len(motifs))
# useMotif = ''
# while useMotif == '':
#     for motif in motifs:
#         answer = input('Use this motif: ' + motif + '? (y/n)')
#         if answer == 'y':
#             useMotif = motif
#             break

#######################################################
# Look for those motifs in the same genes of Blumeria #
#######################################################

# motif='MEME-1'
# print ('Running FIMO...')
# # Run on the enzymes of one control species (A. nidulans)
# call('~/meme/bin/fimo -oc output/' + pathway + '/fimo/control/ -thresh 0.0005 output/' + pathway + '/meme/meme.xml output/' + pathway + '/control.fa', shell=True)
# # Run on the training data (to get a nice TSV of motif results)
# call('~/meme/bin/fimo -oc output/' + pathway + '/fimo/training/ -thresh 0.0005 output/' + pathway + '/meme/meme.xml output/' + pathway + '/upstream.fa', shell=True)
# # Run on Blumeria
# call('~/meme/bin/fimo -oc output/' + pathway + '/fimo/blumeria/ -thresh 0.0005 output/' + pathway + '/meme/meme.xml output/' + pathway + '/bluUpstream.fa', shell=True)

##############################################
# Produce a list of enzymes and their motifs #
##############################################

locs = ['output/' + pathway + '/fimo/blumeria/fimo.tsv','output/' + pathway + '/fimo/control/fimo.tsv','output/' + pathway + '/fimo/training/fimo.tsv']
organisms = {'Af':'A. fumigatus','AN':'A. nidulans','SS':'S. sclerotiorum','NC':'N. crassa','MG':'M. oryzae','BC':'B. cinerea','Bg':'B. graminis'}
results = [] # [[Seq ID, Enzyme name, strand, p-value, matched seq],[Another motif result], ...]
for loc in locs:
    fimoFile = open(loc, 'r')
    for line in fimoFile.readlines():
        sline = line.split('\t')
        if sline[0] == 'motif_id':
            continue # Skip first line
        if len(line) < 4:
            break # Stop before the end
        ide = sline[2]
        strand = sline[5]
        pval = sline[7]
        match = sline[-1].rstrip()
        results.append([ide, enzymeIndex[ide], strand, pval, match])
    fimoFile.close()

resultFile = open('output/' + pathway + '/results.tsv', 'w')
resultFile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format('Organism','Seq ID','Enzyme','Strand','p-value','Match Seq'))
for ide, name, strand, pval, match in results:
    resultFile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(organisms[ide[:2]],ide,name,strand,pval,match))
resultFile.close()

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
