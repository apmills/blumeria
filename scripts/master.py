# Python
# $ master.py <-skipBlast> <-skipDomain> <-modules>
# skipBlast - skips the BLAST step, skipDomain - skips the domain checking step, modules - loads all required modules
# Master control script for the Blumeria comparative genomics pipeline

from subprocess import call
import sys
import glob
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import requests
import time

query = '../test/rnap.fa'
# '../test/facb.fa'

# Check for some parameters
if '-modules' in sys.argv:
    call('module load clustal', shell=True)
    call('module load hmmer', shell=True)

if '-skipBlast' in sys.argv or '-skipDomain' in sys.argv:
    repeatFile = open('output/repeat.txt', 'r')
    repeatData = repeatFile.readlines()
    repeatFile.close()
else:
    newRepeat = open('output/repeat.txt', 'w')

# BLAST query sequences (genes of interest)
def blasterMaster(query):
    searched_genomes = []
    for genome in glob.iglob('../db/*.fasta'):
        searched_genomes.append(genome)
        print (genome.split('/')[-1][0:4] + query.split('/')[-1][:-3])
        # tblastx search for the best hits among the genomes present in /db
        call('tblastx -db ' + genome + ' -query ' + query + ' -out output/' + genome.split('/')[-1][0:4] + query.split('/')[-1][:-3] + '.blast' + ' -num_threads 4 -max_target_seqs 1 -outfmt "7 sseqid evalue"', shell=True)
        #-evalue 0.0001 Shouldn't matter because we're only going to look at the best hit anyway

    matches = []
    print ('Searched ' + ','.join(searched_genomes[:-1]) + ', and ' + searched_genomes[-1] + ' for ' + query)

    # Read through BLAST output and return the seqid of any hits with an e-value less than 1e-4
    for output in glob.iglob('output/*' + query.split('/')[-1][:-3] + '.blast'):
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
    return (matches)

if '-skipBlast' in sys.argv:
    print ('Skipping BLAST search')
    exec('hits = ' + repeatData[0].rstrip())
else:
    print('BLASTing...')
    hits = blasterMaster(query)
    newRepeat.write(str(hits) + '\n')


# Retrieve sequences of best hit in each genome
def fastaFinder(lines, seqid):
    i = 0
    seqs = []
    collect = False
    while i < len(lines):
        if lines[i][0] == '>':
            if lines[i].find(seqid) >= 0:
                collect = True
                print('Found ' + seqid + ' in ' + lines[i])
            elif collect:
                break
        elif collect:
            seqs.append(lines[i].rstrip().upper())
        i += 1

    return ''.join(seqs)

def seqRetriever(seqinfo):
    output = []
    for genome, seqid in seqinfo:
        print ('Retrieving ' + seqid + ' from ' + genome + '...')
        infile = open(genome, 'r')
        seq = fastaFinder(infile.readlines(), seqid)
        output.append([genome, seqid, seq])
    return output

print (hits)
nucSeqs = seqRetriever(hits)

print (nucSeqs)

# Translate sequences

def domainSpotter(seq):
    r = requests.post('https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?', data=(('db','cdd'), ('queries',seq), ('tdata','hits')))

    # print (r.text)

    searchID = r.text.split('#')[2].split()[1].rstrip()

    # print (searchID)

    tick = True
    done = False
    print ('Searching for domains...', end='', flush=True)

    while not done:
        result = requests.post('https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?', data={'cdsid':searchID, 'cddefl':'true'})
        status = result.text.split('#')[4].split()[1].rstrip()
        if status == '3':
            done = False
        else:
            done = True
        time.sleep(2)
        print ('.', end = '', flush=True)
        # if tick:
        #     print ('Tick')
        #     tick = False
        # else:
        #     print ('Tock')
        #     tick = True

    lines = result.text.split('\n')

    tags = []
    for line in lines:
        # if line[0] == '#' or line[0] == '' or line[:5] == 'Query':
        #     continue
        # print (line)
        if line[:2] == 'Q#':
            tags.append(line.split('\t')[8])


    return tags



if '-skipDomain' in sys.argv:
    print ('Skipping Domain search')

    exec('protSeqs = ' + repeatData[1])
else:
    protSeqs = []
    for genome, seqid, seq in nucSeqs:
        bioSeq = Seq(seq, generic_dna)
        protSeq = str(bioSeq.translate())
        for orf in protSeq.split('*'):
            if len(orf) <= 50:
                print ('Skipping short ORF (<50 residues)')
                continue

            tags = domainSpotter(orf)
            if len(tags) == 0:
                print ('No tags found, skipping ORF')
                continue
            print('\nORF of length ' + str(len(orf)))
            print('With domain tags:\n-' + '\n-'.join(tags))
            if input('\nKeep this ORF? (y/n)  ') == 'y':
                protSeqs.append([genome, seqid, orf])
                break
            else:
                print('#' * 80)
        # if genome == '../db/m_o_genes.fasta':


        print (protSeqs)

    newRepeat.write(str(protSeqs) + '\n')
#########################################################################
#   At this point it should really try to find domains in different RFs #
#########################################################################

# Align sequences

print ('Aligning sequences...')
seqFile = open('output/unaligned.fa', 'w')
for genome, seqid, seq in protSeqs:
    seqFile.write('>' + seqid + ' ORF | \n')
    seqFile.write(seq + '\n')
seqFile.close()

call('clustalo -i output/unaligned.fa -o output/aligned.aln --force', shell=True)

# Look for query sequence in Blumeria

call('hmmbuild output/query.hmm output/aligned.aln', shell=True)

call('hmmsearch --tblout output/hmmertbl.tsv -o output/hmmer.txt output/query.hmm ../db/blumeria/b_g_pep.fa', shell=True)

# And counter-search to check if it's really a good match

result = open('output/hmmertbl.tsv', 'r')
lines = result.readlines()
result.close()
bluGene = lines[3].split()[0]
bluFile = open('output/b_g_GOI.fa', 'w')
bluSeq = seqRetriever([['../db/blumeria/b_g_genes.fasta', bluGene]])
print (bluSeq)
bluFile.write('>' + bluGene + '\n')
bluFile.write(bluSeq[0][2])
bluFile.close()

call('tblastx -db ../db/a_n_genes.fasta -query output/b_g_GOI.fa -out output/counterBlast.blast -num_threads 4 -max_target_seqs 1 -outfmt "7 sseqid evalue"', shell=True)


# Identify motifs in sequences


# Make a phylogeny of sequences


print ('Done')
