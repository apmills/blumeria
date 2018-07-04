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
from Bio import Phylo
from Bio.Phylo.Applications import PhymlCommandline
from Bio import AlignIO

query = '../test/rnap.fa'
bluGenome = '../db/blumeria/latest/Bgt_genome_v2_1.fa'
# '../test/facb.fa'

#############################
# Check for some parameters #
#############################
# Loading modules on the msc server
if '-modules' in sys.argv:
    call('module load clustal', shell=True)
    call('module load hmmer', shell=True)
    call('module load phyml/3.1', shell=True)

# Allow user to repeat the same search (mostly for testing)
if '-skipBlast' in sys.argv: #or '-skipDomain' in sys.argv:
    repeatFile = open('output/repeat.txt', 'r')
    repeatData = repeatFile.readlines()
    repeatFile.close()
else:
    newRepeat = open('output/repeat.txt', 'w')

# newRepeat = open('output/repeat.txt', 'w')

#############################################
# BLAST query sequences (genes of interest) #
#############################################
def blasterMaster(query):
    # This function handles the BLAST search (through the CLI) and the results
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
                # Skip over some irrelevant lines output by BLAST
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
    # Use data from output/repeat.txt instead of running BLAST (to save time while testing)
    print ('Skipping BLAST search')
    exec('hits = ' + repeatData[0].rstrip())
else:
    print('BLASTing...')
    hits = blasterMaster(query)
    newRepeat.write(str(hits) + '\n')

#################################################
# Retrieve sequences of best hit in each genome #
#################################################
def fastaFinder(lines, seqid):
    # Collects all of the sequence after a desired FASTA header
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
    # Not really a necessary function, but it calls fastaFinder() on all the items in a list of lists
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

#######################
# Translate sequences #
#######################

def domainSpotter(seq):
    # Searches the NCBI's online CDD for any conserved domains within a given ORF

    r = requests.post('https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?', data=(('db','cdd'), ('queries',seq), ('tdata','hits')))

    # print (r.text)

    searchID = r.text.split('#')[2].split()[1].rstrip()

    # print (searchID)

    tick = True
    done = False
    print ('Searching for domains...', end='', flush=True)

    while not done:
        # Ping the website every 2 seconds to see if it's done
        result = requests.post('https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?', data={'cdsid':searchID, 'cddefl':'true'})
        status = result.text.split('#')[4].split()[1].rstrip()
        if status == '3':
            # A 3 here indicates that the search is still in progress
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
    # Uses data from output/repeat.txt instead of prompting the user to select which ORFs he wants
    print ('Skipping Domain search')

    #exec('protSeqs = ' + repeatData[1])
else:
    protSeqs = []
    for genome, seqid, seq in nucSeqs:
        # Using the Seq class from Biopython to perform the translation
        bioSeq = Seq(seq, generic_dna)
        protSeq = str(bioSeq.translate())
        for orf in protSeq.split('*'):
            # Splitting on stop codons to iterate through each ORF
            if len(orf) <= 50:
                print ('Skipping short ORF (<50 residues)')
                continue

            tags = domainSpotter(orf)
            if len(tags) == 0:
                print ('No tags found, skipping ORF')
                continue
            # Give the user some information and allow him to decide if this ORF looks acceptable
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
#######################################################################
# At this point it should really try to find domains in different RFs #
#######################################################################

###################
# Align sequences #
###################
print ('Aligning sequences...')
seqFile = open('output/unaligned.fa', 'w')
# for genome, seqid, seq in protSeqs:
#     seqFile.write('>' + seqid + ' ORF | \n')      Protein sequence
#     seqFile.write(seq + '\n')

for genome, seqid, seq in nucSeqs:
    seqFile.write('>' + seqid + ' Whole gene | \n')      #Nucleotide sequence
    seqFile.write(seq + '\n')
seqFile.close()

# Use Clustal Omega to align the sequences
call('clustalo -i output/unaligned.fa -o output/aligned.aln --force --outfmt=clu', shell=True)

#######################################
# Look for query sequence in Blumeria #
#######################################
place = 1
thisContig = False
def genomeSlicer(start, end, line, contig):
    # Call this on iterating lines of a genome file and it will return all of the sequence that falls between start and end in the relevant contig/sequence/FASTA header
    # Global variables to persist through multiple calls to the function
    global place
    global thisContig
    if thisContig and place >= end:
        # Cut off the function if we've got all the desired sequence
        thisContig = False
        return
    if line[0] == '>':
        # Check if this FASTA header is the one we want
        if contig in line:
            thisContig = True
        return
    if thisContig:
        # Each line has 50 characters on it, only initialise the list if we're going to need it here
        if place > start - 50:
            output = []
        else:
            place += 50
            return
        for nuc in line.rstrip():
            # Grab every nucleotide that falls between start and end
            if place >= start and place < end:
                output.append(nuc)
            place += 1
        return ''.join(output)
    else:
        return

# First create an HMM profile based on the alignment
print ('Building HMM profile...')
call('hmmbuild -o /dev/null output/query.hmm output/aligned.aln', shell=True)

# Now search it (with nhmmer for nucleotides) against the whole Blumeria genome
# call('hmmsearch --tblout output/hmmertbl.tsv -o output/hmmer.txt output/query.hmm ' + bluGenome, shell=True)
print ('Running HMMER search...')
call('nhmmer --tblout output/hmmertbl.tsv -o output/nhmmer.txt output/query.hmm ' + bluGenome, shell=True)

# The result is saved in this tabular file, we're interested in only the best hit
result = open('output/hmmertbl.tsv', 'r')
lines = result.readlines()
result.close()
if lines[2][0] == '#':
    # If there are no hits, then let the user know and don't run the follow up analysis
    print ('No equivalent found in Blumeria graminis')
    quit()

# Extract the location of the hit within the genome FASTA file
start = int(lines[2].split()[6])
stop = int(lines[2].split()[7])
contig = lines[2].split()[0]
print (start)
print (stop)
bluFile = open('output/b_g_GOI.fa', 'w')
genomeFile = open(bluGenome, 'r')
# bluSeq = seqRetriever([['../db/blumeria/b_g_genes.fasta', bluGene]])
# print (bluSeq)
print ('Extracting best match from Blumeria genome at ' + bluGenome + '...')
output = [genomeSlicer(start, stop, x, contig) for x in genomeFile]
selection = [x for x in output if x is not None]
bluSeq = ''.join(selection)
bluFile.write('>Blumeria graminis candidate gene' + '\n')
bluFile.write(bluSeq)
genomeFile.close()
bluFile.close()

# And counter-search to check if it's really a good match, the top hit should be the query sequence
call('tblastx -db ../db/a_n_genes.fasta -query output/b_g_GOI.fa -out output/counterBlast.blast -num_threads 4 -max_target_seqs 1 -outfmt "7 sseqid evalue"', shell=True)

################################
# Identify motifs in sequences #
################################

def ORFanarium(inSeq):
    # Uses advanced AI (if statements) to select the best ORF in a given sequence
    readingFrames = [Seq.translate(inSeq.seq[i:], table='Standard', stop_symbol='*', to_stop=False, cds=False) for i in range(3)]


#################################
# Make a phylogeny of sequences #
#################################
unal = open('output/unaligned.fa', 'a')
bluFile = open('output/b_g_GOI.fa', 'r')
addition = bluFile.read()
unal.write(addition)
unal.close()
bluFile.close()
call('clustalo -i output/unaligned.fa -o output/alignedAll.aln --force --outfmt=clu', shell=True)
AlignIO.convert('output/alignedAll.aln', 'clustal', 'output/phyAlign.phy', 'phylip-relaxed')
cmdline = PhymlCommandline(input='output/phyAlign.phy', alpha='e', bootstrap=100, sequential=False)
call(str(cmdline), shell=True)
my_tree = Phylo.read("output/phyAlign.phy_phyml_tree.txt", "newick")
Phylo.draw(my_tree, show_confidence=True)

# Got to print 'Done' at the end, that's just common decency
print ('Done')
