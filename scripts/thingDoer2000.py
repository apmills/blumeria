import time
from Bio.Seq import Seq
tic = time.clock()
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

annotation = open('../db/blumeria/latest/Bgt_CDS_2014_20_clean.gff', 'r')
geneGenie = [] # Lives on his back
extract = [] # [[start, stop, contig, name, strand]]

for line in annotation:
    sline = line.split()
    if sline[0] == '###':
        continue
    if sline[2] == 'gene':
        extract.append([int(sline[3]), int(sline[4]), sline[0], sline[-1].split('=')[-1], sline[6]])

total = len(extract)
current = 0
genomeFile = open('../db/blumeria/latest/Bgt_genome_v2_1.fa', 'r')
genomeLines = genomeFile.readlines()
genomeFile.close()
for start, stop, contig, name, strand in extract:
    print(str(round(((current/total) * 100), 4)) + '% Done', end='\r')
    place = 1
    thisContig = False
    output = [genomeSlicer(start, stop, x, contig) for x in genomeLines]
    selection = [x for x in output if x is not None]
    bluSeq = ''.join(selection)
    if strand == '-':
        bluSeqC = Seq(bluSeq)
        bluSeqC = bluSeqC.reverse_complement()
        geneGenie.append([name, str(bluSeqC)])
    else:
        geneGenie.append([name, str(bluSeq)])
    current += 1

# print (upstreamGirl)
outfile = open('bluGenes.fa', 'w')
for name, seq in geneGenie:
    outfile.write('>' + name + '\n')
    outfile.write(seq + '\n')
outfile.close()
toc = time.clock()
print ('And it only took ' + str(toc - tic) + ' seconds')
