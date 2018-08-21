import time
import datetime
from Bio.Seq import Seq

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

annotation = open('../db/blumeria/latest/Bgh.EF2.40.gff3', 'r')
upstreamGirl = [] # Looking for downstream man [name, sequence]
extract = [] # [[start, stop, contig, name, strand]]

for line in annotation:
    sline = line.split()
    if sline[0] == '###':
        continue
    if sline[2] != 'gene':
        continue
    if sline[6] == '+':
        extract.append([int(sline[3])-1001, int(sline[3]) -1, sline[0], sline[-1].split(';')[0].split(':')[1], '+'])
    elif sline[6] == '-':
        extract.append([int(sline[4])+1, int(sline[4]) + 1001, sline[0], sline[-1].split(';')[0].split(':')[1], '-'])

tic = time.clock()
total = len(extract)
current = 0
genomeFile = open('../db/blumeria/latest/Bgh_genome.fa', 'r')
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
        upstreamGirl.append([name, str(bluSeqC)])
    else:
        upstreamGirl.append([name, bluSeq])
    current += 1
    toc = time.clock()
    print('Elapsed time: ' + str(datetime.timedelta(seconds=(toc-tic))) + '| Remaining: ' + str(datetime.timedelta(seconds=(float(total)/current)*(toc-tic)-(toc-tic))))

print (upstreamGirl)
outfile = open('bghUpstream.fa', 'w')
for name, seq in upstreamGirl:
    outfile.write('>' + name + '\n')
    outfile.write(seq + '\n')
outfile.close()
toc = time.clock()
print ('And it only took ' + str(toc - tic) + ' seconds')
