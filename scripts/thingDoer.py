import time
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
upstreamGirl = [] # Looking for downstream man [name, sequence]
extract = [] # [[start, stop, contig, name]]

for line in annotation:
    sline = line.split()
    if sline[0] == '###':
        continue
    if sline[2] == 'gene' and sline[6] == '+':
        extract.append([int(sline[3])-1001, int(sline[3]) -1, sline[0], sline[-1].split('=')[-1]])

total = len(extract)
current = 0
genomeFile = open('../db/blumeria/latest/Bgt_genome_v2_1.fa', 'r')
genomeLines = genomeFile.readlines()
genomeFile.close()
for start, stop, contig, name in extract:
    print(str(round(((current/total) * 100), 4)) + '% Done', end='\r')
    place = 1
    thisContig = False
    output = [genomeSlicer(start, stop, x, contig) for x in genomeLines]
    selection = [x for x in output if x is not None]
    bluSeq = ''.join(selection)
    upstreamGirl.append([name, bluSeq])
    current += 1

print (upstreamGirl)
outfile = open('bluUpstream.fa', 'w')
for name, seq in upstreamGirl:
    outfile.write('>' + name + '\n')
    outfile.write(seq + '\n')
outfile.close()
toc = time.clock()
print ('And it only took ' + str(toc - tic) + ' seconds')
