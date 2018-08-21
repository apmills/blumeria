# Python
# $ pathway.py <Pathway Name> <Verbose>
# Pathway Name arg should be a folder containing fasta files of AA sequences for relevant proteins or 'KEGG' for a KEGG pathway search
# Finds the best matches for a collection of S. cerevisiae genes in the various databases

from subprocess import call
import sys, glob, requests, time, re, os, textwrap, pandas
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import Phylo, AlignIO, SeqIO
from Bio.Phylo.Applications import PhymlCommandline

pathway = sys.argv[1]
if pathway.upper() == 'KEGG':
    fetch = True
elif pathway.upper().strip('-') == 'HELP':
    print (textwrap.dedent("""\
        ---
        Help:

        This is pathway.py.
        A tool for rapidly searching for transcription factor binding site motifs
        within promoter sequences of a set of genes. By default, it is searching within
        a Blumeria graminis f. sp. tritici genome and a set of related fungi.

        To begin, first identify the KEGG pathway or module that you are studying.
        This will have an identifier such as ko12345 or M12345. Then simply run this
        program as '$ python pathway.py KEGG' and enter this identifier when prompted.

        It is also possible to supply your own set of enzymes to query. Simply create
        a directory named e.g. testPathway and run the program as
        '$ python pathway.py testPathway' instead.

        The program will next ask for a binding site motif in IUPAC format. This should
        correspond to the nucleotide motif that your TF of interest binds to.

        Finally, the program will run automatically from this point and produce some
        output. The key output is in output/<Pathway>/results.tsv.
        ---
    """))
    quit()
else:
    fetch = False

if 'verbose' in sys.argv:
    verbose = True
else:
    verbose = False
###########################################
# Collect enzyme data for a given pathway #
###########################################
def getKegg(pathCode):
    # Returns a dictionary of {Enzyme Name : EC number} for a whole KEGG pathway
    try:
        result = requests.get('http://rest.kegg.jp/get/' + pathCode)
    except requests.exceptions.ConnectionError:
        print ('Could not connect to KEGG!')
        print ('Or you used an invalid KEGG pathway ID!')
        return
    lines = result.text.split('\n')
    enzymes = {} # {Enzyme Name: EC #}
    collect = False
    for line in lines:
        if 'NAME' in line:
            pathName = re.search(r'NAME\s+(.+)', line)[1]
            continue
        if 'ORTHOLOGY' in line:
            collect = True
        elif 'COMPOUND' in line:
            collect = False
        elif 'CLASS' in line:
            collect = False
        if collect:
            keggMatch = re.search(r'(K\d+)\s+', line)
            keggNo = keggMatch[1]
            # if not inYeast(keggNo):
            #     print (keggNo + " ain't in yeast")
            #     continue
            nameMatch = re.search(r'K\d+\s\s(.+)\s\[E', line)
            if not nameMatch:
                continue
            ecMatch = re.search(r'\[EC:((\d+\.\d+\.\d+\.\d+)[\s\]])+', line)
            if not ecMatch:
                continue
            name = nameMatch[1]
            multi = False
            if len(name.split(' / ')) > 1:
                name = name.split(' / ')
                multi = True
            else:
                name = [name]
            ec = ecMatch[0].strip('EC:[]')
            if multi:
                ec = ec.split()
            elif len(ec.split()) > 1:
                ec = [ec.split()[0]]
            else:
                ec = [ec]
            for n, e in zip(name, ec):
                if e.strip('0123456789.') == '':
                    # Only include valid EC numbers (no ambiguous ones like 1.3.99.-)
                    enzymes[e] = n
    return enzymes, pathName

def getUniprot(ec):
    # Returns a fungal protein sequence (with FASTA header) for a given EC number
    try:
        url = 'https://www.uniprot.org/uniprot/?query=ec%3A' + ec + '&fil=organism%3A559292&sort=score&format=fasta'
    except requests.exceptions.ConnectionError:
        print('Could not connect to UniProt!')
        return
    result = requests.get(url)
    if len(result.text) < 10:
        # If there are no results for yeast, look more broadly at all Leotiomyceta for info
        url = 'https://www.uniprot.org/uniprot/?query=ec%3A' + ec + '+AND+taxonomy%3A"leotiomyceta+%5B716546%5D"&sort=score&format=fasta'
        result = requests.get(url)
        if len(result.text) < 10:
            return None
    lines = result.text.split('\n')
    active = False
    seqs = []
    for line in lines:
        if '>' in line:
            active = True
            header = line + '\n'
            continue
        if active:
            if '>' in line:
                break
            seqs.append(line)

    seq = ''.join(seqs)
    if len(seq) > 3500:
        # Don't include ludicrously long polypeptides
        return None
    return header + seq

if fetch:
    while pathway[:2].lower() != 'ko' and pathway[0].upper() != 'M':
        pathway = input('Enter a KEGG pathway (ko######) or module (M#####):\n')
        if (pathway[:2].lower() != 'ko' and pathway[0].upper() != 'M') or len(pathway) < 6 or len(pathway) > 7:
            print ('You must enter a vaild KEGG pathway or module ID')
            print ('---')
    if pathway[0] == 'm':
        pathway = pathway.upper()
    elif pathway[0] == 'K':
        pathway = pathway.lower()
    print('Accessing KEGG...')
    enzymeDict, pathName = getKegg(pathway)
    nickName = ''
    sname = pathName.split()
    if len(sname) == 1:
        nickName = sname[0][:6]
    elif len(sname) == 2:
        nickName = sname[0][:4].capitalize()
        nickName += sname[1][:2].capitalize()
    else:
        for w in sname:
            nickName += w[0].upper()
    if verbose:
        print (enzymeDict)
    print ('Found pathway "' + pathName + '"')
    print ('Abbreviated to ' + nickName)
    print ('Collected the following enzymes from the pathway')
    print ('EC Number\tEnzyme Name')
    for key in enzymeDict.keys():
        print ('{} \t{}'.format(key, enzymeDict[key]))
    # garbage = []
    print('\nEnzymes with no suitable hits on UniProt will be removed')
    print('Accessing UniProt...')
    totalPathway = [] # [[Enzyme Name, EC #, FASTA entry, UniProt acc], ...]
    seen = []
    accCount = {}
    for key in enzymeDict.keys():
        test = getUniprot(key)
        if test:
            if verbose:
                print (enzymeDict[key] + ':')
                print (test)
            acc = re.search(r'\|(.+)\|', test)[1]
            if acc not in seen:
                totalPathway.append([enzymeDict[key], key, test, acc])
                seen.append(acc)
                accCount[acc] = 1
            else:
                accCount[acc] = accCount[acc] + 1

    if not os.path.exists('./' + nickName):
        os.makedirs('./' + nickName)
    call('rm -f .' + nickName + '/*.fasta', shell=True)
    indexFile = open('./' + nickName + '/index.csv', 'w')
    indexFile.close()
    for name, ec, fasta, acc in totalPathway:
        outfile = open('./' + nickName + '/' + acc + '.fasta', 'w')  # Here
        outfile.write(fasta)
        outfile.close()
        indexFile = open('./' + nickName + '/index.csv', 'a')
        if accCount[acc] == 1:
            indexFile.write(acc + ',' + name + '\n')
        else:
            indexFile.write(acc + ',' + name + ' (+' + str(accCount[acc] - 1) + ')\n')
    indexFile.close()

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

if pathway.upper() == 'KEGG':
    pathway = nickName

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
        if not os.path.exists('output/' + pathway):
            os.makedirs('output/' + pathway)
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

# Setting reverse complement parameter for FIMO later on
rc = ''
while rc.lower() != 'y' and rc.lower() != 'n':
    rc = input('Search for motif on reverse strand? (y/n): ')

if rc.lower() == 'n':
    norc = '-norc '
else:
    norc = ' '

def blasterMaster(query, blu=False):
    # This function handles the BLAST search (through the CLI) and the results
    searched_genomes = []
    global verbose
    if blu:
        if not os.path.exists('output/' + pathway + '/blu'):
            os.makedirs('output/' + pathway + '/blu')
        call('tblastn -db '+ blu + ' -query ' + query + ' -out output/' + pathway + '/blu/' + query.split('/')[-1][:-3] + '.blast' + ' -num_threads 4 -max_target_seqs 1 -outfmt "7 sseqid evalue"', shell=True)
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
            if verbose:
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
    global verbose
    i = 0
    seqs = []
    collect = False
    while i < len(lines):
        if lines[i][0] == '>':
            if seqid in lines[i]:
                collect = True
                if verbose:
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
print ('Finding matching sequences...')
enzymes = [] # [[Enzyme ID X, Fungus Y genome, Enzyme X gene ID in Y, Enzyme name],...]
indexfile = open(pathway + '/index.csv', 'r')
accIndex = {}
enzymeIndex = {}
for line in indexfile.readlines():
    sline = line.split(',')
    accIndex[sline[0]] = sline[1].rstrip()
for enzymeId in glob.iglob(pathway + '/*.fasta'):
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
print ('Finding matches in Blumeria...')
# This should be a list of [[Genes FASTA file, Genome name, Promoter FASTA file (in same directory as genes)]]
queryGenomes = [['../db/blumeria/latest/bluGenes.fa', 'bgt', 'bluUpstream0.fa'], ['../db/blumeria/latest/bghGenes.fa', 'bgh', 'bghUpstream.fa']]
for blufspGenome, blufsp, fspUp in queryGenomes:
    bluEnzymes = []
    for enzymeId in glob.iglob(pathway + '/*.fasta'):
        acc = enzymeId.split('/')[-1].split('.')[0]
        hit = blasterMaster(enzymeId, blu=blufspGenome)
        enzymeIndex[hit] = accIndex[acc]
        bluEnzymes.append([enzymeId.split('/')[-1].split('.')[0], hit, accIndex[acc]])
    if verbose:
        print (bluEnzymes)

    seqFile = open('output/' + pathway + '/' + blufsp + 'Upstream.fa', 'w')
    for x, ide, fullName in bluEnzymes:
        # upper = genome[:-11] + 'upstream.fasta'
        gnome = open('../db/blumeria/latest/' + fspUp, 'r')
        lines = gnome.readlines()
        i = 0
        seqs = []
        collect = False
        while i < len(lines):
            if lines[i][0] == '>':
                if lines[i].rstrip() == '>' + ide:
                    collect = True
                    if verbose:
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

if verbose:
    print (enzymeIndex)
################################################
# Identify motifs in those upstreams with MEME #
################################################

print ('Running FIMO...')
matrices = [matrix for matrix in glob.glob('output/' + pathway + '/seqMatrix*')]
locs = ['output/' + pathway + '/fimo/blumeria','output/' + pathway + '/fimo/control','output/' + pathway + '/fimo/training']
organisms = {'Af':'A. fumigatus','AN':'A. nidulans','SS':'S. sclerotiorum','NC':'N. crassa','MG':'M. oryzae','BC':'B. cinerea','Bg':'f.sp. tritici','BL':'f.sp. hordei'}
results = [] # [[Seq ID, Enzyme name, strand, p-value, matched seq],[Another motif result], ...]
call('cat output/'+pathway+'/*.fa > output/'+pathway+'/allUpstream.fasta', shell=True)
if not os.path.exists('output/' + pathway + '/fimo'):
    os.makedirs('output/' + pathway + '/fimo')
for matrix in matrices:
    call('cat ' + matrix + ' | ~/meme/libexec/meme-5.0.1/matrix2meme | ~/meme/bin/fimo -oc output/' + pathway + '/fimo ' + norc + '-thresh 0.001 - output/' + pathway + '/allUpstream.fasta', shell=True)
    fimoFile = open('output/' + pathway + '/fimo/fimo.tsv', 'r')
    for line in fimoFile.readlines():
        sline = line.split('\t')
        if sline[0] == 'motif_id':
            continue # Skip first line
        if len(line) < 4:
            break # Stop before the end
        ide = sline[2]
        strand = sline[5]
        pval = sline[7]
        qval = sline[8]
        match = sline[-1].rstrip()
        results.append([ide, enzymeIndex[ide], strand, pval, qval, match])
    fimoFile.close()

resultFile = open('output/' + pathway + '/results.tsv', 'w')
resultFile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('Organism','Seq ID','Enzyme','Strand','p-value','q-value','Match Seq'))
for ide, name, strand, pval, qval, match in results:
    resultFile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(organisms[ide[:2]],ide,name,strand,pval,qval,match))
resultFile.close()
df = pandas.read_table('output/' + pathway + '/results.tsv')
# print (df)
summary = df.groupby(['Organism', 'Enzyme']).size().reset_index(name='Motif Count')
summary.to_csv('output/' + pathway + '/summary.tsv', sep='\t')
print ('Motif Occurences per Enzyme Promoter Sequence per Organism')
print (summary)
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
