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

infile = open('../db/m_o_genes.fasta')

seq = fastaFinder(infile.readlines(), 'MGG_04108')

print (seq)
