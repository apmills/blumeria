# Python
# Messing about with requests

import requests
import time
from Bio import Phylo
from Bio.Phylo.Applications import PhymlCommandline
from Bio import AlignIO
from subprocess import call


#call('clustalo -i ../scripts/output/unaligned.fa -o phyAlign.clu --force --outfmt=clu', shell=True)
#AlignIO.convert('phyAlign.clu', 'clustal', 'phyAlign.phy', 'phylip-relaxed')
#cmdline = PhymlCommandline(input='phyAlign.phy', alpha='e', bootstrap=1, sequential=False)
#call(str(cmdline), shell=True)
my_tree = Phylo.read("phyAlign.phy_phyml_tree.txt", "newick")
Phylo.draw(my_tree, show_confidence=False)

quit()

seqs = ['MANNNSDRQGLEPRVIRTLGSQALSGPSISNRTSSSE','ANPHFSKNVKEAMIKTASPTPLSTPIYRIAQACDRCRSKKTRCDGKRPQCSQCAAVGFECRISDKLLRKAYPKGYTESLEERVRELEAENKRLLALCDIKEQQISLVSQSRPQTSTDNTINGNFKHDLKDAPLNLSSTNIYLLNQTVNKQLQNGKMDGDNSGSAMSPLGAPPPPPHKDHLCDGVSCTNHLHVKPTSTSLNDPTAISFEQDEAPGLPAVKALKSMTTHQRSTQLATLVSLSIPRSTEEILFIPQLLTRIRQIFGFNSKQCLYTVSLLSSLKNRLPAPRLLAPSTSTKLKEKDEDK', 'KLDDDSAFVKRFQSTNLSEFVDLKKFLISLKFNINSFSKQSEKPANDQDDELLSLTEIKELLHLFFKFWSNQVPILNNDHFLIYFNNFVEVVKHLSTENLETNNTTKSTVTTNHEIFALKLLMMLQMGLLVKIKMEKIKYTVPKNPKAKYARLMAYYHQLSLIIPKNPYFLNMSTTSLPSLQLLSLASFYYLNVGDISAIYGVRGRIVSMAQQLRLHRCPSAVLSVHSNPVLQKFEQSERRLLFWAIYYVDVFASLQLGVPRLLKDFDIECALPISDVEYKDQLSMENEKADKKAKKIQLQGQVSSFSLQIIRFAKILGNILDSIFKRGMMDERITSEVALVHENALDNWRNQLPEMYYFQITVNGTVNLDEIRATNQRNTETKFDKKDIILFEKKILLLFYFLAKSMIHLPVIATKPLPKNVDNATKKKQSMFNNDSKGATNQDHMILDVDMTSPAIRTSSSYIILQQATNATLTIFQAINSMYLPLPLNVSRTLIRFSLLCARGSLEYTKGGALFLDNKNLLLDTIKDIENDRLLDLPGIASWHTLKLFDMSINLLLKAPNVKVERLDKFLEKKLNYYNRLMGLPPATTTSLKPLFGSQSKNSLENRQRTPNVKRENPEHEYLYGNDSNNNNNSEAGHSPMTNTTNGNK']

#r = requests.post('https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?', data=(('db','cdd'), ('queries',seq), ('tdata','hits')))

place = 1
thisContig = False
def genomeSlicer(start, end, line, contig):
    global place
    global thisContig
    if thisContig and place >= end:
        thisContig = False
        return
    if line[0] == '>':
        if contig in line:
            thisContig = True
        return
    if thisContig:
        if place > start - 50:
            output = []
        else:
            place += 50
            return
        for nuc in line.rstrip():
            if place >= start and place < end:
                output.append(nuc)
            place += 1
        return ''.join(output)
    else:
        return


genome = open('../db/blumeria/latest/Bgt_genome_v2_1.fa')
print('Starting...')
output = [genomeSlicer(360690, 361973, x, 'Bgt_ctg-10000_consensus') for x in genome]
selection = [x for x in output if x is not None]
print (''.join(selection))
quit()

print (r.text)


searchID = r.text.split('#')[2].split()[1].rstrip()

print (searchID)

tick = True
done = False

while not done:
    result = requests.post('https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?', data={'cdsid':searchID, 'cddefl':'true'})
    status = result.text.split('#')[4].split()[1].rstrip()
    if status == '3':
        done = False
    else:
        done = True
    time.sleep(3)
    if tick:
        print ('Tick')
        tick = False
    else:
        print ('Tock')
        tick = True

lines = result.text.split('\n')

tags = []
for line in lines:
    # if line[0] == '#' or line[0] == '' or line[:5] == 'Query':
    #     continue
    print (line)
    if line[:2] == 'Q#':
        tags.append(line.split('\t')[8])


print (tags)
