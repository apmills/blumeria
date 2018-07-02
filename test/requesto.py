# Python
# Messing about with requests

import requests
import time

seqs = ['MANNNSDRQGLEPRVIRTLGSQALSGPSISNRTSSSE','ANPHFSKNVKEAMIKTASPTPLSTPIYRIAQACDRCRSKKTRCDGKRPQCSQCAAVGFECRISDKLLRKAYPKGYTESLEERVRELEAENKRLLALCDIKEQQISLVSQSRPQTSTDNTINGNFKHDLKDAPLNLSSTNIYLLNQTVNKQLQNGKMDGDNSGSAMSPLGAPPPPPHKDHLCDGVSCTNHLHVKPTSTSLNDPTAISFEQDEAPGLPAVKALKSMTTHQRSTQLATLVSLSIPRSTEEILFIPQLLTRIRQIFGFNSKQCLYTVSLLSSLKNRLPAPRLLAPSTSTKLKEKDEDK', 'KLDDDSAFVKRFQSTNLSEFVDLKKFLISLKFNINSFSKQSEKPANDQDDELLSLTEIKELLHLFFKFWSNQVPILNNDHFLIYFNNFVEVVKHLSTENLETNNTTKSTVTTNHEIFALKLLMMLQMGLLVKIKMEKIKYTVPKNPKAKYARLMAYYHQLSLIIPKNPYFLNMSTTSLPSLQLLSLASFYYLNVGDISAIYGVRGRIVSMAQQLRLHRCPSAVLSVHSNPVLQKFEQSERRLLFWAIYYVDVFASLQLGVPRLLKDFDIECALPISDVEYKDQLSMENEKADKKAKKIQLQGQVSSFSLQIIRFAKILGNILDSIFKRGMMDERITSEVALVHENALDNWRNQLPEMYYFQITVNGTVNLDEIRATNQRNTETKFDKKDIILFEKKILLLFYFLAKSMIHLPVIATKPLPKNVDNATKKKQSMFNNDSKGATNQDHMILDVDMTSPAIRTSSSYIILQQATNATLTIFQAINSMYLPLPLNVSRTLIRFSLLCARGSLEYTKGGALFLDNKNLLLDTIKDIENDRLLDLPGIASWHTLKLFDMSINLLLKAPNVKVERLDKFLEKKLNYYNRLMGLPPATTTSLKPLFGSQSKNSLENRQRTPNVKRENPEHEYLYGNDSNNNNNSEAGHSPMTNTTNGNK']

#r = requests.post('https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?', data=(('db','cdd'), ('queries',seq), ('tdata','hits')))

bort = open('../scripts/output/repeat.txt', 'r')
bortl = bort.readlines()

snort = eval('hits = ' + bortl[0].rstrip())

print(snort)
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
