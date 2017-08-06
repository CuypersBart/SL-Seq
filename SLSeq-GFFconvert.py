import os, argparse

parser = argparse.ArgumentParser(description="This tools converts a standard GFF file to a format that is suited for SL-seq analysis. Briefly, the scripts expands the provided annotations to include the 5' UTRs of the genes, where a significant proportion of the SL-seq reads map.", epilog="For all questions, problems and suggestions please contact the author Bart Cuypers at bcuypers@itg.be")
parser.add_argument('GFF_file', metavar='GFF_file', type=str,
                    help='The input GFF file with gene annotations')
parser.add_argument('-o', dest='outputfile', type=str,metavar='',
                    help="Output filename (default: add 'SL_' prefix to inputfilename)")

args = parser.parse_args()

class Gene:

    def __init__(self, name, chromosome, start, stop, strand, col1,col2,col3,col6,col7,col8,col9):
        self.name = str(name)
        self.fullname=str(name)
        self.chromosome = str(chromosome)
        self.chromosomeSimple = str(chromosome.strip('_v01s1'))
        self.start = int(start)
        self.stop = int(stop)
        self.strand = str(strand)
        self.start_min = min(self.start, self.stop)
        self.stop_max = max(self.start, self.stop)
        self.stop_max2 = max(self.start, self.stop)
        self.fasta_gene = "NA"
        self.fasta_protein = "NA"
        self.col1=str(col1)
        self.col2=str(col2)
        self.col3=str(col3)
        self.col6=str(col6)
        self.col7=str(col7)
        self.col8=str(col8)
        self.col9=str(col9.strip('\n'))

        if self.strand == "+":
            self.IsPlus = 1
        elif self.strand == "-":
            self.IsPlus = 0

    def length(self):
        return self.stop - self.start

LiAbrev = ["tRNA", "snoR", "rRNA", "snRN"]
LiAbberations = ['%2C','%2F','%3A','%28','%29']
DictGeneAnnos = {}
LiChromosomes = set([])
with open(args.GFF_file, 'r') as table:
    for strLine in table:
        if strLine[0:2] != '##':
            LiLine = strLine.split('\t')
            if LiLine[2] == 'gene' or LiLine[2] == 'pseudogene':
                Anno = LiLine[8].split(';')
                StrGeneName = Anno[0].split('ID=')[1]
                if LiLine[0] not in LiChromosomes:
                    LiChromosomes.add(LiLine[0])
                if StrGeneName in DictGeneAnnos:
                    print LiLine
                DictGeneAnnos[StrGeneName]=Gene(StrGeneName,LiLine[0],LiLine[3],LiLine[4],LiLine[6],LiLine[0],LiLine[1],LiLine[2],LiLine[5],LiLine[6],LiLine[7],LiLine[8])

Genelists = {}
for chromosome in LiChromosomes:
    Genelists[chromosome] = []
for gene in DictGeneAnnos:
    gene = DictGeneAnnos[gene]
    Genelists[gene.chromosome].append(gene)

for chromosome in LiChromosomes:
    Genelists[chromosome]=sorted(Genelists[chromosome], key=lambda x: x.start_min)

removelist = []

for chromosome in LiChromosomes:
    length = len(Genelists[chromosome])
    if length == 1:
        gene=Genelists[chromosome][0]
        if gene.strand == '+':
            Genelists[chromosome][0].start_min = 1
        elif gene.strand == '-':
            Genelists[chromosome][0].stop_max = gene.stop_max+500
        print("warning chromosome "+ str(gene.chromosome) + " contains only 1 gene!")
    for i in range(len(Genelists[chromosome])):
        overlap = False
        gene = Genelists[chromosome][i]
        if i == 0:
            if gene.stop_max >= Genelists[chromosome][i+1].start_min:
                removelist.append(gene.name)
                overlap = True
        elif i == int(length-1):
            if gene.start_min <= Genelists[chromosome][i-1].stop_max:
                removelist.append(gene.name)
                overlap = True
        elif gene.stop_max >= Genelists[chromosome][i+1].start_min or gene.start_min <= Genelists[chromosome][i-1].stop_max:
                removelist.append(gene.name)
                overlap = True
        if overlap == False:
            if i == 0 and gene.strand == '+':
                Genelists[chromosome][i].start_min = 1
            elif i == 0 and gene.strand == '-':
                Genelists[chromosome][i].stop_max = Genelists[chromosome][i+1].start_min-1
            elif gene.strand == '+' and Genelists[chromosome][i-1].strand == "-":
                numb=(gene.start_min-Genelists[chromosome][i-1].stop_max2)/2
                Genelists[chromosome][i-1].stop_max = Genelists[chromosome][i-1].stop_max2 + numb
                Genelists[chromosome][i].start_min = Genelists[chromosome][i].start_min - numb
            elif gene.strand == '+':
                Genelists[chromosome][i].start_min = Genelists[chromosome][i-1].stop_max-1
            elif gene.strand == "-" and i == int(length-1):
                 Genelists[chromosome][i].stop_max = gene.stop_max+500
            elif gene.strand == "-":
                Genelists[chromosome][i].stop_max = Genelists[chromosome][i+1].start_min-1
            else:
                print(gene.strand)
                print('Warning gene ' + gene.name + ' lacks strand annotation')

if args.outputfile == None:
    outputfile = 'SL_' + args.GFF_file
else:
    outputfile = args.outputfile

with open(outputfile,'wb') as myfile:
    for chromosome in LiChromosomes:
        for gene in Genelists[chromosome]:
            if gene.name in removelist:
                print 'gene with overlapping region => gene discarded: ',
                print gene.name
            else:
                myfile.write(gene.col1 + '\t' + gene.col2 + '\t' + gene.col3 +'\t' + str(gene.start_min) + '\t' + str(gene.stop_max)  + '\t' + gene.col6 + '\t' + gene.col7 + '\t' + gene.col8 + '\t' + gene.col9 + '\n')
