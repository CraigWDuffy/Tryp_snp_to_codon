import gffutils
from Bio import SeqIO, Seq
import vcf
import sys

#from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import IUPAC
#import Bio.codonalign
#from Bio.codonalign.codonseq import cal_dn_ds

print(sys.argv[1])
print(sys.argv[2])

#For current testing sys.argv[2] = 'seqs' and sys.argv[1] = 'Tcongolense'

# Step 1 - Load in input files
# Input files = GFF, CDSfasta, SNPs, gene list

#Read in the genes to be processed


TvDB = gffutils.create_db("TriTrypDB-28_TvivaxY486.gff", dbfn="TvDB.db", force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)

TcDB = gffutils.create_db("IL3000a.gff", dbfn="TcDB.db", force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)

with open(sys.argv[2]) as f:
#with open("seqs") as f:
    genes = f.read().splitlines()

print(genes)

#Determine speces and read in the pop data and gene seqs
if sys.argv[1] == "Tvivax":
    #if "Tcongolense" == "Tcongolense":
    inDB = gffutils.FeatureDB('TvDB.db', keep_order=True)
    indexedCDS = SeqIO.index("TriTrypDB-28_TvivaxY486_AnnotatedCDSs.fasta", "fasta")
    snps = vcf.Reader(filename="raw_variants_AllTvivax-dp5.recode-biallelic.vcf.gz")
elif sys.argv[1] == "Tcongolense":
    #elif "Tcongolense" == "Tvivax":
    inDB = gffutils.FeatureDB('TcDB.db', keep_order=True)
    indexedCDS = SeqIO.index("TriTrypDB-32_TcongolenseIL3000_AnnotatedCDSs.fasta", "fasta")
    snps = vcf.Reader(filename="Final_vcf_Tcongo.vcf.gz")
else:
    print("Accepted species inputs are Tbrucei, Tcongolense or Tvivax")
    sys.exit()

# Have read in the data - next step is to get the positional information for each gene and extract the snps.

#loop over all genes


outGenes = open("snpsSummary", "w")

for currGene in genes:
    #print(currGene)
    print(currGene, file=outGenes)
    print("Gene Chromosome Strand Pos(chromosomal) RefBase(Genomic) AltBase(Genomic) Pos(genic) RefBase(genic) AltBase(genic) Codon RefAminoAcid AltAminoAcid Syn/Non RefCount HetCount AltCount", file=outGenes)
    cStrand = inDB[currGene].strand
    cStart = inDB[currGene].start
    cEnd = inDB[currGene].end
    cChrom = inDB[currGene].chrom
    #print(cChrom)
    if cStrand == "+":
        cSeq = indexedCDS[currGene].seq
    else:
        cSeq = indexedCDS[currGene].seq.reverse_complement()
    
    SNPRegion = []
    
    if sys.argv[1] == "Tvivax":
        cChromBin = 1 if cChrom == "TvY486_01" else 2 if cChrom == "TvY486_02" else 3 if cChrom == "TvY486_03" else 4 if cChrom == "TvY486_04" else 5 if cChrom == "TvY486_05" else 6 if cChrom == "TvY486_06" else 7 if cChrom == "TvY486_07" else 8 if cChrom == "TvY486_08" else 9 if cChrom == "TvY486_09" else 10 if cChrom == "TvY486_10" else 11 if cChrom == "TvY486_11" else 12 if cChrom == "TvY486_bin" else 0
    
        if cChromBin == 0:
            print("WARNING: Unknown chromosome ID")
    else:
        cChromBin = cChrom
        
    currExons = []
    for i in inDB.children(currGene, featuretype="exon"):
    	currExons.append(i)


    if len(currExons) == 0:
        print ("No exons found! Exiting")
        sys.exit()
        
    try:
        for record in snps.fetch(cChromBin, cStart-1, cEnd):
            SNPRegion.append(record)
    except:
        print(currGene + " - No snps!")
    
    AltSeq = cSeq.tomutable()
    
    RefCount=[]
    HetCount=[]
    AltCount=[]
    
    for cSNP in SNPRegion:
        RefCount.append(str(cSNP.samples).count("GT=0/0"))
        HetCount.append(str(cSNP.samples).count("GT=0/1"))
        AltCount.append(str(cSNP.samples).count("GT=1/1"))
        #print(RefCount, HetCount, AltCount)
        myAdjust=0
        #print(cSNP.REF)
        #print(cSNP.ALT)
        #print(str(cSNP.POS) + " " + str(cStart) + " " + str(cEnd))
        if len(currExons) > 2:
            print("Uh oh, more than 2 exons found and I'm not programmed to handle that yet")
            sys.exit()
        elif len(currExons) == 2:
            if cSNP.POS > currExons[0][4] & cSNP.POS < currExons[1][3]:
                pass
            if cSNP.POS >= currExons[1][3]:
                myAdjust = currExons[1][3] - currExons[0][4] - 1
        
        if cSNP.REF != AltSeq[cSNP.POS - cStart - myAdjust]:
            print(str(cSNP.REF) + " " + str(AltSeq[cSNP.POS - cStart]) + " " + str(cSNP.POS))
            print("PROBLEM!")
            sys.exit()
        #print(AltSeq[cSNP.POS - cStart])
        AltSeq[cSNP.POS - cStart - myAdjust] = str(cSNP.ALT[0])
    AltSeq = AltSeq.toseq()

    if cStrand == "-":
        cSeq = cSeq.reverse_complement()
        AltSeq = AltSeq.reverse_complement()
        SNPRegion.reverse()
        RefCount.reverse()
        HetCount.reverse()
        AltCount.reverse()
        # Can then access the data in savedRegion - need to parse it through to the seqIO files somehow
    i=0
    nSNP=0

    if len(SNPRegion) > 0:
        while i < len(cSeq):
            if cSeq[i] != AltSeq[i]:
                j = i//3
                x = i+1
                z = j+1
                dnds = "Syn"
                if cSeq.translate()[j] != AltSeq.translate()[j]:
                    dnds = "Non"
                print(currGene, SNPRegion[nSNP].CHROM, cStrand, SNPRegion[nSNP].POS, SNPRegion[nSNP].REF, SNPRegion[nSNP].ALT[0], x, cSeq[i], AltSeq[i], z, cSeq.translate()[j], AltSeq.translate()[j], dnds, RefCount[nSNP], HetCount[nSNP], AltCount[nSNP], file=outGenes)
                nSNP=nSNP+1
        
            i = i+1

        print("", file=outGenes)
        print(">" + currGene + "-mRNA-Ref", file=outGenes)
        print(cSeq, file=outGenes)
        print(">" + currGene + "-mRNA-Alt", file=outGenes)
        print(AltSeq, file=outGenes)
        print(">" + currGene + "-Ref", file=outGenes)
        print(cSeq.translate(), file=outGenes)
        print(">" + currGene + "-Alt", file=outGenes)
        print(AltSeq.translate(), file=outGenes)
        print("", file=outGenes)
    else:
        print("0 SNPs\n", file=outGenes)
        
    
    

outGenes.close()

