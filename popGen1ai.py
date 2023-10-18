import gffutils
from Bio import SeqIO
import vcf
import sys

# Check the number of command line arguments
if len(sys.argv) != 3:
    print("Usage: python script.py <species> <gene_list_file>")
    sys.exit(1)

species = sys.argv[1]
gene_list_file = sys.argv[2]

# Define a dictionary to map species to files
species_files = {
    "Tvivax": {
        "gff": "TriTrypDB-28_TvivaxY486.gff",
        "cds_fasta": "TriTrypDB-28_TvivaxY486_AnnotatedCDSs.fasta",
        "snps_vcf": "raw_variants_AllTvivax-dp5.recode-biallelic.vcf.gz",
    },
    "Tcongolense": {
        "gff": "IL3000a.gff",
        "cds_fasta": "TriTrypDB-32_TcongolenseIL3000_AnnotatedCDSs.fasta",
        "snps_vcf": "Final_vcf_Tcongo.vcf.gz",
    },
}

if species not in species_files:
    print("Accepted species inputs are Tvivax or Tcongolense")
    sys.exit(1)

# Load input files based on the selected species
species_info = species_files[species]
db = gffutils.create_db(species_info["gff"], dbfn=f"{species}_DB.db", force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
indexed_cds = SeqIO.index(species_info["cds_fasta"], "fasta")
snps = vcf.Reader(filename=species_info["snps_vcf"])

with open(gene_list_file) as gene_list:
    genes = gene_list.read().splitlines()

# Open the output file for writing
with open("snpsSummary", "w") as out_genes:
    for curr_gene in genes:
        c_strand = db[curr_gene].strand
        c_start, c_end, c_chrom = db[curr_gene].start, db[curr_gene].end, db[curr_gene].chrom

        if c_strand == "+":
            c_seq = indexed_cds[curr_gene].seq
        else:
            c_seq = indexed_cds[curr_gene].seq.reverse_complement()

        c_chrom_bin = c_chrom if species == "Tcongolense" else {"TvY486_01": 1, "TvY486_02": 2, "TvY486_03": 3, "TvY486_04": 4, "TvY486_05": 5, "TvY486_06": 6, "TvY486_07": 7, "TvY486_08": 8, "TvY486_09": 9, "TvY486_10": 10, "TvY486_11": 11, "TvY486_bin": 12}.get(c_chrom, 0)

        if c_chrom_bin == 0:
            print(f"WARNING: Unknown chromosome ID for gene {curr_gene}")
            continue

        curr_exons = list(db.children(curr_gene, featuretype="exon"))
        if not curr_exons:
            print(f"No exons found for gene {curr_gene}! Skipping.")
            continue

        snp_region = list(snps.fetch(c_chrom_bin, c_start - 1, c_end))

        alt_seq = c_seq.tomutable()
        counts = {'Ref': 0, 'Het': 0, 'Alt': 0}

        for c_snp in snp_region:
            sample_str = str(c_snp.samples)
            counts['Ref'] += sample_str.count("GT=0/0")
            counts['Het'] += sample_str.count("GT=0/1")
            counts['Alt'] += sample_str.count("GT=1/1")

            if c_snp.REF != alt_seq[c_snp.POS - c_start]:
                print(f"PROBLEM with SNP at position {c_snp.POS} for gene {curr_gene}!")

        alt_seq = alt_seq.toseq()

        if c_strand == "-":
            c_seq = c_seq.reverse_complement()
            alt_seq = alt_seq.reverse_complement()
            snp_region.reverse()
            counts['Ref'], counts['Het'], counts['Alt'] = counts['Ref'], counts['Het'], counts['Alt']

        n_snp = 0

        with open("snpsSummary", "a") as out_genes:
            for i in range(len(c_seq)):
                if c_seq[i] != alt_seq[i]:
                    j = i // 3
                    x = i + 1
                    z = j + 1
                    dnds = "Syn" if c_seq.translate()[j] == alt_seq.translate()[j] else "Non"
                    snp = snp_region[n_snp]

                    print(curr_gene, snp.CHROM, c_strand, snp.POS, snp.REF, snp.ALT[0], x, c_seq[i], alt_seq[i], z, c_seq.translate()[j], alt_seq.translate()[j], dnds, counts['Ref'], counts['Het'], counts['Alt'], file=out_genes)
                    n_snp += 1

            if n_snp > 0:
                print("", file=out_genes)
                print(f">{curr_gene}-mRNA-Ref", file=out_genes)
                print(c_seq, file=out_genes)
                print(f">{curr_gene}-mRNA-Alt", file=out_genes)
                print(alt_seq, file=out_genes)
                print(f">{curr_gene}-Ref", file=out
