import os
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import statistics
Entrez.email = "jmelnick@luc.edu"

if not os.path.isdir("PipelineProject_Joshua_Melnick"): #make directory if it doesn't exist
    os.system("mkdir PipelineProject_Joshua_Melnick")

os.chdir("PipelineProject_Joshua_Melnick") #this is the working directory

with open("PipelineProject.log", "w") as f: #create empty log file
    pass


#download NC_006273.2 reference genome
if not os.path.isdir("ncbi_dataset/data/GCF_000845245.1"): #no need to download repeatedly
    os.system("datasets download genome accession GCF_000845245.1 --include gff3,rna,cds,protein,genome,seq-report")
    os.system("unzip ncbi_dataset.zip") #unzip

cdsfile = 'ncbi_dataset/data/GCF_000845245.1/cds_from_genomic.fna' #list of CDS comes with the download, but the format is wrong; we want the header rows to just have the CDS name
CDS = 0 #number of CDS
with open('pipeproteins.fasta', 'w') as f: #create our trimmed CDS list
    with open(cdsfile, 'r') as g: #read the old list
        x = g.readlines() #iterate through lines of old list
        for line in x:
            if line[0] != ">": #nucleotide sequence, just copy unchanged
                f.write(line)
            elif line[0] == ">": #header line
                CDS += 1 #increase our CDS count
                splitspace = line.split(" ") #remove everything after the first space
                splitunder = splitspace[0][21:].split("_") #remove everything before the start of the name and then splil by underscore
                f.write(">" + splitunder[0] + splitunder[1] + '\n') #there is one underscore in our trimmed name, remove everything after the second underscorre

with open("PipelineProject.log", "a") as f: #write CDS count to log
    f.write("The HCMV genome (NC_006273.2) has {} CDS.\n".format(str(CDS)))


if not os.path.isfile("kalindex.idx"): #index the transcriptome with kallisto and call it kalindex
    os.system("kallisto index -i kalindex.idx pipeproteins.fasta")

trans = [] #list of transcriptomes
for file in os.listdir():
    if file[0:3] == "SRR" and "_" not in file:
        trans.append(str(file)) #list of transcriptome base names

if not os.path.isdir("results"): #make a results folder
    os.system("mkdir results")

for tran in trans: #quantify the TPM for each transcriptome using a 2 threads and a bootstrap value of 10
    if not os.path.isdir("results/{}".format(tran)):
        os.system("time kallisto quant -i kalindex.idx -o results/{0} -b 10 -t 2 {0}_1.fastq {0}_2.fastq".format(tran))

with open("sample_table.txt", "w") as f: #create sample table for sleuth of sample name, condition (2 or 6 dpi), and path
    f.write("sample condition path\nSRR5660030 2dpi results/SRR5660030\nSRR5660033 6dpi results/SRR5660033\nSRR5660044 2dpi results/SRR5660044\nSRR5660045 6dpi results/SRR5660045")


with open("sample_table.txt", "r") as f: #get lists of samples, conditions, and results folders
    samplelines = f.readlines()[1:]
    for line in range(len(samplelines)):
        samplelines[line] = samplelines[line].strip().split()

donors = {} #dictionary of samples to track donors
by_donor = {"1":[],"3":[]} #dictionary of samples belonging to each donor
for d in range(len(samplelines)):
    if d < 2: #samples 1 and 2 come from donor 1
        donors[samplelines[d][0]] = "1"
        by_donor["1"].append(samplelines[d][0]) #add sample to donor 1 dictionary
    else: #samples 3 and 4 come from donor 3
        donors[samplelines[d][0]] = "3"
        by_donor["3"].append(samplelines[d][0]) #add sample to donor 3 dictionary

with open("PipelineProject.log", "a") as f: #headers
    f.write("\nsample\tcondition\tmin_tpm\tmed_tpm\tmean_tpm\tmax_tpm\n")

for sample in samplelines:
    with open("{}/abundance.tsv".format(sample[2]), "r") as f:
        cols = f.readlines()[1:] #ignore header
        for col in range(len(cols)):
            cols[col] = float(cols[col].split()[4]) #get TPM column
    with open("PipelineProject.log", "a") as f: #get sample, condition, min TPM, median TPM, mean TPM, and max TPM
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(sample[0],sample[1],min(cols),statistics.median(cols),statistics.mean(cols),max(cols)))


rsleuth = '''
#load package
library(sleuth)

#read table describing samples and kallisto output, assign to variable called stab
stab = read.table("sample_table.txt",header=TRUE)

#initialize sleuth object with sleuth_prep function
so = sleuth_prep(stab)

#fit model comparing the two conditions
so = sleuth_fit(so,~condition,'full')

#fit reduced model to compare in likelihood ratio test
so = sleuth_fit(so,~1,'reduced')

#perform likelihood ratio test for differential expression between conditions
so = sleuth_lrt(so,'reduced','full')

#load dplyr package for data frame filtering
library(dplyr)

#extract test results from dplyr object
sleuth_table = sleuth_results(so,'reduced:full','lrt',show_all=FALSE)

#filter most significant results (FDR/qval < 0.05) and sort by pval
sleuth_significant = dplyr::filter(sleuth_table, qval <= 0.05) |> dplyr::arrange(pval)

#write FDR < 0.05 transcripts to file
write.table(dplyr::select(sleuth_significant, target_id, test_stat, pval, qval), file="fdr05_results.txt", quote=FALSE, sep="\t",row.names=FALSE)
'''
with open("rsleuth.R","w") as f: #write the R script
    f.write(rsleuth)

if not os.path.isfile("fdr05_results.txt"): #run the R script if you haven't done it before
    os.system("Rscript rsleuth.R")

with open("PipelineProject.log", "a") as f: #write output to log
    f.write('\n') #separate the section with a blank line
    with open("fdr05_results.txt", "r") as l: #read the output of the R script
        for line in l.readlines(): #loop through the lines of the output from R
            f.write(line) #add them to the log
    f.write('\n') # add another blank line for the next section


#build bowtie2 index called HCMV
if not os.path.isfile("HCMV.1.bt2"):
    os.system("bowtie2-build ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna HCMV")

for samp in samplelines: #map transcriptomes to index
    if not os.path.isfile("{}_mapped_1.fq.gz".format(samp[0])):
        os.system("bowtie2 --quiet -x HCMV -1 {0}_1.fastq -2 {0}_2.fastq -S HCMV_{0}.sam --al-conc-gz {0}_mapped_%.fq.gz".format(samp[0]))
    with open("PipelineProject.log", "a") as f: #output number of reads mapped and unmapped to log
        readtotal = 0
        filtered = 0
        with open("HCMV_{}.sam".format(samp[0]), "r") as s: #sam file
            #readtotal = int(os.system("wc -l {}".format(s))) #count the number of lines for original read number, exclude first 3
            #filtered = int(os.system('grep "YF:" {} | wc -l'.format(s))) #count the number of lines with the YF flag that indicates it was filtered out
            for line in s.readlines()[3:]: #count the number of lines for original read number, exclude first 3
                readtotal += 1
                if "YF:" in line: #count the number of lines with the YF flag that indicates it was filtered out
                    filtered += 1
        remaining = readtotal - filtered #the number of reads after filtration by Bowtie2
        f.write("Donor {0} ({1}) had {2} read pairs before Bowtie2 filtering and {3} read pairs after.\n".format(donors[samp[0]],samp[1],readtotal,remaining))

with open("PipelineProject.log", "a") as f:
    f.write('\n') #separate sections in the log

for patient in by_donor: #assemble in spades
    com = "spades.py -k 77 --only-assembler --pe-1 1 {0}_1.fastq --pe-2 1 {0}_2.fastq --pe-1 2 {1}_1.fastq --pe-2 2 {1}_2.fastq -o donor_{2}_assembly/".format(by_donor[patient][0],by_donor[patient][1],patient)
    with open("PipelineProject.log", "a") as f:
        f.write(com + '\n')
    if not os.path.isdir("donor_{}_assembly".format(patient)):
        os.system(com)

#obtain sequences for Betaherpesvirinae subfamily
if not os.path.isfile("ncbi_dataset/data/genomic.fna"):
    os.system("datasets download virus genome taxon betaherpesvirinae --refseq --include genome")
    os.system("unzip ncbi_dataset.zip")

#build Betaherpesvirinae database
if not os.path.isfile("betaherpesvirinae.nhr"):
    os.system("makeblastdb -in ncbi_dataset/data/genomic.fna -out betaherpesvirinae -dbtype nucl")

for don in by_donor: #loop through our donors
    contigs = list(SeqIO.parse("donor_{}_assembly/contigs.fasta".format(don), "fasta")) #list of contigs
    longid = "" #name of the longest contig
    longseq = "" #sequence of the longest contig
    i = 0 #index
    longin = 0 #index tracker
    for c in contigs: #loop through that donor's contigs
        if len(c.seq) > len(longseq): #if this is the longest contig so far, update our variables
            longid = c.id
            longseq = c.seq
            longin = i
        i += 1
    SeqIO.write(contigs[longin], "d{}_contig.fasta".format(don), "fasta") #write longest contig to a file
    #run blast+
    blastcommand = 'blastn -query d{0}_contig.fasta -db betaherpesvirinae -out blast{0}.txt -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle"'.format(don)
    os.system(blastcommand)
    with open("PipelineProject.log", "a") as f: #add the output to the log
        f.write("\nDonor{}:\n".format(don))
        f.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n") #header
        with open("blast{}.txt".format(don), "r") as j:
            for line in j.readlines()[0:10]: #top 10 results
                f.write(line)

