# PipelineProject
Python wrapper to analyze transcriptomes of Human herpesvirus 5, also known as Human cytomegalovirus or HCMV. (See Cheng et al., 2017: https://www.ncbi.nlm.nih.gov/pubmed/29158406). 

This pipeline uses transcriptomes from two donors at 2- and 6- days post-infection (dpi). You can find the data here manually: 

Donor 1 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896360

Donor 1 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896363

Donor 3 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896374

Donor 3 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896375



For simplicity, just enter this code to download the following transcriptomes for testing: 

wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030

wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033

wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044

wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045


Then uncompress with fasterq-dump

fasterq-dump SRR5660030

fasterq-dump SRR5660033

fasterq-dump SRR5660044

fasterq-dump SRR5660045


After this step, you can implement the pipeline in Python. You will need Biopython for it to run properly, as well as the R programming language with the sleuth and dplyr packages. You will also need Kallisto, SPAdes, Bowtie2, and BLAST+. 
