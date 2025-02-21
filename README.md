# PipelineProject
Python wrapper to analyze transcriptomes of Human herpesvirus 5, also known as Human cytomegalovirus or HCMV. (See Cheng et al., 2017: https://www.ncbi.nlm.nih.gov/pubmed/29158406). 

This pipeline uses transcriptomes from two donors at 2- and 6- days post-infection (dpi). You can find the data here manually: 

Donor 1 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896360

Donor 1 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896363

Donor 3 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896374

Donor 3 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896375

## Obtaining the sample datasets

### Method 1: Full data

For simplicity, just enter this code to download the following transcriptomes for testing: 

```
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030

wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033

wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044

wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045
```


Then uncompress with fasterq-dump

```
fasterq-dump SRR5660030

fasterq-dump SRR5660033

fasterq-dump SRR5660044

fasterq-dump SRR5660045
```

### Method 2: Partial data

Alternatively, instead of using wget and fasterq-dump for the full datasets, you can use the files provided in this repository, which contain only the first 10000 reads. 

## Running the script

Implement the Python script in the command line by calling `python PipelineProject.py`. It will generate a folder called `PipelineProject_Joshua_Melnick` if the folder doesn't already exist. It should try to move your data into that folder. The program has several checkpoints to avoid redoing the more complex commands by checking if the file it would have generated already exists. However, that means that if you want to run it with modified data, you should do it in a fresh folder. 

The output of this pipeline will be found in `PipelineProject.log`. I have provided the log produced by running this script with the full dataset. 

## Dependencies

You will need Biopython for it to run properly, as well as the R programming language with the sleuth and dplyr packages. You will also need Kallisto, SPAdes, Bowtie2, and BLAST+. 
