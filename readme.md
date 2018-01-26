# AgRenSeq

## Description
AgRenSeq is a pipeline to identify candidate resistance (_R_) genes in plants directly from a diversity panel. The diversity panel needs to be sequenced (_R_ gene enrichment sequencing - RenSeq) and phenotyped. Phenotype scores need to be converted to AgRenSeq scores that assign positive values to resistance and negative values to suscetibility. An intermediate phenotype should have an AgRenSeq score close to zero.

For RenSeq you will need a bait library that targets _R_ genes in your plant species. A bait library for _Aegilops tauschii_ can be found _here_ (in preparation). We reccomend [Arbor biosciences](http://www.arborbiosci.com/) for synthesis of baits. They also offer the [enrichment service](http://www.arborbiosci.com/products/myreads-ngs-services-for-targeted-sequencing/). 

More about this method can be found in the manuscript [http://biorxiv.org/cgi/content/short/248146v1](http://biorxiv.org/cgi/content/short/248146v1)




## Pre-requisites
### JRE 1.6
Make sure you have the Java Runtime Environments 1.6 or higher. Download from [http://java.com](http://java.com)

### Sequence Quality Trimmer
We recommend to quality trim your sequences before the k-mer counting. Check your quality with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). Use a tool such as Trimmomatic ([http://www.usadellab.org/cms/?page=trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)) for read preprocessing.

### Jellyfish
Download and install jellyfish from [http://www.genome.umd.edu/jellyfish.html](http://www.genome.umd.edu/jellyfish.html)

### De novo assembly software
At one point we need assemblies of RenSeq data. We reccomend using CLC assembly cell ([https://www.qiagenbioinformatics.com/products/clc-assembly-cell/](https://www.qiagenbioinformatics.com/products/clc-assembly-cell/)), but free software, such as MaSuRCA ([http://www.genome.umd.edu/masurca.html](http://www.genome.umd.edu/masurca.html)) works as well.


### NLR Parser
Download NLR-Parser from [github.org/MutantHunter](https://github.com/steuernb/MutantHunter/releases/download/1/NLR-Parser.jar)

To run it, you will also need the [meme.xml](https://github.com/steuernb/MutantHunter/blob/master/meme.xml) containing the definitions of NLR associated motifs.

### R
For visualization, we use R. Donwload from [https://cran.r-project.org/](https://cran.r-project.org/)


## Pipeline

### 1. Preprocess reads (e.g. with Trimmomatic)

### 2. Count k-mers from read files for each accession

```
zcat accession1_R?.fastq.gz | jellyfish count -C -m 51 -s 3G -o accession1.jf /dev/fd/0
jellyfish dump -L 10 -ct accession1.jf > accession1.dump.txt
```

### 3. Create configuration file for Presence/Absense Matrix

This is a simple tab separated file with the accession names in the first column and the paths to the jellyfish dumps in the second column.

```
accession1	path/to/accession1.dump.txt
accession2	path/to/accession2.dump.txt
...
accessionN	path/to/accessionN.dump.txt
```

### 4. Create the Presence/Absense Matrix

```
java -jar AgRenSeq_CreatePresenceMatrix.jar -i accessions.txt -o AgRenSeq_k51_presencematrix.txt -t 3 -n 10
```

#### Parameters

Parameter | Argument | Description
--- | --- | ---
-i | accessions.txt | Mandatory. The path to the configuration file created in step 3.
-o | outputMatrix.txt | Mandatory. The path to the output file that will contain the matrix.
-n | integer | Default 10. The minimum kmer count for a k-mer to be considered present.
-t | integer | Default 3. A k-mer present in less accessions than _this value_ or present in all but _this value_ accessions will not be printed.

### 5. Create phenotype file

This is a tab separated file with accession names in the first colum. The following columns contain AgRenSeq scores. The recoreded score will be the average of all scores in one line. For AgRenSeq, the scores need to be negative for susceptible and positive for resistant. 

This is an example of the conversion for Stackman's IT (for wheat stem rust) to AgRenSeq scores.

Stackman's IT | AgRenSeq score
--- | ---
0 | 2; | 1.671- | 1.331 | 11+ | 0.672- | 0.332 | 02+ | -0.333- | -0.673 | -13+ | -1.334 | -2

### 6. _De novo_ assembly of resistant accession

Pick an accession where you expect _R_ genes to be (according to your phenotype). Run a de novo assembly on the RenSeq data. 

We have good experience with CLC assembly cell.  


### 7. Run NLR-Parser on assembly of resistant accession

This will select the contigs in the _de novo_ assembly that are associated with NLRs and in this way gets rid of off-target contigs.

```
java -jar NLR-Parser.jar -t <number of threads> -y <path/to/meme/bin/mast> -x <path/to/meme.xml> -i <sub-seqeunces.fasta> -o <output.nlr.txt>
```


### 9. Generate association scores of k-mers and project those onto the denovo assembly

This process will sum up the AgRenSeq scores from accessions where a k-mer is present and assigns the sum as an association score to a k-mer. In a second step, all association scores from k-mers within a contig from the _de novo_ assembly will be recorded in a tab separated file. For each contig, one line per unique association score is written as well as the number of k-mers that have been assigned with that score. Column 1 is the contig identifier, column 2 is a running number that increases with each contig, column 3 is the association score, column 4 is the number of k-mers in that contig that have been assigned with that score.

```
java -jar AgRenSeq_RunAssociation.jar -i prenseceMatrix -p phenotype -o AgRenSeqResult.txt

```


### 10. Plot dot-columns

Plot the result from step 9. using R. A simple script for R will look similar to this:

```
file<- read.table("AgRenSeqResult.txt", sep="\t")

v <- file$V1
x<-  file$V2
y <- file$V3
z<- file$V4

plot( x, y, pch = 20, cex = 0.5, main="TTKSK (Ug99)", ylab="score", xlab="NLR contigs BW_01077")
points( x, y, pch = 20, cex = 0.5)
points( x[z>25], y[z>25], pch = 20, cex = 1)
points( x[z>50], y[z>50], pch = 20, cex = 1.5)
points( x[z>70], y[z>70], pch = 20, cex = 2)
points( x[z>100], y[z>100], pch = 20, cex = 2.5)
points( x[z>125], y[z>125], pch = 20, cex = 3)
```










