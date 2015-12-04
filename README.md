# IPED
A highly efficient denoising tool for Illumina paired-end 16S rRNA amplicon sequencing data.
The development of high-throughput sequencing technologies has revolutionized the field of 16S RNA amplicon sequencing approaches, aimed at assessing microbial diversity. An important step in preprocessing these amplicon sequencing data is the removal of sequencing errors, often referred to as denoising. However most implementations are primarily focused on denoising 454 pyrosequencing data. We have developed a dedicated algorithm (IPED) for correcting sequencing errors in paired-end reads obtained from the Illumina MiSeq platform. The machine learning method developed in this work is able to predict positions in the sequencing reads potentially containing errors. Subsequently this information is used to group those error-prone reads with correct reads resulting in error-free consensus reads by masking potentially erroneous positions during this clustering. Benchmarking our algorithm on real-life mock data showed an average improvement of 50% compared with the second best algorithm, and for amplicons where paired read are not completely overlapping almost an order of magnitude reduction of the error rate has been observed. Reducing the error rate had a positive effect on the clustering of reads in operational taxonomic units, with an almost perfect correlation between the number of clusters and the number of species present in the mock communities. The IPED software can be downloaded here (http://science.sckcen.be/~/media/Files/Science/EHS/Bioinformatics_expertise/IPED/IPED_unix.zip?la=en) and downloaded here (http://science.sckcen.be/~/media/Files/Science/EHS/Bioinformatics_expertise/IPED/IPED_mac.zip?la=en).
# Installation Requirement
Only Perl needs to be installed to run IPED. All other software packages that are required to run IPED are included in the downloaded file (IPED.run). In case you are interested in the source code of IPED, this is also included in the downloaded file. Only in case you want to run the source code, you will need to install those software components separately, and adapt the source code referring to those software components accordingly. In all other cases, we encourage the end-user to use the IPED.run command.

# Included Software
Software listed below is used by the IPED algorithm. However you do NOT need to install it separately as these software modules are included in the IPED software.

Mothur v.1.33.3:
Available at http://www.mothur.org/wiki/Download_mothur. 
Note about changes made in the mothur package integrated in this package:
The command called "pre.cluster" is modified to be compatible with the IPED algorithm.
The command called "make.contigs" is modified to produce an additional IPED-formatted quality file.
WEKA 3.7.11: 
Available online at http://www.cs.waikato.ac.nz/ml/weka/.
Both WEKA and mothur are distributed under the GNU licence.
# Syntax:
!!! make sure you use an underscore (and not a hyphen) to specify the option you want to set.
!!! make sure to use the complete path when describing files
When running IPED the first time on your sample, you are obliged to first run the IPED version of make.contigs in order to produce an IPED adapted quality file.  IPED can be run following two different approaches:

In two consecutive commands: first running the IPED version of make.contigs, and secondly run the denoising (IPED) step.
in one single command, combining the IPED adapted make.contigs command and the denoising (IPED) step.
The description below follows the first approach, but an example of how the command line looks like in case of the second approach, is given in the "testing" section.

Step 1: Creating contigs

Mandatory Options:
_F Forward fastq
_R Reverse fastq
_o Output Path.
Non mandatory Options
_D Deltaq default =6 (check mothur documentation for more info)
_I Insert default = 20 (check mothur documentation for more info)
Example command: ./IPED.run _F ./sample.forward.fastq _R ./sample.reverse.fastq _o OUTPUT_PATH

Step 2: Denoising (IPED)

Mandatory Options:
_n Name file with the redundancy
Tab separated file with the read ID on the first column and the IDs of identical reads on the second column [separated by comma]).
Can be produced after Dereplication of the reads via the mothur command [unique.seqs(fasta = file.fasta)]

_f Aligned sequences file
Aligned multi-fasta files without any ambiguity, only "AGTC" bases (can be created using the align.seqs and filter.seqs commands in mothur).

_o Output Path              
_c Contigs fasta file
Fasta file created in the first step (make.contigs).

_q Qual file of contigs
Quality file created in the first step (make.contigs)

Non mandatory Options
_p number of processors, [default 1]
_g group file (in case of having different sample-groups)
 Tab separated file, with read ID in the first column and read sample (group) name on the second, as produced by mothur.

_d differences tolerated
[Default: will be automatically calculated with a cut-off leading to clustering of reads with 98% similarity. For example, a length of 400 bp will results in a default difference cutoff of 4 (it is recommended to increase the differences by "1" each additional 100 bp in the read average length i.e. 98% similarity). Check mothur for more information.

_i log ID
fill in to continue a previous run. If left empty, a new run will be created with a random number (default: random number)
Example command: ./IPED.run _f ./sample.fasta _n ./sample.names _c /your/complete/path/sample.trim.contigs.fasta _q /your/complete/path/sample.contigs.qual _o OUTPUT_PATH
To see the program help, type IPED.run (with no parameters).

# Output Files
The IPED program generates different text output files distributed over two folders "IPED_Final" and "IPED_Temp"
Inside each of them another folder can be found, having as name in both folders the same random number (i.e. the process ID)

IPED_Final:
Contains the final alignment and name file after denoising and preclustering named:

  Results.IPED.names
  Results.IPED.fasta
  In case of using _F & _R options:
    Forward_raw_fastq_name.trim.contig.fasta
    Forward_raw_fastq_name.contig.qual
IPED_Temp
Split-up alignment files 
  (0.fasta, 1.fasta ,etc) and name files (0.names ,1.names, etc) according to the number of processors used.
IPED model (classifier) input and output named 0.Test.arff, 0.Test.Final, etc respectively.
Split-up alignment files with the marked (i.e. potentially erroneous) positions (0.Results, 1.Results, etc).
Merged aligment file (Results.fasta) and names (Results.names)
output from modified Pre-cluster algorithm: Results.precluster.fasta, Result.precluster.names
Final output after converting the masked positions back to their original nucleotide (Results.precluster.pick.fasta, Results.precluster.pick.names)
logfile
The program produces a log file of the steps being run, making it possible to monitor the different steps, and track down possible errors:
IPED: OUTPUT_PATH/IPED_Temp/****/IPED_****.logfile 
(**** represents a random process ID).

# Testing
Type:

IPED.run _f ./sample.fasta _n ./sample.names _F ./sample.forward.fastq _R ./sample.reverse.fastq _o OUTPUT_PATH
The different input file (sample.fasta, sample.names, sample.forward.fastq, sample.reverse.fastq) are included in the software. It will produce in the output path a file containing the results i.e. two files named Sample.IPED.fasta and Sample.IPED.names. The command given above integrates the two steps (running make.contigs and running IPED) in one command. In case you want to re-analyze a sample, and contigs have already been created for that sample, you can omit the make.contig step and input directly in the command line the contig-file and the IPED adapted quality file. The Test dataset is a part of the data published in Kozich et al. (2013)  and downloaded from http://www.mothur.org/MiSeqDevelopmentData.html.
# Citing:
If you are going to use IPED, please cite it with the included software (mothur, WEKA):

M.Mysara, J. Raes, N.Leys, P.Monsieurs (2014) IPED: A highly efficient denoising tool for Illumina paired-end 16S rRNA amplicon sequencing data.
Schloss PD, Westcott SL, Ryabin T, Hall JR, Hartmann M, Hollister EB, et al. (2009). Introducing mothur: open-source, platform-independent, community-supported software for describing and comparing microbial communities. Applied and environmental microbiology 75:7537–41.
Hall M, National H, Frank E, Holmes G, Pfahringer B, Reutemann P, et al. (2009). The WEKA Data Mining Software?: An Update. SIGKDD Explorations 11:10–18.
Contact:
For questions, bugs and suggestions, please refer to mohamed.mysara@sckcen.be & pieter.monsieurs@sckcen.be
Developed by M.Mysara et al. 2014
