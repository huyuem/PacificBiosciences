# pacbio_isoseq_pipeline
This is pipeline wrapping different steps in PacBio Iso-Seq analysis from RS II cell (i.e. A01_1) to HTML reports ready to be presented.<br>
It is mainly developped for internal use in pacBio and has functionalities overlapping with [isoaux](https://github.com/bowhan/isoaux). <br>

## Install
The following software/tools are needed: <br>

-  **Linux enviroment with smrtanalysis_2.3.0.140936 installed**<br>
Specify the home of smrtanalysis in the `config/config.sh` file, such as 
```bash
declare -x SMRT_HOME='/pbi/analysis/smrtportal/beta/install/smrtanalysis_2.3.0.140936'
```

-  **SGE job submission envirment**<br>
Specify the job submission command in the `config/config.sh` file, such as
```bash
declare -x SUBMIT_CMD='qsub -V -cwd -pe smp 8 -l h_rt=240:0:0,mem_free=8G'
# -V: load enviromental variable
# -cwd: use current working directory 
# -pe smp 8: use smp 8 as the parallel enviroment; use `qconf -spl` to list all parallel enviroment
# -l: limitations, consult your SGE admin; a typical job runs overnight 
```
*Currently only SGE is supported*

- User specific directories set-up
Edit `config/config.sh` to specify the directory related enviromental variables, such as <br> 
```bash
# R libraries will be installed in this directory if not already exist
declare -x RLIBDIR='/home/bhan/Rlib'
# the directory storing genome sequence, index, and annotation files
declare -x ANNOTATION_DIR='/home/bhan/annotation'
```

-  **Python3**<br>
Packages used: biopython
```bash
# to install
pip3 install biopython
```

-  **R**<br>
Relatively new version; missing packages will be installed automatically during the pipeline run.<br>
But it is highly recommended to install them beforehand. 
```R
# to install those packages
install.packages(c("readr", "dplyr", "tidyr", "rmarkdown", "corrplot", "ggplot2", "ggthemes", "plotly", "shiny"))
```

- **MySQL**<br>
If using the `annotation` pipeline to prepare for annotation files, MySQL is required to obtain information from UCSC.

- Other tools

[trim_isoseq_polyA](https://github.com/bowhan/trim_isoseq_polyA)

[mrna_size_from_gff and colmerge](https://github.com/bowhan/isoaux)

[deepTools](https://github.com/fidelram/deepTools)

[picard](https://broadinstitute.github.io/picard/)

[kent tools](https://github.com/ENCODE-DCC/kentUtils)

[gmap and gmap_build](http://research-pub.gene.com/gmap/)

## Usage
```bash
# 1. Preparation
# edit config file to specify your smrt analysis home and qsub enviroment
vi pacbio_isoseq_pipeline/config/config.sh 
# prepare genome sequence, annotation files using annotation pipeline
isoseq.sh anno -g hg19
isoseq.sh anno -g hg38
isoseq.sh anno -g mm10
isoseq.sh anno -g dm3
# the pipeline submits jobs to build gmap index, so please wait until those jobs finish before running analysis!

# 2. Prepare the sample information
# copy example.csv 
cp pacbio_isoseq_pipeline/sample/example.csv my_isoseq_sample.csv
# edit my_isoseq_sample.csv by filling it in your sample
vi my_isoseq_sample.csv
```
The user will need to fill the CSV file which specify the information needed, including:

1.  genome<br>
The name of the genome whose annotation has been obtained using the `annotation` pipeline.<br>
Specifically, the genome fasta and genes.gtf files need to be presented in the `${ANNOTATION_DIR}` folder, which is defined in `config/config.sh`. <br>
For example, if `hg19` is used, the pipeline will expect two files called `hg19.fa` and `hg19.genes.gtf` inside of the `${ANNOTATION_DIR}` folder.<br>
Additionally, a gmap index can be put in the ``${ANNOTATION_DIR}`/index/gmap_index` folder. <br>

- The "annotation" pipeline `isoseq.sh anno` prepares those files for popular genome assemblies using UCSC srouce.<br>

2. tissue<br>
Libraries of different tissues are compared separately. <br>

3. treatment<br>
Libraries of different treatment are compared directly.<br>

4. size bin<br>
The size bin used when the libraries are generated, like "1-3". <br>
Only libraries in the same size bins are compared directly. <br>

5. Path to the RS II cells<br>
The complete, ABSOLUTE path to the RS II cells (usually named [A-Z]0[12]_1). <br>
At least three `bax.h5` files are expected in the `Analysis_Results` subdirectory.<br>

6. Path to the custom barcode sequence file<br>
The absolute path to a fasta file with custom barcode sequence. A barcode fasta file looks like:<br>
```
>F0
GCAGTCGAACATGTAGCTGACTCAGGTCACTCAGACGATGCGTCATGGTAGTAAGCAGTGGTATCAACGCAGAGTAC
>R0
GTACTCTGCGTTGATACCACTGCTTACTACCATGACGCATCGTCTGAGTGACCTGAGTCAGCTACATGTTCGACTGC
>F1
GCAGTCGAACATGTAGCTGACTCAGGTCACCTGCGTGCTCTACGACGGTAGTAAGCAGTGGTATCAACGCAGAGTAC
>R1
GTACTCTGCGTTGATACCACTGCTTACTACCGTCGTAGAGCACGCAGGTGACCTGAGTCAGCTACATGTTCGACTGC
>F2
GCAGTCGAACATGTAGCTGACTCAGGTCACCATAGCGACTATCGTGGGTAGTAAGCAGTGGTATCAACGCAGAGTAC
>R2
GTACTCTGCGTTGATACCACTGCTTACTACCCACGATAGTCGCTATGGTGACCTGAGTCAGCTACATGTTCGACTGC
```
With forward and reverse primer intervened. If default primers are used (no barcode), fill `NA`. <br>

7. Which barcode is used<br>
The position of the barcode inside of the custom barcode fasta file.<br> 
**Start counting from 0**.
It has nothing to do with the names of PacBio barcode products. <br>
If no barcode is used, fill `NA`.<br>

-   Please use the template in the `sample/example.csv` file. Do **NOT** change the header.<br>
-   Avoid using space, `.`, `-`. Use **underscore `_`** to separate words.

```bash
# run the pipeline
isoseq.sh all -c my_isoseq_sample.csv -J myjobname -E my@emailaddress.com
```

