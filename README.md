# pacbio_isoseq_pipeline
This is a bash based pipeline to wrap all the different steps of PacBio Iso-Seq pipeline from RS II cell (like A01_1) to reports ready to be presented.

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

-  **Python3**<br>
Packages used: biopython
```bash
# to install
pip3 install biopython
```

- User specific directories set-up
Edit `config/config.sh` to specify the directory related enviromental variables, such as <br> 
```bash
# R libraries will be installed in this directory if not already exist
declare -x RLIBDIR='/home/bhan/Rlib'
# the user will need to put genome sequence and annotations in this directory
declare -x ANNOTATION_DIR='/home/bhan/annotation'
```

-  **R**<br>
Relatively new version; missing packages will be installed automatically

- trim_isoseq_polyA<br>
You can obtain it from https://github.com/bowhan/trim_isoseq_polyA 

## Usage
#### The user will need to fill a CSV file which specify the information needed, including:

1.  genome<br>
The genome fasta and genes.gtf files need to be presented in the `${ANNOTATION_DIR}` folder, which is defined in `config/config.sh`. <br>
For example, if `hg19` is used, the pipeline will expect two files called `hg19.fa` and `hg19.genes.gtf` inside of the `${ANNOTATION_DIR}` folder.<br>
Additionally, a gmap index can be put in the ``${ANNOTATION_DIR}`/index/gmap_index` folder. If it is not present, the pipeline will build the index on the first run using this genome. <br>

2. tissue<br>
Libraries of different tissues are compared separately. <br>

3. treatment<br>
Libraries of different treatment are compared directly.<br>

4. size bin<br>
The size bin used when the libraries are generated, like "1-3". <br>
Only libraries in the same size bins are compared directly. <br>

5. Path to the RS II cells<br>
The complete, absolute path to the cells (usually named [A-Z]0[12]_1). <br>
At least three `bax.h5` files are expected in the `Analysis_Results` subdirectory.<br>

6. Path to the custom barcode sequence file<br>
The absolute path to a fasta file with custom barcode sequence.<br>
If default primers are used, use `NA`. <br>

7. Which barcode is used<br>
The position of the barcode inside of the custom barcode fasta file.<br> 
**Start counting from 0**.
It has nothing to do with the names of PacBio barcode products. <br>
If no barcode is used, use `NA`.<br>

-   Please use the template in the `sample/example.csv` file. Do **NOT** change the header.<br>
-   Avoid using space, `.`, `-`. Use **underscore `_`** to separate words.

#### Invoke the pipeline

```bash
# edit config file to specify your smrt analysis home and qsub enviroment
vi pacbio_isoseq_pipeline/config/config.sh 

# link the genome sequence and transcriptome annotation file
# I recommend to use the UCSC annotation from iGenome
cd pacbio_isoseq_pipeline/annotation
ln -s /path/to/your/annotation/dir/hg19.fa hg19.fa
ln -s /path/to/your/annotation/dir/hg19.genes.gtf hg19.genes.gtf
# now you can use "hg19" in the sample.csv file 

# copy example.csv 
cp pacbio_isoseq_pipeline/sample/example.csv my_isoseq_sample.csv
# edit my_isoseq_sample.csv by filling it in your sample
vi my_isoseq_sample.csv

# run the pipeline
pacbio_isoseq_pipeline/isoseq.sh -c my_isoseq_sample.csv -J myjobname -E my@emailaddress.com
```
