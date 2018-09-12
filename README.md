# Michigan Imputation Helper

Here are two tools to 1) help prepare PLINK format genetic data for imputation using the [Michigan Imputation Server](https://imputationserver.sph.umich.edu) and 2) process the output back into PLINK format. 

MichPrep.r prepares PLINK files for uploading onto the Michigan Imputation Server. It carries out the Michigan Imputation Server [Data Preparation guidelines](https://imputationserver.readthedocs.io/en/latest/prepare-your-data/). It takes the genome-wide binary PLINK files, creates a frequency file, runs them through the recommended [Will Rayner authored scripts](http://www.well.ox.ac.uk/~wrayner/tools/), converts the PLINK files to per chromosome .vcf format, sorts and compresses the .vcf files, and then checks them using the recommended [CheckVCF.py](https://github.com/zhanxw/checkVCF) script.

PostMich.r downloads the imputed files from the Michigan Imputation Server and processes them back into PLINK format. The script downloads the files, decrypts them, converts them to PLINK format, extracts variants that are biallelic, have unique IDs and are well imputed (Rsq >= 0.3), and then merges the per chromosome files together.



## Getting Started with MichPrep.r

### Prerequisites

* Genome-wide genetic data in binary PLINK format (.bed/.bim/.fam) and that has undergone appropriate quality control.
* R and the required packages:

```R
install.packages(c('optparse','foreach','doMC'))
```

* [PLINK 1.9 software](https://www.cog-genomics.org/plink2)
* [perl](https://www.perl.org/get.html)
* [bcftools](https://samtools.github.io/bcftools/)
* Michigan Imputation Server preparation [script](http://www.well.ox.ac.uk/~wrayner/tools/) written by Will Rayner
* Michigan Imputation Server [CheckVCF.py](https://github.com/zhanxw/checkVCF) script

### Parameters

| Flag                | Description                                      |
| :------------------ | ------------------------------------------------ |
| --Input             | Prefix of genome-wide PLINK binaries             |
| --Input_dir         | Path to directory containing PLINK binaries      |
| --Output            | Prefix for output files                          |
| --Output_dir        | Directory for output files                       |
| --plink             | Path to plink software binary                    |
| --perl              | Path to perl software binary                     |
| --perl_library      | Path to perl library [optional]                  |
| --bcftools          | Path to bcftools software binary                 |
| --plink_checker     | Path to Will Rayner perl script                  |
| --plink_checker_ref | Path to reference data required by plink_checker |
| --vcf_checker       | Path to checkVCF.py                              |
| --vcf_checker_ref   | Path to reference data required by vcf_checker   |
| --ncores            | Number of cores available                        |
| --memory            | Available Memory in MB                           |

### Output

| Name    | Description                                                  |
| ------- | ------------------------------------------------------------ |
| .log    | Log file.                                                    |
| .vcf.gz | Per chromosome files formatted for Michigan Imputation Server |



## Getting Started with PostMich.r

### Prerequisites

- R and the required packages:

```r
install.packages('optparse')
```

- [PLINK 1.9 software](https://www.cog-genomics.org/plink2)
- File listing wget commands provided by the Michigan Imputation Server
- Password provided by the Michigan Imputation Server

### Parameters

| Flag          | Description                                                  |
| :------------ | ------------------------------------------------------------ |
| --Mich_wget   | Path to file listing wget commands from the Michigan Imputation Server website |
| --PostImp_dir | Path to directory where the Michigan Imputation Server files will be stored |
| --Output      | Prefix for output files                                      |
| --Output_dir  | Directory for processed PLINK files                          |
| --plink       | Path to plink software binary                                |
| --password    | Password to decrypt the genetic data                         |

### Output

The PostImp_dir will contain the files provided by the Michigan Imputation Server:

| Name           | Description                                                  |
| -------------- | ------------------------------------------------------------ |
| .zip           | Per chromosome password protected folder containing imputed genetic data in vcf format |
| .log           | Per chromosome log files from the imputation process         |
| qcreport.html  | Information on the imputation process                        |
| statistics.txt | SNP-statistics from the imputation process                   |

The Output_dir will contain:

| Name                            | Description                                              |
| ------------------------------- | -------------------------------------------------------- |
| -PostMich.log                   | Log file                                                 |
| .bed/.bim/.fam                  | Genome-wide PLINK binaries                               |
| -MichPrep-Non-biallelic.snplist | List of variants excluded because they weren't biallelic |
| -MichPrep-Duplicate.snplist     | List of variants excluded because they had duplicate IDs |




## Examples

##### Preparing data for HRC imputation:

```shell
qsub -l h_vmem=25G,mem_free=25G -b y Rscript MichPrep.r \
	--Input Test \
	--Input_dir /home/mpmop/Data \
	--Output Test-PreImp \
	--Output_dir /home/mpmop/Data-PreImp \
	--plink /share/apps/plink2 \
	--perl /opt/perl/bin/perl \
	--perl_library /opt/perl/lib/5.14.2 \
	--bcftools /home/mpmop/Software/bcftools/bcftools \
	--plink_checker /home/mpmop/Software/MichiganImp/HRC-1000G-check-bim-OPEdit.pl \
	--plink_checker_ref /home/mpmop/Software/MichiganImp/HRC.r1-1.GRCh37.wgs.mac5.sites.tab \
	--vcf_checker /home/mpmop/Software/MichiganImp/checkVCF.py \
	--vcf_checker_ref /home/mpmop/Software/MichiganImp/hs37d5.fa \
	--ncores 5 \
	--memory 25000
```

##### Downloading and processing the imputed data:

```shell
qsub -b y Rscript PostMich.r \
	--Mich_wget /home/mpmop/wget_commands.txt \
	--PostImp_dir /home/mpmop/Data-PostImp \
	--Output Test-PostImp-PreQC \
	--Output_dir /home/mpmop/Data-PostImp-PreQC \
	--plink /share/apps/plink2 \
	--password Bjhb785Gvkk7875
```



## Help

This script was written by Dr Oliver Pain under the supervision of Dr Richard Anney whilst at the MRC Centre for Neuropsychiatric Genetics and Genomics, Cardiff University.

If you have any questions or comments email me (paino@cardiff.ac.uk).







