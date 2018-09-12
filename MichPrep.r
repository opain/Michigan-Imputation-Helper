#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at Cardiff University under the supervision of Richard Anney.
suppressMessages(library("optparse"))

option_list = list(
make_option("--Input", action="store", default=NA, type='character',
		help="Prefix of genome-wide PLINK binaries [required]"),
make_option("--Input_dir", action="store", default=NA, type='character',
		help="Path to directory containing PLINK binaries [required]"),
make_option("--Output", action="store", default=NA, type='character',
		help="Prefix for output files [required]"),
make_option("--Output_dir", action="store", default=NA, type='character',
		help="Directory for output files [required]"),
make_option("--plink", action="store", default=NA, type='character',
		help="Path to plink software binary [required]"),
make_option("--perl", action="store", default=NA, type='character',
		help="Path to perl software binary [required]"),
make_option("--perl_library", action="store", default=NA, type='character',
		help="Path to perl library [optional]"),
make_option("--bcftools", action="store", default=NA, type='character',
		help="Path to bcftools software binary [required]"),
make_option("--plink_checker", action="store", default=NA, type='character',
		help="Path to Will Rayner perl script [required]"),
make_option("--plink_checker_ref", action="store", default=NA, type='character',
		help="Path to reference data required by plink_checker [required]"),
make_option("--vcf_checker", action="store", default=NA, type='character',
		help="Path to checkVCF.py [required]"),
make_option("--vcf_checker_ref", action="store", default=NA, type='character',
		help="Path to reference data required by vcf_checker [required]"),
make_option("--ncores", action="store", default=NA, type='numeric',
		help="Number of cores available [required]"),
make_option("--memory", action="store", default=NA, type='numeric',
		help="Available Memory in MB [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

opt$temp<-paste0(opt$Output_dir,'/temp')

library(foreach)
library(doMC)
registerDoMC(opt$n_cores)

########################################
# 0. Create directory for the pre-imputation files to be stored in.
########################################

system(paste0('mkdir -p ',opt$Output_dir))
system(paste0('mkdir -p ',opt$temp))

sink(file = paste0(opt$Output_dir,'/', opt$Output,'MichPrep.log'), append = F)
cat('#############################################
# Michigan Imputation Server Preparation
#############################################

Specified Options:\n',sep='')

sink()

########################################
# 1. Check the bim files against the HRC imputation panel using William Rayner's perl script.
########################################

# Change the working directory
setwd(opt$temp)

# Create a frequency file for each dataset
sink(file = paste0(opt$Output_dir,'/', opt$Output,'MichPrep.log'), append = T)
paste0('Creating .frq file...')
sink()
system(paste0(opt$plink,' --freq --bfile ',opt$Input_dir,'/',opt$Input,' --out ',opt$Output))

# Run perl script
sink(file = paste0(opt$Output_dir,'/', opt$Output,'MichPrep.log'), append = T)
cat('Running plink_checker script...',sep='')
sink()

if(!is.na(opt$perl_library)){
	system(paste0('export PERL5LIB=',opt$perl_library))
}
system(paste0(opt$perl,' ',opt$plink_checker,' -b ',opt$Input_dir,'/',opt$Input,'.bim -f ',opt$Output,'.frq -r ', opt$plink_checker_ref,' -h'))

sink(file = paste0(opt$Output_dir,'/', opt$Output,'MichPrep.log'), append = T)
cat('Done!\n',sep='')
sink()

# Edit R-plink.sh file to match the location of the files.
# Run plink script
plink_script<-read.table('Run-plink.sh', sep='*', header=F, stringsAsFactors=F)
plink_script$V1[1]<-paste0('plink --bfile ',opt$Input_dir,'/',opt$Input,' --exclude Exclude-',opt$Input,'-HRC.txt --make-bed --out TEMP1')
plink_script$V1<-gsub('plink', opt$plink, plink_script$V1)
write.table(plink_script,'Run-plink-OPedit.sh', col.names=F, row.names=F, quote=F) 
system('sh Run-plink-OPedit.sh')

########################################
# 2. Convert genome-wide PLINK files into per chromosome, bp sorted, compressed vcf files.
########################################

# Convert to vcf format
sink(file = paste0(opt$Output_dir,'/', opt$Output,'MichPrep.log'), append = T)
cat('Converting plink files to vcf...',sep='')
sink()

for(i in 1:22){
	system(paste0(opt$plink,' --bfile ',opt$Input,'-updated-chr',i,' --recode vcf --out ',opt$Input,'-updated-chr',i))
}

sink(file = paste0(opt$Output_dir,'/', opt$Output,'MichPrep.log'), append = T)
cat('Done!\n',sep='')
sink()

# Sort and compress vcf files
sink(file = paste0(opt$Output_dir,'/', opt$Output,'MichPrep.log'), append = T)
cat('Sorting and compressing vcf files...',sep='')
sink()

setwd(opt$Output_dir)
foreach(i=1:22, .combine=c) %dopar% {
	system(paste0(opt$bcftools,' sort temp/',opt$Input,'-updated-chr',i,'.vcf -Oz -o ',opt$Output,'-chr',i,'.vcf.gz --max-mem ',floor((opt$memory*1000*0.8)/opt$ncores)))
}

sink(file = paste0(opt$Output_dir,'/', opt$Output,'MichPrep.log'), append = T)
cat('Done!\n',sep='')
sink()

########################################
# 3. Check compressed vcf files using the Michigan Imputation Server script checkVCF.py.
########################################

# Run checkVCF.py across each chromosome.
sink(file = paste0(opt$Output_dir,'/', opt$Output,'MichPrep.log'), append = T)
cat('Running vcf_checker script...',sep='')
sink()

foreach(i=1:22, .combine=c) %dopar% {
	system(paste0(opt$vcf_checker,' -r ',opt$vcf_checker_ref,' -o temp/',opt$Output,'-chr',i,' ',opt$Output,'-chr',i,'.vcf.gz'))
}

sink(file = paste0(opt$Output_dir,'/', opt$Output,'MichPrep.log'), append = T)
cat('Done!\n',sep='')
sink()

# Check there are no invalid genotypes
sink(file = paste0(opt$Output_dir,'/', opt$Output,'MichPrep.log'), append = T)
for(i in 1:22){
	log<-read.table(paste0('temp/',opt$Output,'-chr',i,'.check.log'), header=F,stringsAsFactors=F, sep='*')
	grepl('invalid genotypes', log$V1)
	n_inval<-as.numeric(unlist(strsplit(log$V1[grepl('invalid genotypes', log$V1)],' '))[4])
	if(n_inval > 0){
		cat('Chromosome ', i,' contains ', n_inval,' invalid genotypes. See ',unlist(strsplit(log$V1[18],' '))[15],'.\n',sep='')
		error<-1
	}
}
sink()

##########################################
# 4. Delete temporary files
##########################################

if(!exists('error')){
	system('rm -r temp')
}

sink(file = paste0(opt$Output_dir,'/', opt$Output,'MichPrep.log'), append = T)
cat('Finished!',sep='')
sink()
