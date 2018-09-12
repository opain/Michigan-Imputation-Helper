#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at Cardiff University under the supervision of Richard Anney.
suppressMessages(library("optparse"))

option_list = list(
make_option("--Mich_wget", action="store", default=NA, type='character',
		help="Path to file containing wget commands from the Michigan Imputation Server website [required]"),
make_option("--PostImp_dir", action="store", default=NA, type='character',
		help="Path to directory containing downloaded [required]"),
make_option("--Output", action="store", default=NA, type='character',
		help="Prefix for output files [required]"),
make_option("--Output_dir", action="store", default=NA, type='character',
		help="Directory for output files [required]"),
make_option("--plink", action="store", default=NA, type='character',
		help="Path to plink software binary [required]"),
make_option("--password", action="store", default=NA, type='character',
		help="Password to decrypt the genetic data [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

########################################
# 0. Create directory for the pre-imputation files to be stored in.
########################################

system(paste0('mkdir -p ',opt$PostImp_dir))
system(paste0('mkdir -p ',opt$Output_dir))
system(paste0('mkdir -p ',opt$Output_dir,'/temp'))

sink(file = paste0(opt$Output_dir,'/', opt$Output,'-PostMich.log'), append = F)
cat('#############################################
# Post Michigan Imputation Processing
#############################################

Specified Options:\n',sep='')

opt

sink()

######
# 1. Download the files
######
# The website provides a series of wget commands to download the results.
# Download them to same folders folders ending -PostImp

Mich_wget<-read.table(opt$Mich_wget, header=F, sep='*', stringsAsFactors=F)

sink(file = paste0(opt$Output_dir,'/', opt$Output,'-PostMich.log'), append = T)
cat('Downloading files from server...',sep='')
sink()

for(i in Mich_wget$V1){
	system(paste0(i,' -P ',opt$PostImp_dir))
}

sink(file = paste0(opt$Output_dir,'/', opt$Output,'-PostMich.log'), append = T)
cat('Done!\n',sep='')
sink()

######
# 2. Unzip the password protected genetic files.
######

genetic_files<-list.files(path = opt$PostImp_dir, pattern = '*.zip')
genetic_files<-gsub('chr_','',genetic_files)
genetic_files<-as.numeric(gsub('.zip','',genetic_files))

sink(file = paste0(opt$Output_dir,'/', opt$Output,'-PostMich.log'), append = T)
cat('Decrypting ',length(genetic_files),' files containing genetic data...',sep='')
sink()

for(i in genetic_files){
	system(paste0('7za x ',opt$PostImp_dir,'/chr_',i,".zip -p'",opt$password,"' -o",opt$Output_dir,'/temp/'))
}

sink(file = paste0(opt$Output_dir,'/', opt$Output,'-PostMich.log'), append = T)
cat('Done!\n',sep='')
sink()

######
# 3. Convert to plink and remove variants that are non-biallelic, duplicates, or have Rsq less than 0.3.
######

sink(file = paste0(opt$Output_dir,'/', opt$Output,'-PostMich.log'), append = T)
cat('Converting .vcf files to plink format...',sep='')
sink()

for(i in genetic_files){
	system(paste0('gunzip -c ',opt$Output_dir,'/temp/chr',i,'.info.gz > ',opt$Output_dir,'/temp/chr',i,'.info'))
	system(paste0(opt$plink,' --vcf ',opt$Output_dir,'/temp/chr',i,'.dose.vcf.gz --make-bed --out ',opt$Output_dir,'/temp/chr',i,'.init'))
	system(paste0(opt$plink,' --bfile ',opt$Output_dir,'/temp/chr',i,'.init --bmerge ',opt$Output_dir,'/temp/chr',i,'.init --merge-mode 6 --out ',opt$Output_dir,'/temp/chr',i))
	system(paste0(opt$plink,' --bfile ',opt$Output_dir,'/temp/chr',i,'.init --exclude ',opt$Output_dir,'/temp/chr',i,'.missnp --make-bed --out ',opt$Output_dir,'/temp/chr',i,'.bi'))
	system(paste0(opt$plink,' --bfile ',opt$Output_dir,'/temp/chr',i,'.bi --list-duplicate-vars ids-only suppress-first --out ',opt$Output_dir,'/temp/chr',i,'.bi'))
	system(paste0(opt$plink,' --bfile ',opt$Output_dir,'/temp/chr',i,'.bi --exclude ',opt$Output_dir,'/temp/chr',i,'.bi.dupvar --make-bed --out ',opt$Output_dir,'/temp/chr',i,'.bi.noDup'))
	system(paste0(opt$plink,' --bfile ',opt$Output_dir,'/temp/chr',i,'.bi.noDup --qual-scores ',opt$Output_dir,'/temp/chr',i,'.info 7 1 1 --qual-threshold 0.3 --make-bed --out ',opt$Output_dir,'/temp/chr',i,'.bi.noDup.Rsq3'))
}

sink(file = paste0(opt$Output_dir,'/', opt$Output,'-PostMich.log'), append = T)
cat('Done!\n',sep='')
sink()

# Create files containing variants excluded due to each exclusion criteria.
system(paste0('cat ',opt$Output_dir,'/temp/*.missnp > ',opt$Output_dir,'/',opt$Output,'-MichPrep-Non-biallelic.snplist'))
system(paste0('cat ',opt$Output_dir,'/temp/*.dupvar > ',opt$Output_dir,'/',opt$Output,'-MichPrep-Duplicate.snplist'))

# Count how many variants were were to start with
N_SNP_init<-system(paste0('wc -l ',opt$Output_dir,'/temp/chr*.init.bim'), intern=T)
N_SNP_init<-N_SNP_init[length(N_SNP_init)]
N_SNP_init<-gsub(' ', '', N_SNP_init)
N_SNP_init<-as.numeric(gsub('total', '', N_SNP_init))

# Count how many SNPs were non-biallelic
N_SNP_bi<-system(paste0('wc -l ',opt$Output_dir,'/temp/chr*.bi.bim'), intern=T)
N_SNP_bi<-N_SNP_bi[length(N_SNP_bi)]
N_SNP_bi<-gsub(' ', '', N_SNP_bi)
N_SNP_bi<-as.numeric(gsub('total', '', N_SNP_bi))
N_SNP_nonBi<-N_SNP_init-N_SNP_bi

# Count how many SNPs were removed due to duplication
N_SNP_bi_noDup<-system(paste0('wc -l ',opt$Output_dir,'/temp/chr*.bi.noDup.bim'), intern=T)
N_SNP_bi_noDup<-N_SNP_bi_noDup[length(N_SNP_bi_noDup)]
N_SNP_bi_noDup<-gsub(' ', '', N_SNP_bi_noDup)
N_SNP_bi_noDup<-as.numeric(gsub('total', '', N_SNP_bi_noDup))
N_SNP_bi_Dup<-N_SNP_bi-N_SNP_bi_noDup

# Count how many SNPs were removed due to low Rsq
N_SNP_bi_noDup_Rsq3<-system(paste0('wc -l ',opt$Output_dir,'/temp/chr*.bi.noDup.Rsq3.bim'), intern=T)
N_SNP_bi_noDup_Rsq3<-N_SNP_bi_noDup_Rsq3[length(N_SNP_bi_noDup_Rsq3)]
N_SNP_bi_noDup_Rsq3<-gsub(' ', '', N_SNP_bi_noDup_Rsq3)
N_SNP_bi_noDup_Rsq3<-as.numeric(gsub('total', '', N_SNP_bi_noDup_Rsq3))
N_SNP_bi_lowRsq<-N_SNP_bi_noDup-N_SNP_bi_noDup_Rsq3

sink(file = paste0(opt$Output_dir,'/', opt$Output,'-PostMich.log'), append = T)
cat('Data originally contained ',N_SNP_init,' variants.\n', sep='')
cat(N_SNP_nonBi," variants were not biallelic (",opt$Output_dir,'/',opt$Output,"-MichPrep-Non-biallelic.snplist).\n", sep='')
cat(N_SNP_bi_Dup," variants did not have unique IDs i.e.chr:bp (",opt$Output_dir,'/',opt$Output,"-MichPrep-Duplicate.snplist).\n", sep='')
cat(N_SNP_bi_lowRsq," variants had an Rsq < 0.3.\n", sep='')
cat('After exclusions ',N_SNP_bi_noDup_Rsq3,' variants remain.\n',sep='')
sink()

######
# 4. Merge per chromosome files
######

sink(file = paste0(opt$Output_dir,'/', opt$Output,'-PostMich.log'), append = T)
cat('Merging per chromosome files...',sep='')
sink()

# Create merge.list.
system(paste0('ls ',opt$Output_dir,"/temp/*.Rsq3.fam | sed 's/.fam//g' > ",opt$Output_dir,'/temp/merge.list'))

# Merge the per chromosome files.
system(paste0(opt$plink,' --merge-list ',opt$Output_dir,'/temp/merge.list --make-bed --out ',opt$Output_dir,'/',opt$Output))

sink(file = paste0(opt$Output_dir,'/', opt$Output,'-PostMich.log'), append = T)
cat('Done!\n',sep='')
sink()

# Delete temporary files.
system(paste0('rm -r ',opt$Output_dir,'/temp'))

sink(file = paste0(opt$Output_dir,'/', opt$Output,'-PostMich.log'), append = T)
cat('Finished!\n',sep='')
sink()









