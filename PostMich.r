#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at Cardiff University under the supervision of Richard Anney.
start.time <- Sys.time()
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
		help="Password to decrypt the genetic data [required]"),
make_option("--HRC_ref", action="store", default=NA, type='character',
		help="Directory to HRC.r1-1.GRCh37.wgs.mac5.sites.tab [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

suppressMessages(library(data.table))

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

if(is.na(opt$Mich_wget) == T){
	cat('Error: --Mich_wget is not specified.\n.')
	q()
}
if(is.na(opt$PostImp_dir) == T){
	cat('--PostImp_dir is not specified.\n.')
	q()
}
if(is.na(opt$Output) == T){
	cat('--Output is not specified.\n.')
	q()
}
if(is.na(opt$Output_dir) == T){
	cat('--Output_dir is not specified.\n.')
	q()
}
if(is.na(opt$plink) == T){
	cat('--plink is not specified.\n.')
	q()
}
if(is.na(opt$password) == T){
	cat('--password is not specified.\n.')
	q()
}
if(is.na(opt$HRC_ref) == T){
	cat('--HRC_ref is not specified.\n.')
	q()
}

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
cat('Converting .vcf files to plink format and inserting RSIDs...',sep='')
sink()

for(i in genetic_files){
	# Decompress info file
	system(paste0('gunzip -c ',opt$Output_dir,'/temp/chr',i,'.info.gz > ',opt$Output_dir,'/temp/chr',i,'.info'))
	
	# Convert to vcf format
	system(paste0(opt$plink,' --vcf ',opt$Output_dir,'/temp/chr',i,'.dose.vcf.gz --make-bed --keep-allele-order --out ',opt$Output_dir,'/temp/chr',i,'.init'))
	
	# Merge with self and remove any variants that have more than 2+ alleles.
	system(paste0(opt$plink,' --bfile ',opt$Output_dir,'/temp/chr',i,'.init --bmerge ',opt$Output_dir,'/temp/chr',i,'.init --merge-mode 6 --out ',opt$Output_dir,'/temp/chr',i,'.init'))
	system(paste0(opt$plink,' --bfile ',opt$Output_dir,'/temp/chr',i,'.init --exclude ',opt$Output_dir,'/temp/chr',i,'.init.missnp --make-bed --out ',opt$Output_dir,'/temp/chr',i,'.init.bi'))
	
	# Remove any variants with duplicate IDs.
	system(paste0(opt$plink,' --bfile ',opt$Output_dir,'/temp/chr',i,'.init.bi --list-duplicate-vars ids-only suppress-first --out ',opt$Output_dir,'/temp/chr',i,'.init.bi'))
	system(paste0(opt$plink,' --bfile ',opt$Output_dir,'/temp/chr',i,'.init.bi --exclude ',opt$Output_dir,'/temp/chr',i,'.init.bi.dupvar --make-bed --out ',opt$Output_dir,'/temp/chr',i,'.init.bi.noDup'))
	
	# Extract variants with an R2 greater than 0.8.
	system(paste0(opt$plink,' --bfile ',opt$Output_dir,'/temp/chr',i,'.init.bi.noDup --qual-scores ',opt$Output_dir,'/temp/chr',i,'.info 7 1 1 --qual-threshold 0.8 --make-bed --out ',opt$Output_dir,'/temp/chr',i,'.init.bi.noDup.Rsq3'))

	# Update the variant IDs to RSIDs where possible.
	bim<-data.frame(fread(paste0(opt$Output_dir,'/temp/chr',i,'.init.bi.noDup.Rsq3.bim')))
	names(bim)<-paste0('bim_',names(bim))
	system(paste0("awk '$1 == \"",i,"\" { print $1\" \"$2\" \"$3 }' ",opt$HRC_ref," > ", opt$Output_dir,"/temp/HRC_ref_",i))
	HRC_ref<-data.frame(fread(paste0(opt$Output_dir,"/temp/HRC_ref_",i)))
	names(HRC_ref)<-paste0('HRC_ref_',names(HRC_ref))
	HRC_ref$HRC_ref_ID<-paste(HRC_ref[,1],HRC_ref[,2],sep=':')
	bim_HRC_ref<-merge(bim,HRC_ref, by.x='bim_V2', by.y='HRC_ref_ID')
	if(dim(bim_HRC_ref)[1] != dim(bim)[1]){
		sink(file = paste0(opt$Output_dir,'/', opt$Output,'-PostMich.log'), append = T)
		cat("Error when updating to RSIDs for chr",i,". Reference and target do not merge 1:1.\n",sep='')
		sink()
	}
	ID_update<-bim_HRC_ref[c('bim_V2','HRC_ref_V3')]
	ID_update<-ID_update[ID_update$HRC_ref_V3 != '.',]
	fwrite(ID_update, paste0(opt$Output_dir,"/temp/ID_update.",i,".txt"), col.names=F, sep=' ')
	system(paste0(opt$plink,' --bfile ',opt$Output_dir,'/temp/chr',i,'.init.bi.noDup.Rsq3 --update-name ',opt$Output_dir,"/temp/ID_update.",i,".txt --make-bed --out ",opt$Output_dir,'/temp/chr',i,'.init.bi.noDup.Rsq3.rsid'))
	rm(bim,HRC_ref,bim_HRC_ref,ID_update)
}

sink(file = paste0(opt$Output_dir,'/', opt$Output,'-PostMich.log'), append = T)
cat('Done!\n',sep='')
sink()

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
N_SNP_noDup<-system(paste0('wc -l ',opt$Output_dir,'/temp/chr*.noDup.bim'), intern=T)
N_SNP_noDup<-N_SNP_noDup[length(N_SNP_noDup)]
N_SNP_noDup<-gsub(' ', '', N_SNP_noDup)
N_SNP_noDup<-as.numeric(gsub('total', '', N_SNP_noDup))
N_SNP_Dup<-N_SNP_bi-N_SNP_noDup

# Count how many SNPs were removed due to low Rsq
N_SNP_Rsq3<-system(paste0('wc -l ',opt$Output_dir,'/temp/chr*.Rsq3.bim'), intern=T)
N_SNP_Rsq3<-N_SNP_Rsq3[length(N_SNP_Rsq3)]
N_SNP_Rsq3<-gsub(' ', '', N_SNP_Rsq3)
N_SNP_Rsq3<-as.numeric(gsub('total', '', N_SNP_Rsq3))
N_SNP_lowRsq<-N_SNP_noDup-N_SNP_Rsq3

# Count how many SNPs RSIDs available
N_SNP_RSID<-system(paste0('wc -l ',opt$Output_dir,"/temp/ID_update.*.txt"), intern=T)
N_SNP_RSID<-N_SNP_RSID[length(N_SNP_RSID)]
N_SNP_RSID<-gsub(' ', '', N_SNP_RSID)
N_SNP_RSID<-as.numeric(gsub('total', '', N_SNP_RSID))

# Create files containing variants excluded due to each exclusion criteria.
if(N_SNP_nonBi > 0){
	system(paste0('cat ',opt$Output_dir,'/temp/*.missnp > ',opt$Output_dir,'/',opt$Output,'-MichPrep-Non-biallelic.snplist'))
}
if(N_SNP_Dup > 0){
	system(paste0('cat ',opt$Output_dir,'/temp/*.dupvar > ',opt$Output_dir,'/',opt$Output,'-MichPrep-Duplicate.snplist'))
}

sink(file = paste0(opt$Output_dir,'/', opt$Output,'-PostMich.log'), append = T)
cat('Data originally contained ',N_SNP_init,' variants.\n', sep='')
if(N_SNP_nonBi > 0){
	cat(N_SNP_nonBi," variants were not biallelic (",opt$Output_dir,'/',opt$Output,"-MichPrep-Non-biallelic.snplist).\n", sep='')
}
if(N_SNP_Dup > 0){
	cat(N_SNP_Dup," variants did not have unique IDs i.e.chr:bp:a1:a2 (",opt$Output_dir,'/',opt$Output,"-MichPrep-Duplicate.snplist).\n", sep='')
}
cat(N_SNP_lowRsq," variants had an Rsq < 0.8.\n", sep='')
cat('After exclusions ',N_SNP_Rsq3,' variants remain.\n',sep='')
cat("IDs were updated to RSIDs for ",N_SNP_RSID," variants.\n", sep='')
sink()

######
# 4. Merge per chromosome files
######

sink(file = paste0(opt$Output_dir,'/', opt$Output,'-PostMich.log'), append = T)
cat('Merging per chromosome files...',sep='')
sink()

# Create merge.list.
system(paste0('ls ',opt$Output_dir,"/temp/*.rsid.fam | sed 's/.fam//g' > ",opt$Output_dir,'/temp/merge.list'))

# Merge the per chromosome files.
system(paste0(opt$plink,' --merge-list ',opt$Output_dir,'/temp/merge.list --make-bed --out ',opt$Output_dir,'/',opt$Output))

sink(file = paste0(opt$Output_dir,'/', opt$Output,'-PostMich.log'), append = T)
cat('Done!\n',sep='')
sink()

# Delete temporary files.
system(paste0('rm -r ',opt$Output_dir,'/temp'))

# Count final number of variants
N_SNP_final<-system(paste0('wc -l ',opt$Output_dir,'/',opt$Output,".bim"), intern=T)
N_SNP_final<-as.numeric(unlist(strsplit(N_SNP_final,' '))[1])

end.time <- Sys.time()
time.taken <- end.time - start.time

sink(file = paste0(opt$Output_dir,'/', opt$Output,'-PostMich.log'), append = T)
cat('After merging ',N_SNP_final,' variants remained.\n',sep='')
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,3)),attr(time.taken, 'units'),'\n')
sink()








