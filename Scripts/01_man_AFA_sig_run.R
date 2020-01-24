#Import libraries
library(data.table)
library(dplyr)
library(qqman)
library(colorspace)
library(tibble)

#Create function to paste in drug name
"%&%" = function(a,b) paste(a,b,sep="")

#Create color vector
colors<-sequential_hcl(4,"SunsetDark")

#Read in file and filter for only significant p-vals (threshol chosen: 5e-20)
AFA_pqtl<-fread("/home/ashley/TopMed/trans_eQTLs_AFA_pqtl_trans_10e-4.txt")
AFA_pqtl_sig<-filter(AFA_pqtl, AFA_pqtl$pvalue<5e-20)

#Generate list of proteins with significnat p-vals
for (protein in AFA_pqtl_sig$gene){
  if(exists("proteins_list")){
    proteins_list<-append(proteins_list, protein)}
  else{
    proteins_list<-c(protein)}
}

#Set for unique list
protein_set<-unique(proteins_list)

#Using set, subset original file to get all snps for proteins with sig hits
for(protein in protein_set){
  if(exists("AFA_pqtl_man")){
    AFA_pqtl_man<-rbind(AFA_pqtl_man, filter(AFA_pqtl, AFA_pqtl$gene==protein))
  }
  else {
    AFA_pqtl_man<-filter(AFA_pqtl, AFA_pqtl$gene==protein)
  }
}

#Split snp col into two cols, chr and bp
for(snp in AFA_pqtl_man$snps){
  snpID<-strsplit(snp, ":")
  chrom<-snpID[[1]][1]
  bp<-snpID[[1]][2]
  bp<-as.numeric(bp)
  chro<-strsplit(chrom, "r")
  chr<-chro[[1]][2]
  chr<-as.numeric(chr)
  if(exists("chr_list")){
    chr_list<-append(chr_list, chr)}
  else{
    chr_list<-c(chr)}
  if(exists("bp_list")){
    bp_list<-append(bp_list, bp)}
  else{
    bp_list<-c(bp)}
}  

#Add cols
AFA_pqtl_man<-add_column(AFA_pqtl_man, chr = chr_list, .before = "gene")
AFA_pqtl_man<-add_column(AFA_pqtl_man, bp = bp_list, .before = "gene")

#Make manhattan plots for each sig protein
for(protein in protein_set){
  AFA_pqtl_manh<-filter(AFA_pqtl_man, AFA_pqtl_man$gene==protein)
  png(filename = "/home/ashley/TopMed/Sig_Man_Plots/AFA_pqtl_" %&% protein %&% ".manplot.png", res=100)
  manhattan(AFA_pqtl_manh, chr = "chr", bp = "bp", p = "pvalue", col = colors)
  dev.off()
}



#Identidfy possible trans proteins:
#Create new dataframe with all snps with pval<5e-8
for(protein in protein_set){
  AFA_pqtl_manh<-filter(AFA_pqtl_man, AFA_pqtl_man$gene==protein)
  for (pval in AFA_pqtl_manh$pvalue) {
    if (pval<5e-8){
      if(exists("AFA_pqtl_manha")){
        AFA_pqtl_manha<-rbind(AFA_pqtl_manha, filter(AFA_pqtl_manh, AFA_pqtl_manh$pvalue==pval))
      }
      else {
        AFA_pqtl_manha<-filter(AFA_pqtl_manh, AFA_pqtl_manh$pvalue==pval)
      }
    }
  }
}

#Create new dataframe with only chr and protein id, keep only unique rows (ie chr/proteins combos)
AFA_peaks_chrs<-distinct(select(AFA_pqtl_manha, 2,4))
#Keep only duplicate proteins (ie proteins will more than one chr with hit at pval<5e-8, possibly trans)
AFA_peaks_chrs_trans<-AFA_peaks_chrs[AFA_peaks_chrs$gene %in% AFA_peaks_chrs$gene[duplicated(AFA_peaks_chrs$gene)],]
fwrite(AFA_peaks_chrs_trans, "/home/ashley/TopMed/trans_proteins_AFA.txt", sep = "\t")


#Get list of possible trans proteins
for (protein in AFA_peaks_chrs_trans$gene){
  if(exists("trans_proteins_list")){
    trans_proteins_list<-append(trans_proteins_list, protein)}
  else{
    trans_proteins_list<-c(protein)}
}

#Make set
trans_protein_set<-unique(trans_proteins_list)

#Make manhattan plots for each possibly trans protein
for(protein in trans_protein_set){
  AFA_pqtl_manhatt<-filter(AFA_pqtl_man, AFA_pqtl_man$gene==protein)
  png(filename = "/home/ashley/TopMed/Sig_Trans_Man_Plots/AFA_pqtl_" %&% protein %&% ".manplot.png", res=100)
  manhattan(AFA_pqtl_manhatt, chr = "chr", bp = "bp", p = "pvalue", col = colors)
  dev.off()
}

