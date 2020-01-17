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

#Read file
HIS_pqtl<-fread("/home/ashley/TopMed/trans_eQTLs_HIS_pqtl_trans_10e-4.txt")

#Split snp col into two cols, chr and bp
for(snp in HIS_pqtl$snps){
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
HIS_pqtl<-add_column(HIS_pqtl, chr = chr_list, .before = "gene")
HIS_pqtl<-add_column(HIS_pqtl, bp = bp_list, .before = "gene")

#Make list of proteins
for (protein in HIS_pqtl$gene){
  if(exists("proteins_list")){
    proteins_list<-append(proteins_list, protein)}
  else{
    proteins_list<-c(protein)}
}

#Unique proteins only
protein_set<-unique(proteins_list)

#Make manhattan plot for each protein
for(protein in protein_set){
  HIS_pqtl_man<-filter(HIS_pqtl, HIS_pqtl$gene==protein)
  png(filename = "HIS_pqtl_" %&% protein %&% ".manplot.png", res=100)
  manhattan(HIS_pqtl_man, chr = "chr", bp = "bp", p = "pvalue", col = colors)
  dev.off()
}

