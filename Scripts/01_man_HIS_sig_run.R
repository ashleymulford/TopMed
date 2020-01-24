library(data.table)
library(dplyr)
library(qqman)
library(colorspace)
library(tibble)

#Create function to paste in drug name
"%&%" = function(a,b) paste(a,b,sep="")

#Create color vector
colors<-sequential_hcl(4,"SunsetDark")

HIS_pqtl<-fread("/home/ashley/TopMed/trans_eQTLs_HIS_pqtl_trans_10e-4.txt")
HIS_pqtl_sig<-filter(HIS_pqtl, HIS_pqtl$pvalue<5e-20)

for (protein in HIS_pqtl_sig$gene){
  if(exists("proteins_list")){
    proteins_list<-append(proteins_list, protein)}
  else{
    proteins_list<-c(protein)}
}

protein_set<-unique(proteins_list)

for(protein in protein_set){
  if(exists("HIS_pqtl_man")){
    HIS_pqtl_man<-rbind(HIS_pqtl_man, filter(HIS_pqtl, HIS_pqtl$gene==protein))
  }
  else {
    HIS_pqtl_man<-filter(HIS_pqtl, HIS_pqtl$gene==protein)
  }
}


for(snp in HIS_pqtl_man$snps){
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

HIS_pqtl_man<-add_column(HIS_pqtl_man, chr = chr_list, .before = "gene")
HIS_pqtl_man<-add_column(HIS_pqtl_man, bp = bp_list, .before = "gene")


for(protein in protein_set){
  HIS_pqtl_manh<-filter(HIS_pqtl_man, HIS_pqtl_man$gene==protein)
  png(filename = "/home/ashley/TopMed/Sig_Man_Plots/HIS_pqtl_" %&% protein %&% ".manplot.png", res=100)
  manhattan(HIS_pqtl_manh, chr = "chr", bp = "bp", p = "pvalue", col = colors)
  dev.off()
}



for(protein in protein_set){
  HIS_pqtl_manh<-filter(HIS_pqtl_man, HIS_pqtl_man$gene==protein)
  for (pval in HIS_pqtl_manh$pvalue) {
    if (pval<5e-8){
      if(exists("HIS_pqtl_manha")){
        HIS_pqtl_manha<-rbind(HIS_pqtl_manha, filter(HIS_pqtl_manh, HIS_pqtl_manh$pvalue==pval))
      }
      else {
        HIS_pqtl_manha<-filter(HIS_pqtl_manh, HIS_pqtl_manh$pvalue==pval)
      }
    }
  }
}


HIS_peaks_chrs<-distinct(select(HIS_pqtl_manha, 2,4))
HIS_peaks_chrs_trans<-HIS_peaks_chrs[HIS_peaks_chrs$gene %in% HIS_peaks_chrs$gene[duplicated(HIS_peaks_chrs$gene)],]
fwrite(HIS_peaks_chrs_trans, "/home/ashley/TopMed/trans_proteins_HIS.txt", sep = "\t")



for (protein in HIS_peaks_chrs_trans$gene){
  if(exists("trans_proteins_list")){
    trans_proteins_list<-append(trans_proteins_list, protein)}
  else{
    trans_proteins_list<-c(protein)}
}

trans_protein_set<-unique(trans_proteins_list)


for(protein in trans_protein_set){
  HIS_pqtl_manhatt<-filter(HIS_pqtl_man, HIS_pqtl_man$gene==protein)
  png(filename = "/home/ashley/TopMed/Sig_Trans_Man_Plots/HIS_pqtl_" %&% protein %&% ".manplot.png", res=100)
  manhattan(HIS_pqtl_manhatt, chr = "chr", bp = "bp", p = "pvalue", col = colors)
  dev.off()
}
