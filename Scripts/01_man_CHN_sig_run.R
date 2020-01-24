library(data.table)
library(dplyr)
library(qqman)
library(colorspace)
library(tibble)

#Create function to paste in drug name
"%&%" = function(a,b) paste(a,b,sep="")

#Create color vector
colors<-sequential_hcl(4,"SunsetDark")

CHN_pqtl<-fread("/home/ashley/TopMed/trans_eQTLs_CHN_pqtl_trans_10e-4.txt")
CHN_pqtl_sig<-filter(CHN_pqtl, CHN_pqtl$pvalue<5e-12)

for (protein in CHN_pqtl_sig$gene){
  if(exists("proteins_list")){
    proteins_list<-append(proteins_list, protein)}
  else{
    proteins_list<-c(protein)}
}

protein_set<-unique(proteins_list)

for(protein in protein_set){
  if(exists("CHN_pqtl_man")){
    CHN_pqtl_man<-rbind(CHN_pqtl_man, filter(CHN_pqtl, CHN_pqtl$gene==protein))
  }
  else {
    CHN_pqtl_man<-filter(CHN_pqtl, CHN_pqtl$gene==protein)
  }
}


for(snp in CHN_pqtl_man$snps){
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

CHN_pqtl_man<-add_column(CHN_pqtl_man, chr = chr_list, .before = "gene")
CHN_pqtl_man<-add_column(CHN_pqtl_man, bp = bp_list, .before = "gene")


for(protein in protein_set){
  CHN_pqtl_manh<-filter(CHN_pqtl_man, CHN_pqtl_man$gene==protein)
  png(filename = "/home/ashley/TopMed/Sig_Man_Plots/CHN_pqtl_" %&% protein %&% ".manplot.png", res=100)
  manhattan(CHN_pqtl_manh, chr = "chr", bp = "bp", p = "pvalue", col = colors)
  dev.off()
}


for(protein in protein_set){
  CHN_pqtl_manh<-filter(CHN_pqtl_man, CHN_pqtl_man$gene==protein)
  for (pval in CHN_pqtl_manh$pvalue) {
    if (pval<5e-8){
      if(exists("CHN_pqtl_manha")){
        CHN_pqtl_manha<-rbind(CHN_pqtl_manha, filter(CHN_pqtl_manh, CHN_pqtl_manh$pvalue==pval))
      }
      else {
        CHN_pqtl_manha<-filter(CHN_pqtl_manh, CHN_pqtl_manh$pvalue==pval)
      }
    }
  }
}


CHN_peaks_chrs<-distinct(select(CHN_pqtl_manha, 2,4))
CHN_peaks_chrs_trans<-CHN_peaks_chrs[CHN_peaks_chrs$gene %in% CHN_peaks_chrs$gene[duplicated(CHN_peaks_chrs$gene)],]
fwrite(CHN_peaks_chrs_trans, "/home/ashley/TopMed/trans_proteins_CHN.txt", sep = "\t")



for (protein in CHN_peaks_chrs_trans$gene){
  if(exists("trans_proteins_list")){
    trans_proteins_list<-append(trans_proteins_list, protein)}
  else{
    trans_proteins_list<-c(protein)}
}

trans_protein_set<-unique(trans_proteins_list)

for(protein in trans_protein_set){
  CHN_pqtl_manhatt<-filter(CHN_pqtl_man, CHN_pqtl_man$gene==protein)
  png(filename = "/home/ashley/TopMed/Sig_Trans_Man_Plots/CHN_pqtl_" %&% protein %&% ".manplot.png", res=100)
  manhattan(CHN_pqtl_manhatt, chr = "chr", bp = "bp", p = "pvalue", col = colors)
  dev.off()
}
