library(TwoSampleMR)
library(data.table)
library(tidyverse)
library(dplyr)
FileNames <-list.files(paste0(getwd()),pattern=".gz")
exp_dat_ids <- FileNames
exps <- FileNames
out<-fread("finngen_R7_I9_PHLETHROMBDVTLOW.gz",header = T)
out$trait <- 'DVT'  
outcomeid <- out
rm(out)
head(outcomeid)
outcomeid$N <- 274098
outcome<-format_data(outcomeid,type="outcome",
                     snp_col = "rsids",
                     phenotype_col = "trait",
                     effect_allele_col = "alt",
                     other_allele_col = "ref",
                     beta_col = "beta",
                     se_col = "sebeta",
                     samplesize_col = "N",
                     pval_col = "pval",
                     eaf_col = "af_alt",
                     chr_col = "chrom",
                     pos_col = "pos"
)
rm(outcomeid)


for (qaq in 2:length(exp_dat_ids)) { # 
  exp_dat_id <- exp_dat_ids[qaq]
  exp <- exps[qaq]
  

  exposure<- try(fread(paste0(getwd(),"/",FileNames[qaq]),fill=TRUE),silent = T)
  
  exposure$Phenotype <- FileNames[qaq]
  
  head(exposure)
  

  exposure<-format_data(exposure,type="exposure",
                        snp_col = "MarkerName",
                        phenotype_col = "Phenotype",
                        effect_allele_col = "Allele1",
                        other_allele_col = "Allele2",
                        beta_col = "Effect",
                        se_col = "StdErr",
                        samplesize_col = "TotalSampleSize",
                        pval_col = "`P-value`",
                        eaf_col = "Freq1",
                        chr_col = "Chr",
                        pos_col = "Pos")
  
  
  head(exposure)
  head(outcome)
  
  
  exposure$id.exposure<-FileNames[qaq]
  exposure$exposure<-FileNames[qaq]
  outcome$id.outcome<-outcome$outcome
  
  
 
  ld="D:1000G" 
  wld="D:1000G"
  
  
  LDSC_rg<-function(expo,outcome,an,sample_prev=NA,
                    population_prev=NA,ld,wld,chr_filter=c(1:22),n_blocks=200){
    id.o<-outcome$id.outcome[1]
    id.e<-expo$id.exposure[1]
    
    expo<-expo%>%mutate(Z=beta.exposure/se.exposure)
    expo<-expo%>%select(SNP=SNP,N=samplesize.exposure,Z=Z
                        ,A1=effect_allele.exposure
                        ,A2=other_allele.exposure)
    expo<-as_tibble(expo)
    
    outcome<-outcome%>%mutate(Z=beta.outcome/se.outcome)
    outcome<-outcome%>%select(SNP=SNP,N=samplesize.outcome,Z=Z
                              ,A1=effect_allele.outcome
                              ,A2=other_allele.outcome)
    outcome<-as_tibble(outcome)
    
    
    dat<-list(expo,outcome)
    names(dat)<-c(id.e,id.o)
    
    rm(expo,outcome)
    
    
    res<-try(ldscr::ldsc_rg(dat,ancestry = an,sample_prev=sample_prev,
                            population_prev=population_prev,ld=ld,wld=wld,
                            n_blocks=n_blocks,chr_filter=chr_filter))
    
    return(res)
    
  }
  
  LDSC_res<-LDSC_rg(exposure,outcome,ld=ld,wld=wld)
  
  h2 <- LDSC_res[["h2"]]
  rg <- LDSC_res[["rg"]]
  
  
  
  write.csv(h2,file = paste0("mendelian/",exp,"_h2.csv"), row.names = FALSE)
  write.csv(rg,file = paste0("mendelian/",exp,"_rg.csv"), row.names = FALSE)
  
}