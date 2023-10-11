library(TwoSampleMR)
TL<-extract_instruments(outcomes = 'ebi-a-GCST90017020',p1=1e-05,clump = TRUE,access_token = NULL)
dim(TL)
outcome_dat<-extract_outcome_data(snps = TL$SNP,outcomes = 'ieu-a-988',proxies=FALSE,)
dat<-harmonise_data(TL,outcome_dat,action = 3)   
head(dat)

snp_add_eaf <- function(dat, build = "37", pop = "EUR")
{
  stopifnot(build %in% c("37","38"))
  stopifnot("SNP" %in% names(dat))
  
  server <- ifelse(build == "37","http://grch37.rest.ensembl.org","http://rest.ensembl.org")
  pop <- paste0("1000GENOMES:phase_3:",pop)
  
  snp_reverse_base <- function(x)
  {
    x <- str_to_upper(x)
    stopifnot(x %in% c("A","T","C","G"))
    switch(x,"A"="T","T"="A","C"="G","G"="C")
  }
  
  res_tab <- lapply(1:nrow(dat), function(i)
  {
    print(paste0("seaching for No.", i, " SNP"))
    dat_i <- dat[i,]
    
    ext <- paste0("/variation/Homo_sapiens/",dat_i$SNP, "?content-type=application/json;pops=1")
    url <- paste(server, ext, sep = "")
    res <- httr::GET(url)
    
    httr::stop_for_status(res)
    
    res <- httr::content(res)
    res_pop <- jsonlite::fromJSON(jsonlite::toJSON(res))$populations
  
    res_pop <- try(res_pop[res_pop$population == pop,])
    if("try-error" %in% class(res_pop))
    {
      print(paste0("There is not information for population ",pop))
      queried_effect_allele <- "NR"
      queried_other_allele <- "NR"
      queried_eaf <- -1
    }
    else
    {
      if(nrow(res_pop)==0)
      {
        print(paste0("There is not information for population ",pop))
        queried_effect_allele <- "NR"
        queried_other_allele <- "NR"
        queried_eaf <- -1
      }
      else
      {
        queried_effect_allele <- res_pop[1,"allele"][[1]]
        queried_other_allele <- res_pop[2,"allele"][[1]]
        queried_eaf <- res_pop[1,"frequency"][[1]]    
      }
    }
    
    effect_allele <- ifelse("effect_allele.exposure" %in% names(dat),
                            dat_i$effect_allele.exposure,
                            dat_i$effect_allele)
    
    other_allele <- ifelse("effect_allele.exposure" %in% names(dat),
                           dat_i$other_allele.exposure,
                           dat_i$other_allele)
    
    if("effect_allele.exposure" %in% names(dat))
    {
      name_output <- unique(c(names(dat), "eaf.exposure","reliability.exposure"))
    }
    else
    {
      name_output <- unique(c(names(dat), "eaf","reliability.exposure"))
    }
    
    len_effect_allele <- nchar(effect_allele)
    len_other_allele <- nchar(other_allele)
    
    if(len_effect_allele==1&len_other_allele==1)
    {
      if((queried_effect_allele==effect_allele & queried_other_allele==other_allele)|
         (queried_effect_allele==other_allele & queried_other_allele==effect_allele))
      {
        dat_i$eaf.exposure <- ifelse(effect_allele == queried_effect_allele,
                                     queried_eaf,
                                     1-queried_eaf)
        dat_i$eaf <- dat_i$eaf.exposure 
        dat_i$reliability.exposure <- "high"
      }
      else
      {
        r_queried_effect_allele <- snp_reverse_base(queried_effect_allele)
        r_queried_other_allele <- snp_reverse_base(queried_other_allele)
        if((r_queried_effect_allele==effect_allele & r_queried_other_allele==other_allele)|
           (r_queried_effect_allele==other_allele & r_queried_other_allele==effect_allele))
        {
          dat_i$eaf.exposure <- ifelse(effect_allele == r_queried_effect_allele,
                                       queried_eaf,
                                       1-queried_eaf)
          dat_i$eaf <- dat_i$eaf.exposure 
          dat_i$reliability.exposure <- "high"
        }
        else
        {
          dat_i$eaf.exposure <- ifelse(effect_allele == queried_effect_allele,
                                       queried_eaf,
                                       1-queried_eaf)
          dat_i$eaf <- dat_i$eaf.exposure 
          dat_i$reliability.exposure <- "low"
        }
      }
    }
    
    else
    {
 
      short_allele <- ifelse(len_effect_allele==1,
                             effect_allele,
                             other_allele)
      short_allele_eaf <- ifelse(short_allele == queried_effect_allele, 
                                 queried_eaf, 
                                 1-queried_eaf)
      dat_i$eaf.exposure <- ifelse(effect_allele == short_allele,
                                   short_allele_eaf,
                                   1-short_allele_eaf)
      dat_i$eaf <- dat_i$eaf.exposure 
      dat_i$reliability.exposure <- "low"
    }
    
    dat_i[name_output]
  })
  
  return(do.call(rbind, res_tab))
}



dat <- snp_add_eaf(dat)

dat$EAF2 <- (1 - dat$eaf.exposure)
dat$MAF <- pmin(dat$eaf.exposure, dat$EAF2)
PVEfx <- function(BETA, MAF, SE, N){
  pve <- (2*(BETA^2)*MAF*(1 - MAF))/((2*(BETA^2)*MAF*(1 - MAF)) + ((SE^2)*2*N*MAF*(1 - MAF)))
  return(pve) 
}
dat$PVE <- mapply(PVEfx, dat$beta.exposure, dat$MAF, dat$se.exposure, N = dat$samplesize.exposure)
dat$FSTAT <- ((dat$samplesize.exposure - 1 - 1)/1)*(dat$PVE/(1 - dat$PVE))  #F值
results<-mr(dat)  


DAT<-read.csv('datF.csv',header = T)
heterogeneity <- mr_heterogeneity(DAT)  
heterogeneity
mr(DAT,method_list=c('mr_ivw_mre')) 
mr
pleio <- mr_pleiotropy_test(DAT)

library(MRPRESSO)  
mr_presso(BetaOutcome = 'beta.outcome',
          BetaExposure = 'beta.exposure', 
          SdOutcome = 'se.outcome', 
          SdExposure = 'se.exposure', 
          data = DAT, OUTLIERtest = TRUE, 
          DISTORTIONtest = TRUE, SignifThreshold = 0.05, NbDistribution = 5000, seed = NULL)
mrpresso_result <- mr_presso(BetaOutcome = 'beta.outcome',
                             BetaExposure = 'beta.exposure', 
                             SdOutcome = 'se.outcome', 
                             SdExposure = 'se.exposure', 
                             data = DAT, 
                             OUTLIERtest = TRUE, 
                             DISTORTIONtest = TRUE, 
                             SignifThreshold = 0.05, 
                             NbDistribution = 5000, 
                             seed = NULL)
library(MRcML)
summary(dat)

cML_result = mr_cML(dat$beta.exposure,
                    dat$beta.outcome,
                    dat$se.exposure,
                    dat$se.outcome,
                    n = 7008,
                    random_start = 100,
                    random_seed = 1,
)
cML_result
cat("cML-MA: beta.hat=", cML_result$MA_BIC_theta, " se=", cML_result$MA_BIC_se, "pval=", cML_result$MA_BIC_p, "\n")
result <- data.frame(beta.hat = cML_result$MA_BIC_theta,
                     se = cML_result$MA_BIC_se,
                     pval = cML_result$MA_BIC_p)
