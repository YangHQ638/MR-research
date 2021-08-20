
library(TwoSampleMR)

library(MRPRESSO)
cholesterol<-extract_instruments(outcomes = 'ieu-a-301',clump = TRUE,
                                   r2=0.01,kb=5000,access_token = NULL)
attach(cholesterol)
fstatistic.exposure<-beta.exposure^2/se.exposure^2
r_squared<-(2*(1-eaf.exposure)*eaf.exposure*beta.exposure)/(se.exposure*187365^0.5)
cholesterol<-cbind(cholesterol,fstatistic.exposure,r_squared)
sum(r_squared,na.rm = TRUE)
detach(cholesterol)

dim(cholesterol)
#write.csv(cholesterol,file = "D:/mr/t2d and cholecystitis/新建文件夹/cholesterol exopusre.csv")
cholecystitis<-extract_outcome_data(snps = cholesterol$SNP,
                                    outcomes = 'ukb-d-K81',
                                    proxies = TRUE,
                                    maf_threshold = 0.01,
                                    access_token = NULL)
dim(cholecystitis)

mydata<-harmonise_data(exposure_dat = cholesterol,
                       outcome_dat = cholecystitis,
                       action = 2)
mydata<-subset(mydata, SNP!='rs2954029')
mydata<-subset(mydata, SNP!='rs964184')

mydata<-subset(mydata, SNP!='rs7412')
mydata<-subset(mydata, SNP!='rs10468017')
mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = mydata, NbDistribution = 5000,  
          SignifThreshold = 0.05)
mydata<-mydata[-c(96,88,99),]

sum(mydata$r_squared,na.rm = TRUE)
write.csv(mydata[,c('SNP','effect_allele.exposure','other_allele.exposure','effect_allele.outcome','other_allele.outcome','beta.exposure','beta.outcome','se.outcome',"se.exposure",'pval.outcome','pval.exposure','fstatistic.exposure')],
          file = "D:/mr/t2d and cholecystitis/新建文件夹/cholesterol identified SNPs.csv")
result<-mr(mydata,method_list = c("mr_ivw_mre",
                                  "mr_egger_regression",
                                  "mr_ivw_fe",
                                  "mr_weighted_median"
))
#计算OR值和OR的可信区间
attach(result)
newresult<-transform(result,
                     OR=exp(b),
                     upper=exp(b+1.96*se*b),
                     lower=exp(b-1.96*se*b))
detach(result)

write.csv(newresult[,c('method','pval','OR','upper','lower')],
          file = "D:/mr/t2d and cholecystitis/新建文件夹/result of cholesterol.csv")

het<-mr_heterogeneity(mydata,method_list = c("mr_ivw_mre",
                                             "mr_egger_regression",
                                             "mr_ivw_fe",
                                             "mr_weighted_median"))
write.csv(het[,c('method','Q_pval','Q')],
          file = "D:/mr/t2d and cholecystitis/新建文件夹/heterogeneity test of cholesterol.csv")

pleio<-mr_pleiotropy_test(mydata)
write.csv(pleio[,c('egger_intercept','pval')],
          file = "D:/mr/t2d and cholecystitis/新建文件夹/pleiotropy test of cholesterol.csv")

single<-mr_leaveoneout(mydata)
mr_leaveoneout_plot(single)
loo<-function (leaveoneout_results) 
{
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("plyr", quietly = TRUE)
  res <- plyr::dlply(leaveoneout_results, c("id.exposure", 
                                            "id.outcome"), function(d) {
                                              d <- plyr::mutate(d)
                                              if (sum(!grepl("All", d$SNP)) < 3) {
                                                return(blank_plot("Insufficient number of SNPs"))
                                              }
                                              d$up <- d$b + 1.96 * d$se
                                              d$lo <- d$b - 1.96 * d$se
                                              d$tot <- 1
                                              d$tot[d$SNP != "All"] <- 0.01
                                              d$SNP <- as.character(d$SNP)
                                              nom <- d$SNP[d$SNP != "All"]
                                              nom <- nom[order(d$b)]
                                              d <- rbind(d, d[nrow(d), ])
                                              d$SNP[nrow(d) - 1] <- ""
                                              d$b[nrow(d) - 1] <- NA
                                              d$up[nrow(d) - 1] <- NA
                                              d$lo[nrow(d) - 1] <- NA
                                              d$SNP <- ordered(d$SNP, levels = c("All", "", 
                                                                                 nom))
                                              ggplot2::ggplot(d, ggplot2::aes(y = SNP, x = b)) + ggplot2::geom_vline(xintercept = 0, 
                                                                                                                     linetype = "dotted") + ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, 
                                                                                                                                                                                 xmax = up, size = as.factor(tot), colour = as.factor(tot)), 
                                                                                                                                                                    height = 0) + ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot))) + 
                                                ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% 
                                                                                                      "")), colour = "grey") + ggplot2::scale_colour_manual(values = c("black", 
                                                                                                                                                                       "red")) + ggplot2::scale_size_manual(values = c(0.3, 
                                                                                                                                                                                                                       1)) + ggplot2::theme(legend.position = "none", 
                                                                                                                                                                                                                                            axis.text.y = ggplot2::element_text(size = 8), axis.ticks.y = ggplot2::element_line(size = 0), 
                                                                                                                                                                                                                                            axis.title.x = ggplot2::element_text(size = 8)) + 
                                                ggplot2::labs(y = "", x = "MR leave-one-out sensitivity analysis for\n total cholesterol on cholecystitis")
                                            })
  res
}
loo(single)

mr_scatter_plot(result,mydata)
res_single<-mr_singlesnp(mydata)
mr_forest_plot(res_single)
fp<-function (singlesnp_results, exponentiate = FALSE) 
{
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("plyr", quietly = TRUE)
  res <- plyr::dlply(singlesnp_results, c("id.exposure", 
                                          "id.outcome"), function(d) {
                                            d <- plyr::mutate(d)
                                            if (sum(!grepl("All", d$SNP)) < 2) {
                                              return(blank_plot("Insufficient number of SNPs"))
                                            }
                                            levels(d$SNP)[levels(d$SNP) == "All - Inverse variance weighted"] <- "All - IVW"
                                            levels(d$SNP)[levels(d$SNP) == "All - MR Egger"] <- "All - Egger"
                                            am <- grep("All", d$SNP, value = TRUE)
                                            d$up <- d$b + 1.96 * d$se
                                            d$lo <- d$b - 1.96 * d$se
                                            d$tot <- 0.01
                                            d$tot[d$SNP %in% am] <- 1
                                            d$SNP <- as.character(d$SNP)
                                            nom <- d$SNP[!d$SNP %in% am]
                                            nom <- nom[order(d$b)]
                                            d <- rbind(d, d[nrow(d), ])
                                            d$SNP[nrow(d) - 1] <- ""
                                            d$b[nrow(d) - 1] <- NA
                                            d$up[nrow(d) - 1] <- NA
                                            d$lo[nrow(d) - 1] <- NA
                                            d$SNP <- ordered(d$SNP, levels = c(am, "", nom))
                                            xint <- 0
                                            if (exponentiate) {
                                              d$b <- exp(d$b)
                                              d$up <- exp(d$up)
                                              d$lo <- exp(d$lo)
                                              xint <- 1
                                            }
                                            ggplot2::ggplot(d, ggplot2::aes(y = SNP, x = b)) + ggplot2::geom_vline(xintercept = xint, 
                                                                                                                   linetype = "dotted") + ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, 
                                                                                                                                                                               xmax = up, size = as.factor(tot), colour = as.factor(tot)), 
                                                                                                                                                                  height = 0) + ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot))) + 
                                              ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% 
                                                                                                    "")), colour = "grey") + ggplot2::scale_colour_manual(values = c("black", 
                                                                                                                                                                     "red")) + ggplot2::scale_size_manual(values = c(0.3, 
                                                                                                                                                                                                                     1)) + ggplot2::theme(legend.position = "none", 
                                                                                                                                                                                                                                          axis.text.y = ggplot2::element_text(size = 8), axis.ticks.y = ggplot2::element_line(size = 0), 
                                                                                                                                                                                                                                          axis.title.x = ggplot2::element_text(size = 8)) + 
                                              ggplot2::labs(y = "", x = "MR effect size for\n total cholesterol on cholecystitis")
                                          })
  res
}
fp(res_single)
mr_funnel_plot(res_single)

library(dplyr)
variants<-mydata$SNP
library(ieugwasr)
phenotype<-phewas(variants)
phenotype<-phenotype%>%filter(p<5e-8)
by_trait<-phenotype%>%group_by(trait)
trait_num<-by_trait%>%tally(sort = TRUE)
by_snp<-phenotype%>%group_by(rsid)
snp_group<-by_snp%>%tally()
pleio_snp<-phenotype%>%filter(trait!='Triglycerides')
pleio_snp_bysnp<-pleio_snp%>%group_by(rsid)
pleio_snp_unoverlapped<-pleio_snp_bysnp%>%tally()


