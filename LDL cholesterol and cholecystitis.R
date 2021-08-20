
library(TwoSampleMR)

library(MRPRESSO)
LDL<-extract_instruments(outcomes = 'ieu-a-300',clump = TRUE,
                                   r2=0.01,kb=5000,access_token = NULL)
attach(LDL)
fstatistic.exposure<-beta.exposure^2/se.exposure^2
r_squared<-(2*(1-eaf.exposure)*eaf.exposure*beta.exposure)/(se.exposure*173082^0.5)
LDL<-cbind(LDL,fstatistic.exposure,r_squared)
sum(r_squared,na.rm = TRUE)
detach(LDL)
dim(LDL)
#write.csv(LDL,file = "D:/mr/t2d and cholecystitis/新建文件夹/LDL exopusre.csv")
cholecystitis<-extract_outcome_data(snps = LDL$SNP,
                                    outcomes = 'ukb-d-K81',
                                    proxies = TRUE,
                                    maf_threshold = 0.01,
                                    access_token = NULL)
dim(cholecystitis)

mydata<-harmonise_data(exposure_dat = LDL,
                       outcome_dat = cholecystitis,
                       action = 2)
mydata<-subset(mydata, SNP!='rs2954029')
mydata<-subset(mydata, SNP!='rs519113')
mydata<-subset(mydata, SNP!='rs964184')

mydata<-subset(mydata, SNP!='rs11591147')
mydata<-subset(mydata, SNP!='rs1801689')
mydata<-subset(mydata, SNP!='rs13277801')
mydata<-subset(mydata, SNP!='rs10195252')
mydata<-subset(mydata, SNP!='rs6544713')
mydata<-subset(mydata, SNP!='rs4970834')
mydata<-subset(mydata, SNP!='rs1510226')
mydata<-subset(mydata, SNP!='rs6016381')
mydata<-subset(mydata, SNP!='rs10903129')
mydata<-subset(mydata, SNP!='rs1475701')
mydata<-subset(mydata, SNP!='rs688')
mydata<-subset(mydata, SNP!='rs2710642')
mydata<-subset(mydata, SNP!='rs8017377')
mydata<-subset(mydata, SNP!='rs17800760')
mydata<-subset(mydata, SNP!='rs7640978')
mydata<-subset(mydata, SNP!='rs2287019')
mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = mydata, NbDistribution = 5000,  
          SignifThreshold = 0.05)
#mydata<-mydata[-c(2,67),]
sum(mydata$r_squared,na.rm = TRUE)
write.csv(mydata[,c('SNP','effect_allele.exposure','other_allele.exposure','effect_allele.outcome','other_allele.outcome','beta.exposure','beta.outcome','se.outcome',"se.exposure",'pval.outcome','pval.exposure','fstatistic.exposure')],file = "D:/mr/t2d and cholecystitis/新建文件夹/LDL identified SNPs.csv")
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
          file = "D:/mr/t2d and cholecystitis/新建文件夹/result of LDL.csv")

het<-mr_heterogeneity(mydata,method_list = c("mr_ivw_mre",
                                             "mr_egger_regression",
                                             "mr_ivw_fe",
                                             "mr_weighted_median"))
write.csv(het[,c('method','Q_pval','Q')],
          file = "D:/mr/t2d and cholecystitis/新建文件夹/heterogeneity test of LDL.csv")

pleio<-mr_pleiotropy_test(mydata)
write.csv(pleio[,c('egger_intercept','pval')],
          file = "D:/mr/t2d and cholecystitis/新建文件夹/pleiotropy test of LDL.csv")
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
                                                ggplot2::labs(y = "", x = "MR leave-one-out sensitivity analysis for\n LDL cholesterol on cholecystitis")
                                            })
  res
}
loo(single)

mr_scatter_plot(result,mydata)
fscatter<-function (mr_results, dat) 
{
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("plyr", quietly = TRUE)
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), 
                       function(d) {
                         d <- plyr::mutate(d)
                         if (nrow(d) < 2 | sum(d$mr_keep) == 0) {
                           return(blank_plot("Insufficient number of SNPs"))
                         }
                         d <- subset(d, mr_keep)
                         index <- d$beta.exposure < 0
                         d$beta.exposure[index] <- d$beta.exposure[index] * 
                           -1
                         d$beta.outcome[index] <- d$beta.outcome[index] * 
                           -1
                         mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & 
                                           id.outcome == d$id.outcome[1])
                         mrres$a <- 0
                         if ("MR Egger" %in% mrres$method) {
                           temp <- mr_egger_regression(d$beta.exposure, 
                                                       d$beta.outcome, d$se.exposure, d$se.outcome, 
                                                       default_parameters())
                           mrres$a[mrres$method == "MR Egger"] <- temp$b_i
                         }
                         if ("MR Egger (bootstrap)" %in% mrres$method) {
                           temp <- mr_egger_regression_bootstrap(d$beta.exposure, 
                                                                 d$beta.outcome, d$se.exposure, d$se.outcome, 
                                                                 default_parameters())
                           mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
                         }
                         ggplot2::ggplot(data = d, ggplot2::aes(x = beta.exposure, 
                                                                y = beta.outcome)) + ggplot2::geom_errorbar(ggplot2::aes(ymin = beta.outcome - 
                                                                                                                           se.outcome, ymax = beta.outcome + se.outcome), 
                                                                                                            colour = "grey", width = 0) + ggplot2::geom_errorbarh(ggplot2::aes(xmin = beta.exposure - 
                                                                                                                                                                                 se.exposure, xmax = beta.exposure + se.exposure), 
                                                                                                                                                                  colour = "grey", height = 0) + ggplot2::geom_point(ggplot2::aes(text = paste("SNP:", 
                                                                                                                                                                                                                                               SNP))) + ggplot2::geom_abline(data = mrres, ggplot2::aes(intercept = a, 
                                                                                                                                                                                                                                                                                                        slope = b, colour = method), show.legend = TRUE) + 
                           ggplot2::scale_colour_manual(values = c("#a6cee3", 
                                                                   "#1f78b4", "#b2df8a", "#33a02c", 
                                                                   "#fb9a99", "#e31a1c", "#fdbf6f", 
                                                                   "#ff7f00", "#cab2d6", "#6a3d9a", 
                                                                   "#ffff99", "#b15928")) + ggplot2::labs(colour = "MR Test", 
                                                                                                          x = 'SNP effect on LDL cholesterol', 
                                                                                                          y = 'SNP effect on cholecystitis') + 
                           ggplot2::theme(legend.position = "top", 
                                          legend.direction = "vertical") + ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2))
                       })
  mrres
}
fscatter(result,mydata)
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
                                              ggplot2::labs(y = "", x = "MR effect size for\n LDL cholesterol on cholecystitis")
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
