
library(TwoSampleMR)

library(MRPRESSO)
Triglycerides<-extract_instruments(outcomes = 'ieu-a-302',clump = TRUE,
                         r2=0.01,kb=5000,access_token = NULL)
attach(Triglycerides)
fstatistic.exposure<-beta.exposure^2/se.exposure^2

Triglycerides<-cbind(Triglycerides,fstatistic.exposure)

dim(Triglycerides)
detach(Triglycerides)


cholecystitis<-extract_outcome_data(snps = Triglycerides$SNP,
                                    outcomes = 'ukb-d-K81',
                                    proxies = TRUE,
                                    maf_threshold = 0.01,
                                    access_token = NULL)
dim(cholecystitis)

mydata<-harmonise_data(exposure_dat = Triglycerides,
                       outcome_dat = cholecystitis,
                       action = 2)


mydata<-subset(mydata, SNP!='rs7005265')
mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = mydata, NbDistribution = 1000,  
          SignifThreshold = 0.05)

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


het<-mr_heterogeneity(mydata,method_list = c("mr_ivw_mre",
                                             "mr_egger_regression",
                                             "mr_ivw_fe",
                                             "mr_weighted_median"))


pleio<-mr_pleiotropy_test(mydata)

single<-mr_leaveoneout(mydata)
mr_leaveoneout_plot(single)
