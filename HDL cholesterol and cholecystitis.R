
library(TwoSampleMR)

library(MRPRESSO)
HDL<-extract_instruments(outcomes = 'ieu-a-299',clump = TRUE,
                                   r2=0.01,kb=5000,access_token = NULL)
attach(HDL)
fstatistic.exposure<-beta.exposure^2/se.exposure^2

HDL<-cbind(HDL,fstatistic.exposure,r_squared)

detach(HDL)
dim(HDL)



cholecystitis<-extract_outcome_data(snps = HDL$SNP,
                                    outcomes = 'ukb-d-K81',
                                    proxies = TRUE,
                                    maf_threshold = 0.01,
                                    access_token = NULL)
dim(cholecystitis)

mydata<-harmonise_data(exposure_dat = HDL,
                       outcome_dat = cholecystitis,
                       action = 2)
mydata<-subset(mydata, SNP!='rs6589581')
mydata<-subset(mydata, SNP!='rs964184')

mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = mydata, NbDistribution = 1000,  
          SignifThreshold = 0.05)



result<-mr(mydata,method_list = c("mr_ivw_mre",
                                  "mr_egger_regression",
                                  "mr_ivw_fe",
                                  "mr_weighted_median"
                                  ))
attach(result)
newresult<-transform(result,
                     OR=exp(b))
detach(result)


het<-mr_heterogeneity(mydata,method_list = c("mr_ivw_mre",
                                             "mr_egger_regression",
                                             "mr_ivw_fe",
                                             "mr_weighted_median"))


pleio<-mr_pleiotropy_test(mydata)


single<-mr_leaveoneout(mydata)
mr_leaveoneout_plot(single)
