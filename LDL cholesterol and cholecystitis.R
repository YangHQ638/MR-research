
library(TwoSampleMR)

library(MRPRESSO)
LDL<-extract_instruments(outcomes = 'ieu-a-300',clump = TRUE,
                                   r2=0.01,kb=5000,access_token = NULL)
attach(LDL)
fstatistic.exposure<-beta.exposure^2/se.exposure^2

LDL<-cbind(LDL,fstatistic.exposure)

detach(LDL)
dim(LDL)
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
