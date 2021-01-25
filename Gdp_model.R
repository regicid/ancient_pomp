library (panelPomp) 
library (tidyr)
library(foreach)
library(doParallel)
registerDoParallel(48)
getDoParWorkers()

Countries = c("Arab","China","Europe","India")
Results = read.csv("./Results.csv")
Results$gdp = scale(Results$gdp - min(Results$gdp,na.rm = TRUE),center = FALSE)
Results = Results[-1]
#Results$gdp[is.na(Results$gdp)] = 0
#Results$Nobs = log(Results$Nobs+1) 

PARAM = c("a","c","z","sigma","sigma_obs","N_0")
rwsd = rw.sd(a=.1,z = .2,sigma=.1,sigma_obs = .1,N_0=ivp(.1),c=.1)
Csnippet("double eps = rnorm(0,pow(sigma,2));
         N = pow(z,2)*N + a*gdp  + c + eps;
         ") -> evol_diff

Csnippet("
         N = N_0;
         ") -> rinit
Csnippet("
         Nobs = rnorm(N,pow(sigma_obs,2));
	") -> rmeas


Csnippet("
         lik = dnorm(Nobs,N,pow(sigma_obs,2),give_log);
         ") -> dmeas
Pomps = list()

for(country in Countries){
  Data = dplyr::filter(Results,countries==country)
  z = min(which(!is.na(Data$gdp)))
  Gdp = Data[z:nrow(Data),c(2,10)]
  
  #Gdp$diff = rowSums(R/untable(Distances[country,],12)**2)
  Data = Data[z:nrow(Data),c(2,19)]
  Gdp$cum = vector("numeric",length = nrow(Gdp))
  Gdp$cum[2:nrow(Gdp)] = cumsum(Data$Nobs)[1:nrow(Gdp)-1]
  
  Data %>%
    pomp(
      times="date", t0=Data$date[1],
      rinit=rinit,
      covar = covariate_table(Gdp, times= "date"),
      rprocess=discrete_time(evol_diff,delta.t = 100),
      dmeasure = dmeas,
      obsnames = c("Nobs"),
      statenames=c("N"),
      rmeasure=rmeas,
      paramnames=PARAM,
      covarnames = c("gdp","cum")
    ) -> Pomps[[country]]
}

p = rep(0,length(PARAM))
names(p) = PARAM

Model_diff = panelPomp(Pomps,shared = p)
start = c(a = .5,sigma = .2,N_0 = 0,sigma_obs = .3,z = .5,c=0)
Model_diff %>%
  panelPomp::mif2(
    shared.start=unlist(start),
    specific.start = Model_diff@specific,
    Np=2000,
    Nmif=50,
    cooling.fraction.50=0.5,
    cooling.type="hyperbolic",
    rw.sd= rwsd,
    pars = PARAM
  ) -> mf1


#mf1 = mifs2[47][[1]]
#Data = data.frame(mf1@pconv.rec)
#Data$iteration = as.numeric(row.names(mf1@pconv.rec))
#Data = gather(Data,"key","value",-iteration)
#ggplot(Data,aes(iteration,value)) + geom_line()+ facet_wrap(vars(key),scales = "free")

#for(i in 1:length(mifs)){
#  Data = data.frame(mifs[[2]]@pconv.rec)
#  Data$iteration = as.numeric(row.names(mifs[[1]]@pconv.rec))
#  Data = gather(Data,"key","value",-iteration)
#}

lower = c(a = 0,sigma=1,N_0 = 0,sigma_obs=0.7,z = .8,d = -0.6,b = -0.5,c =-.5,e=-0.3,f=-0.3)
upper = c(a = 0.6,sigma=1.4, N_0 = .3,sigma_obs=1,z = 1,d = 1,b = 0.5,c = .5,e=-0.3,f=-0.3)
sobolDesign(lower = lower[PARAM], upper = upper[PARAM], nseq = 48) -> guesses

foreach (guess=iter(guesses,"row"),
         .combine=c, .packages=c("panelPomp"),
         .errorhandling="remove", .inorder=FALSE) %dopar% {
           
           mf1 %>% panelPomp::mif2(shared.start=unlist(guess),
                                   specific.start = Model_diff@specific,
                                   Np=2000,Nmif=4000,cooling.fraction.50=0.3,

                                   cooling.type="hyperbolic",rw.sd= rwsd,pars = PARAM)
           
           
         } -> mifs4


saveRDS(mifs4,"./mifs_pomp")
