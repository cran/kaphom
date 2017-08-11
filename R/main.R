fleisshom<-function(pp, pm, mm)
{
  j <- length(pp)

  tot <- NULL
  prevs <- NULL
  kappas <- NULL

  for(i in 1:j) {
    tot <- append(tot,(pp[i]+pm[i]+mm[i]))
  }

  for(i in 1:j) {
    prev <- (2*pp[i]+pm[i])/(2*tot[i])
    kap <- ((4*pp[i]*mm[i])-(pm[i]*pm[i]))/((2*pp[i]+pm[i])*(2*mm[i]+pm[i]))
    prevs <- append(prevs,prev)
    kappas <- append(kappas,kap)
  }

  ###########Fleiss Homogeneity Test (Different Prevelences)
  vars <- NULL
  weis <- NULL

  for(i in 1:j) {
    variance <- ((1-kappas[i])/tot[i])*(((1-kappas[i])*(1-2*kappas[i]))+((kappas[i]*(2-kappas[i]))/(2*prevs[i]*(1-prevs[i]))))
    vars <- append(vars,variance)
  }

  for(i in 1:j) {
    wei <- (1/vars[i])
    weis <- append(weis,wei)
  }

  sumwk <- 0
  sumw <- 0

  for(i in 1:j) {
    sumwk <- sumwk + (weis[i]*kappas[i])
    sumw <- sumw + weis[i]
  }

  kw <- sumwk/sumw

  diffleiss <- 0

  for(i in 1:j) {
    diffleiss <- diffleiss + weis[i]*((kappas[i]-kw)*(kappas[i]-kw))
  }

  pdiffleiss <- pchisq(diffleiss,(j-1), lower.tail=FALSE)

  if(pdiffleiss > 0.05){
    print("The difference between the given kappa statistics is NOT statistically significant without the equal prevelences assumption!")
    print(paste0("Fleiss Homogeneity Test Statistic : ", diffleiss))
    print(paste0("P-value : ", pdiffleiss))
  }else {
    print("The difference between the given kappa statistics is statistically significant!")
    print(paste0("Fleiss Homogeneity Test Statistic : ", diffleiss))
    print(paste0("P-value : ", pdiffleiss))
  }
  ##################################################################
  ###########Fleiss Homogeneity Test (Equal Prevelences)

  compp <- 0
  compm <- 0
  commm <- 0
  comtot <- 0

  for(i in 1:j) {
    compp <- compp + pp[i]
    compm <- compm + pm[i]
    commm <- commm + mm[i]
    comtot <- comtot + tot[i]
  }

  prevcom <- (2*compp+compm)/(2*comtot)
  kapcom <- ((4*compp*commm)-(compm*compm))/((2*compp+compm)*(2*commm+compm))

  vareq <- (1-kapcom)*(1-kapcom)*((1-(2*kapcom))+((kapcom*(2-kapcom))/((2*prevcom)*(1-prevcom)*(1-kapcom))))

  weieqs <- NULL

  for(i in 1:j) {
    weieq <- (tot[i]/vareq)
    weieqs <- append(weieqs,weieq)
  }

  eqfleiss <- 0

  for(i in 1:j) {
    eqfleiss <- eqfleiss + weieqs[i]*(kappas[i]-kapcom)*(kappas[i]-kapcom)
  }

  peqfleiss <- pchisq(eqfleiss,(j-1), lower.tail=FALSE)

  if(peqfleiss > 0.05){
    print("The difference between the given kappa statistics is NOT statistically significant with the assumption of equal prevelence!")
    print(paste0("Fleiss Homogeneity Test Statistic : ", eqfleiss))
    print(paste0("P-value : ", peqfleiss))
  }else {
    print("The difference between the given kappa statistics is statistically significant!")
    print(paste0("Fleiss Homogeneity Test Statistic : ", eqfleiss))
    print(paste0("P-value : ", peqfleiss))
  }

}

donnerhom<-function(pp, pm, mm)
{
  j <- length(pp)

  tot <- NULL
  prevs <- NULL
  kappas <- NULL

  for(i in 1:j) {
    tot <- append(tot,(pp[i]+pm[i]+mm[i]))
  }

  for(i in 1:j) {
    prev <- (2*pp[i]+pm[i])/(2*tot[i])
    kap <- ((4*pp[i]*mm[i])-(pm[i]*pm[i]))/((2*pp[i]+pm[i])*(2*mm[i]+pm[i]))
    prevs <- append(prevs,prev)
    kappas <- append(kappas,kap)
  }

  comk <- 0
  comkpay <- 0
  comkpayda <- 0
  p2s <- NULL
  p1s <- NULL
  p0s <- NULL
  p2 <- 0
  p1 <- 0
  p0 <- 0

  for(i in 1:j) {
    comkpay <- comkpay + tot[i]*prevs[i]*(1- prevs[i])*kappas[i]
    comkpayda <- comkpayda + tot[i]*prevs[i]*(1- prevs[i])
  }

  comk <- comkpay / comkpayda

  for(i in 1:j) {
    p2 <- (prevs[i]*prevs[i]) + (prevs[i]*(1-prevs[i])*kappas[i])
    p2s <- append(p2s,signif(p2, digits = 2))
    p1 <- (2*prevs[i])*(1 - prevs[i])*(1 - kappas[i])
    p1s <- append(p1s,signif(p1, digits = 2))
    p0 <- ((1 - prevs[i])*(1 - prevs[i])) + (prevs[i]*(1 - prevs[i])*kappas[i])
    p0s <- append(p0s,signif(p0, digits = 2))
  }

  donnergof <- 0

  for(i in 1:j) {
    donnergof <- donnergof + ((mm[i] - (tot[i]*p0s[i]))*(mm[i] - (tot[i]*p0s[i]))/(tot[i]*p0s[i]))
    donnergof <- donnergof + ((pm[i] - (tot[i]*p1s[i]))*(pm[i] - (tot[i]*p1s[i]))/(tot[i]*p1s[i]))
    donnergof <- donnergof + ((pp[i] - (tot[i]*p2s[i]))*(pp[i] - (tot[i]*p2s[i]))/(tot[i]*p2s[i]))
  }

  pdonnergof <- pchisq(donnergof,(j-1), lower.tail=FALSE)

  if(pdonnergof > 0.05){
    print("The difference between the given kappa statistics is NOT statistically significant!")
    print(paste0("Donner GOF Homogeneity Test Statistic : ", donnergof))
    print(paste0("P-value : ", pdonnergof))
  }else {
    print("The difference between the given kappa statistics is statistically significant!")
    print(paste0("Donner GOF Homogeneity Test Statistic : ", donnergof))
    print(paste0("P-value : ", pdonnergof))
  }

}


lscorehom<-function(pp, pm, mm)
{
  j <- length(pp)

  tot <- NULL
  prevs <- NULL
  kappas <- NULL

  for(i in 1:j) {
    tot <- append(tot,(pp[i]+pm[i]+mm[i]))
  }

  for(i in 1:j) {
    prev <- (2*pp[i]+pm[i])/(2*tot[i])
    kap <- ((4*pp[i]*mm[i])-(pm[i]*pm[i]))/((2*pp[i]+pm[i])*(2*mm[i]+pm[i]))
    prevs <- append(prevs,prev)
    kappas <- append(kappas,kap)
  }

  compp <- 0
  compm <- 0
  commm <- 0
  comtot <- 0

  for(i in 1:j) {
    compp <- compp + pp[i]
    compm <- compm + pm[i]
    commm <- commm + mm[i]
    comtot <- comtot + tot[i]
  }

  prevcom <- (2*compp+compm)/(2*comtot)
  kapcom <- ((4*compp*commm)-(compm*compm))/((2*compp+compm)*(2*commm+compm))

  lsskors <- NULL
  lsvars <- NULL
  lszskors <- NULL
  lsskor <- 0
  lsvar <- 0

  for(i in 1:j) {
    lsskor <- ((pp[i]/(prevs[i] + ((1-prevs[i])*kapcom))) + (mm[i]/((1-prevs[i]) + (prevs[i]*kapcom))) - tot[i])/(1-kapcom)
    lsskors <- append(lsskors,lsskor)
    lsvar <- tot[i]/((1-kapcom)*((1-kapcom)*(1-(2*kapcom))+((kapcom*(2-kapcom))/(2*prevs[i]*(1-prevs[i])))))
    lsvars <- append(lsvars,lsvar)
    lszskors <- append(lszskors,(lsskor/sqrt(lsvar)))
  }

  lsdif <- 0

  for(i in 1:j) {
    lsdif <- lsdif + (lszskors[i]*lszskors[i])
  }

  plsdif <- pchisq(lsdif,(j-1), lower.tail=FALSE)

  if(plsdif > 0.05){
    print("The difference between the given kappa statistics is NOT statistically significant without the assumption of equal prevelance!")
    print(paste0("Likelihood Score Homogeneity Test Statistic : ", lsdif))
    print(paste0("P-value : ", plsdif))
  }else {
    print("The difference between the given kappa statistics is statistically significant without the assumption of equal prevelance!")
    print(paste0("Likelihood Score Homogeneity Test Statistic : ", lsdif))
    print(paste0("P-value : ", plsdif))
  }

  ##################################################################
  ###########Likelihood Score Test Equal Prevelances

  eqlsskors <- NULL
  eqlsvars <- NULL
  eqlsskor <- 0
  eqlsvar <- 0

  for(i in 1:j) {
    eqlsskor <- ((pp[i]/(prevcom + ((1-prevcom)*kapcom))) + (mm[i]/((1-prevcom) + (prevcom*kapcom))) - tot[i])/(1-kapcom)
    eqlsskors <- append(eqlsskors,eqlsskor)
    eqlsvar <- tot[i]/((1-kapcom)*((1-kapcom)*(1-(2*kapcom))+((kapcom*(2-kapcom))/(2*prevcom*(1-prevcom)))))
    eqlsvars <- append(eqlsvars,eqlsvar)
  }

  lseq <- 0

  for(i in 1:j) {
    lseq <- lseq + ((eqlsskors[i]*eqlsskors[i])/eqlsvars[i])
  }

  plseq <- pchisq(lseq,(j-1), lower.tail=FALSE)

  if(plseq > 0.05){
    print("The difference between the given kappa statistics is NOT statistically significant with the assumption of equal prevelance!")
    print(paste0("Likelihood Score Homogeneity Test Statistic : ", lseq))
    print(paste0("P-value : ", plseq))
  }else {
    print("The difference between the given kappa statistics is statistically significant with the assumption of equal prevelance!")
    print(paste0("Likelihood Score Homogeneity Test Statistic : ", lseq))
    print(paste0("P-value : ", plseq))
  }

}


mlscorehom<-function(pp, pm, mm)
{

  j <- length(pp)

  tot <- NULL
  prevs <- NULL
  kappas <- NULL

  for(i in 1:j) {
    tot <- append(tot,(pp[i]+pm[i]+mm[i]))
  }

  for(i in 1:j) {
    prev <- (2*pp[i]+pm[i])/(2*tot[i])
    kap <- ((4*pp[i]*mm[i])-(pm[i]*pm[i]))/((2*pp[i]+pm[i])*(2*mm[i]+pm[i]))
    prevs <- append(prevs,prev)
    kappas <- append(kappas,kap)
  }

  comk <- 0
  comkpay <- 0
  comkpayda <- 0

  for(i in 1:j) {
    comkpay <- comkpay + tot[i]*prevs[i]*(1- prevs[i])*kappas[i]
    comkpayda <- comkpayda + tot[i]*prevs[i]*(1- prevs[i])
  }

  comk <- comkpay / comkpayda

  mlsskors <- NULL
  mlsvars <- NULL
  mlsskor <- 0
  mlsvar <- 0
  comk <- signif(comk, digits = 3)

  for(i in 1:j) {
    mlsskor <- ((pp[i]/(prevs[i] + ((1-prevs[i])*comk))) + (mm[i]/((1-prevs[i]) + (prevs[i]*comk))) - tot[i])/(1-comk)
    mlsskors <- append(mlsskors,mlsskor)
    mlsvar <- tot[i]/((1-comk)*((1-comk)*(1-(2*comk))+((comk*(2-comk))/(2*prevs[i]*(1-prevs[i])))))
    mlsvars <- append(mlsvars,mlsvar)
  }

  mls <- 0
  mlstotskor <- 0
  mlstotvar <- 0

  for(i in 1:j) {
    mlstotskor <- mlstotskor + (mlsskors[i]*mlsskors[i])
    mlstotvar <- mlstotvar + mlsvars[i]
  }

  mlssecond <- mlstotskor/mlstotvar

  for(i in 1:j) {
    mls <- mls + ((mlsskors[i]*mlsskors[i])/mlsvars[i])
  }

  mls <- mls - mlssecond

  pmls <- pchisq(mls,(j-1), lower.tail=FALSE)

  if(pmls > 0.05){
    print("The difference between the given kappa statistics is NOT statistically significant!")
    print(paste0("Modified Likelihood Score Homogeneity Test Statistic : ", mls))
    print(paste0("P-value : ", pmls))
  }else {
    print("The difference between the given kappa statistics is statistically significant!")
    print(paste0("Modified Likelihood Score Homogeneity Test Statistic : ", mls))
    print(paste0("P-value : ", pmls))
  }

}


pearsonhom<-function(pp, pm, mm)
{

  j <- length(pp)

  tot <- NULL
  prevs <- NULL
  kappas <- NULL

  for(i in 1:j) {
    tot <- append(tot,(pp[i]+pm[i]+mm[i]))
  }

  for(i in 1:j) {
    prev <- (2*pp[i]+pm[i])/(2*tot[i])
    kap <- ((4*pp[i]*mm[i])-(pm[i]*pm[i]))/((2*pp[i]+pm[i])*(2*mm[i]+pm[i]))
    prevs <- append(prevs,prev)
    kappas <- append(kappas,kap)
  }
  totpp <- 0
  totpm <- 0
  totmm <- 0
  totn <- 0

  for(i in 1:j) {
    totpp <- totpp + pp[i]
    totpm <- totpm + pm[i]
    totmm <- totmm + mm[i]
  }

  totn <- totpp + totpm + totmm

  pcom <- (2*totpp+totpm)/(2*totn)
  kcom <- ((4*totpp*totmm)-(totpm*totpm))/((2*totpp+totpm)*(2*totmm+totpm))
  qcom <- 1-pcom

  pgof1 <- 0
  pgof2 <- 0
  pgof3 <- 0

  for(i in 1:j) {
    pgof1 <- pgof1 + (pp[i]*pp[i])/tot[i]
  }
  pgof1 <- pgof1/(pcom*(pcom+(qcom*kcom)))

  for(i in 1:j) {
    pgof2 <- pgof2 + (pm[i]*pm[i])/tot[i]
  }
  pgof2 <- pgof2/(2*pcom*qcom*(1-kcom))

  for(i in 1:j) {
    pgof3 <- pgof3 + (mm[i]*mm[i])/tot[i]
  }
  pgof3 <- pgof3/(qcom*(qcom+(pcom*kcom)))

  pgof <- pgof1 + pgof2 + pgof3 - totn

  ppgof <- pchisq(pgof,2*(j-1), lower.tail=FALSE)

  if(ppgof > 0.05){
    print("The difference between the given kappa statistics is NOT statistically significant!")
    print(paste0("Pearson GOF Homogeneity Test Statistic : ", pgof))
    print(paste0("P-value : ", ppgof))
  }else {
    print("The difference between the given kappa statistics is statistically significant!")
    print(paste0("Pearson GOF Homogeneity Test Statistic : ", pgof))
    print(paste0("P-value : ", ppgof))
  }

}
