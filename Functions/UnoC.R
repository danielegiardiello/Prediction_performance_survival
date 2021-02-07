# Uno C
# Function to calculate Uno's C-index
Uno_C<-function(tfup,status,
                linear.predictor,pred.surv,
                thorizon) {
  
  # Create the data frame
  dt<-data.frame(time=tfup,
                 status=status,
                 lp=linear.predictor)
  
  
  # Save in the data the predicted survival and
  # the predicted mortality at time t (e.g. 5 years)
  dt$predsurv<-pred.surv
  dt$predmort<-1-pred.surv
  names(dt)<-c("time","status","lp","predsurv","predmort")
  
  # Censoring probability
  sfit<-dt %>% survfit(Surv(time,status==0)~1,data=.)
  # Survival probability
  sfit1<-dt %>% survfit(Surv(time,status==1)~1,data=.)
  # Save the inverse censoring weighting
  sfit_df<-data.frame(time=c(0,sfit$time),surv=c(1,sfit$surv), 
                      weight=c(1,sfit$surv),
                      weight2=c(1,sfit$surv^2))
  # probc<-min(sfit_df$surv[sfit_df$time<thorizon])  # save censoring prob at time horizon
  # probt<-min(sfit1$surv[sfit1$time<thorizon]) # save surv prob at time horizon
  
  # Add weights in the data with time, status and linear predictor
  dt2<-sqldf('select a.*, b.weight, b.weight2
            from dt as a
            left join sfit_df as b
            on a.time=b.time')
  dt2$weight[dt2$time==5]<-NA  # Assing weight=NA
  dt2$weight2[dt2$time==5]<-NA  # Assing weight=NA
  
  # Order by time
  dt3<-dt2 %>% arrange(time)
  ch<-0
  dh<-0
  N<-nrow(dt3)
  for (i in 1:N) {
    for (j in 1:N){
      if(dt3$time[[i]]<dt3$time[[j]] &
         dt3$lp[[i]]>dt3$lp[[j]] &
         dt3$time[[i]]<thorizon) {
        ch<-ch+(dt3$status[[i]])/(dt3$weight2[[i]])}
      
      if(dt3$time[[i]]<dt3$time[[j]] & 
         dt3$lp[[i]]==dt3$lp[[j]] &
         dt3$time[[i]]<thorizon) {
        ch<-ch+0.5*(dt3$status[[i]])/(dt3$weight2[[i]])}
      
      if(dt3$time[[i]]<dt3$time[[j]] &
         dt3$time[[i]]<4.95) {
        dh<-dh+(dt3$status[[i]])/(dt3$weight2[[i]])
        
      }
    }
  }
  
  c_hat<-ch/dh
  return(c_hat)
}
