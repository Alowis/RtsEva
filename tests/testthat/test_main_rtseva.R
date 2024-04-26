##########################################################################################
###################   DEMO of TSEVA for different locations   ############################
##########################################################################################
##########################################################################################


#load input data
data=ArdecheStMartin
catch="Ardeche"
#Station_data_IDs <- as.vector(t(Q_sim[1, ]))[-1]

#tqx=as.POSIXct((Q_sim[,1]*3600*24), origin="1979-01-01 00:00:00")
#txx=tqx[-1]
#######################  Arguments selection #############################
haz = "drought"
tail="low"
var = "dis"
season="allyear"
######################################################################################


start_time <- Sys.time()
timeStamps=data$time
dt = difftime(timeStamps[2],timeStamps[1],units="days")
dt= as.numeric(dt)
percentile=95
names(data)=c("date","Qs")

if (haz=="drought"){
  #seasonality divide: frost vs non frost
  if (!exists("trans")){trans="rev"}
  print(paste0(trans," transformation used for low flows"))
  #compute 7days moving average
  data$Q7=tsEvaNanRunningMean(data$Qs,7/dt)
  timeAndSeries=data.frame(data$date,data$Q7)

}else if (haz=="flood"){
  percentile=95
  timeAndSeries <- max_daily_value(timeAndSeries)

}

#new check for timestamps
names(timeAndSeries)=c("timestamp","dis")
dt1=min(diff(timeAndSeries$timestamp),na.rm=T)
dt=as.numeric(dt1)
tdim=attributes(dt1)$units
if (tdim=="hours") dt=dt/24
if (dt==1){
  timeDays=timeAndSeries$timestamp
}else{
  timeDays=unique(as.Date(timeAndSeries$timestamp))
}

bounds=c(year(timeDays[1]),year(timeDays[length(timeDays)]))
tbound=c(as.Date(paste0(bounds[1],"-12-31")),as.Date(paste0(bounds[2],"-12-31")))
Impdates=seq(tbound[1],tbound[2],by="1 year")
tindexes=match(Impdates,timeDays)
nv=length(unique(timeAndSeries$dis))

names(timeAndSeries)=c("timestamp","data")
tsm=1/dt
series=timeAndSeries[,2]
timeWindow = 365.25*30; #time windows in days, the correction is done within the functions
windowSize=366
minPeakDistanceInDays=30
timeStamps=timeAndSeries$timestamp
transftypes=c("rev","inv","lninv")
trendtypes=c("trend","trendPeaks")
trendtrans=expand.grid(transftypes,trendtypes)
#choose transformation
tt=4
plot(timeAndSeries)
Nonstat<-TsEvaNs(timeAndSeries, timeWindow, transfType=trendtrans[tt,2], ciPercentile= 90, minPeakDistanceInDays = minPeakDistanceInDays, tail=tail,TrendTh=0.85, trans=trendtrans[tt,1])
nonStationaryEvaParams=Nonstat[[1]]
stationaryTransformData=Nonstat[[2]]

ExRange= c(min(nonStationaryEvaParams$potObj$parameters$peaks),max(nonStationaryEvaParams$potObj$parameters$peaks))

if (haz=="flood") wr2 <- c(seq(min(ExRange),max(ExRange),length.out=700))
if (haz=="drought") wr2 <- c(seq(1.1*min(ExRange),0.1*max(ExRange),length.out=700))

Plot1= tsEvaPlotGPDImageScFromAnalysisObj(wr2, nonStationaryEvaParams, stationaryTransformData, minYear = '1950',trans=trendtrans[tt,1])
timeIndex=1
Plot2 = tsEvaPlotReturnLevelsGPDFromAnalysisObj(nonStationaryEvaParams, stationaryTransformData, timeIndex, trans=trendtrans[tt,1],ylabel="Discharge (m3/s)")
Plot3 = tsEvaPlotReturnLevelsGEVFromAnalysisObj(nonStationaryEvaParams, stationaryTransformData, timeIndex, trans=trendtrans[tt,1],ylabel="Discharge (m3/s)")
stationaryTransformData$timeStampsDay=unique(as.Date(stationaryTransformData$timeStamps))

if (dt==1){
  timeDays=stationaryTransformData$timeStamps
}else{
  timeDays=stationaryTransformData$timeStampsDay
}

bounds=c(year(timeDays[1]),year(timeDays[length(timeDays)]))
tbound=c(as.Date(paste0(bounds[1],"-12-31")),as.Date(paste0(bounds[2],"-12-31")))
Impdates=seq(tbound[1],tbound[2],by="1 years")

#For drought case, 7 days timestamps
datex=yday(timeDays)
dtect=c(diff(datex),-1)
last_days <- timeDays[which(dtect<0)]
tindexes=match(last_days,timeDays)

# Compute return periods and levels
RPgoal=100
timeIndex=tindexes[1]
RLevs100=tsEvaComputeRLsGEVGPD(nonStationaryEvaParams, RPgoal, timeIndex, trans=trendtrans[tt,1])

if (RLevs100$Fit == "No fit") {

  RLgpd <- nonStationaryEvaParams$gevObj$parameters$annualMax
  names(RLgpd) <- year(Impdates)

  RLgev <- nonStationaryEvaParams$gevObj$parameters$annualMax
  names(RLgev) <- year(Impdates)

  nRPgev <- rep(NA, length(Impdates))
  names(nRPgev) <- year(Impdates)

  nRPgpd <- rep(NA, length(Impdates))
  names(nRPgpd) <- year(Impdates)

  params=data.frame(matrix(ncol=17,nrow=(length(Impdates)-1)))
  params[,1]=rep(catch,7)
  params[,2]=year(Impdates)[-1]
  params[,4]=rep(interflag,7)
  if (is.null(colnames(parlist))){
    colnames(params)=rep("nom",17)
  }else{
    colnames(params)=colnames(parlist)
  }
}else{
  #####
  RLgev=RLevs100$ReturnLevels[2]
  RLgpd=RLevs100$ReturnLevels[3]
  ERgev=RLevs100$ReturnLevels[4]
  ERgpd=RLevs100$ReturnLevels[5]
  nRPgev=nRPgpd=100
  params=c()
  for (t in 2:length(Impdates)){
    timeIndex=tindexes[t]
    RLevs100i=tsEvaComputeRLsGEVGPD(nonStationaryEvaParams, RPgoal, timeIndex, trans=trendtrans[tt,1])
    params=c(catch,year(Impdates[t]),timeIndex,RLevs100i$Params)
    names(params)[1:3]=c("catchment","Year","timeIndex")
    Rper=tsEvaComputeTimeRP(params=RLevs100i$Params,RPiGEV=RLevs100$ReturnLevels[2],RPiGPD=RLevs100$ReturnLevels[3])
    nRPgpd=c(nRPgpd,Rper[2])
    nRPgev=c(nRPgev,Rper[1])
    RLgev=cbind(RLgev,RLevs100i$ReturnLevels[2])
    RLgpd=cbind(RLgpd,RLevs100i$ReturnLevels[3])
    ERgev=cbind(ERgev,RLevs100i$ReturnLevels[4])
    ERgpd=cbind(ERgpd,RLevs100i$ReturnLevels[5])
  }

  RLgev=as.data.frame(RLgev)
  names(RLgev)=year(Impdates)
  rownames(RLgev)=RPgoal

  RLgpd=as.data.frame(RLgpd)
  names(RLgpd)=year(Impdates)
  rownames(RLgpd)=RPgoal

  nRPgev=as.data.frame(t(nRPgev))
  names(nRPgev)=year(Impdates)

  nRPgpd=as.data.frame(t(nRPgpd))
  names(nRPgpd)=year(Impdates)
}
