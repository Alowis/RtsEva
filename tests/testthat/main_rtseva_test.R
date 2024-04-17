##########################################################################################
###################   DEMO of TSEVA for different locations   ############################
##########################################################################################
#Import functions from TSEVA
source("~/LFRuns_utils/TSEVA_demo/demo_functions.R")

##########################################################################################

main_path = 'D:/tilloal/Documents/06_Floodrivers/'
valid_path = paste0(main_path,'HERA/')
dis_path<-paste0(main_path,'dis/calibrated/filtered/Histo/')
setwd(valid_path)


#load input data
load(file="HERA6h_H07_19502020.Rdata")
Station_data_IDs <- as.vector(t(Q_sim[1, ]))[-1]

tqx=as.POSIXct((Q_sim[,1]*3600*24), origin="1979-01-01 00:00:00")
txx=tqx[-1]
#######################  Arguments importation #############################
haz = "drought"
var = "dis"
outlets="Hybas07"
outletname = "outletsv8_hybas07_01min"
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")

######################################################################################
library(xts)

#load frost days
load(file=paste0(hydroDir,"/Drought/catchment_frost.Rdata"))
#Hybas07
Catchmentrivers7=read.csv(paste0(hydroDir,"/Catchments/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)
outhyb07=outletopen(hydroDir,outletname)
catmatch=match(outhyb07$outlets,Catchmentrivers7$pointid)
mycat=Catchmentrivers7[catmatch,]

hybas07 <- read_sf(dsn = paste0(hydroDir,"/Catchments/hydrosheds/hybas_eu_lev07_v1c.shp"))
hybasf7=fortify(hybas07)
Catamere07=inner_join(hybasf7,mycat,by= "HYBAS_ID")
Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ")
Catf7=Catamere07

st_geometry(Catf7)=NULL

#remove locations without discharge

Q_sim2=Q_sim[,-which(is.na(Q_sim[2,]))]
#choose the timesere to analyse
id=1905
outlet=Q_sim2[1,id]
df.dis=Q_sim2[,id][-1]

#Use this until i finish to generate csv of hybas outlets
#remove 1950
rmv=which(year(txx)==1950)
df.dis=data.frame(txx,df.dis)
names(df.dis)[c(1,2)]=c("time","dis")
#remove 1950 which is not reliable
df.dis=df.dis[-rmv,]

haz="drought"
season="nonfrost"
start_time <- Sys.time()
print(paste0(id,"/",length(Station_data_IDs)))
catch=Station_data_IDs[id]
timeStamps=txx
dt = difftime(timeStamps[2],timeStamps[1],units="days")
dt= as.numeric(dt)

if (haz=="drought"){
  #seasonality divide: frost vs non frost
  percentile=95
  catmat=Catf7[which(Catf7$pointid==catch),]
  Tcatmat=mycat[which(mycat$HYBAS_ID==catmat$HYBAS_ID),]
  Tcatchment=which(colnames(frostcat)==Tcatmat$pointid)
  data=df.dis
  names(data)=c("date","Qs")
  intermit=interid(data, WindowSize=7)
  interflag=intermit$flags[2]
  if (!exists("trans")){trans="rev"}
  print(paste0(trans," transformation used for low flows"))
  series=data.frame(txx[-rmv],intermit$trdis$Q7)
  #remove frost timesteps, this can be modified to do the anlysis only on frost moments
  if (length(Tcatchment)>0){
    frostserie=data.frame(frostcat[,1],frostcat[,Tcatchment])
    names(frostserie)=c("time","Ta")
    frosttime=which(frostserie[,2]<0)
  }else{
    frosttime=NA
  }
}else if (haz=="flood"){
  percentile=95
  df.disX=disNcopenloc(filename,hydroDir,outhybas,idfix)
  series=data.frame(txx,df.disX$outlets)
  interflag=0
}

names(series)=c("timestamp","dis")
dt1=min(diff(series$timestamp),na.rm=T)
dt=as.numeric(dt1)
tdim=attributes(dt1)$units
if (tdim=="hours") dt=dt/24
if (dt==1){
  timeDays=series$timestamp
}else{
  timeDays=unique(as.Date(series$timestamp))
}

bounds=c(year(timeDays[1]),year(timeDays[length(timeDays)]))
realbound=10*ceiling(bounds/10)
tbound=c(as.Date(paste0(realbound[1],"-12-31")),as.Date(paste0(realbound[2],"-12-31")))
Impdates=seq(tbound[1],tbound[2],by="1 year")
tindexes=match(Impdates,timeDays)
names(series)=c("timestamp","dis")
nv=length(unique(series$dis))

timeAndSeries=series
names(timeAndSeries)=c("timestamp","data")

if (haz=="drought" && length(!is.na(frosttime))>1) {
  if (season=="nonfrost") {
    print("nonfrost season")
    timeAndSeries$data[frosttime]=NA
  } else if (season=="frost") {
    print("frost season")
    timeAndSeries$data[-frosttime]=NA
  } else {
    print("season must be frost or nonfrost")
  }
}

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
Nonstat<-TsEvaNs(timeAndSeries, timeWindow, transfType=trendtrans[tt,2], ciPercentile= 90, minPeakDistanceInDays = minPeakDistanceInDays, mode=haz,TrendTh=0.85, trans=trendtrans[tt,1])
nonStationaryEvaParams=Nonstat[[1]]

stationaryTransformData=Nonstat[[2]]

ExRange= c(min(nonStationaryEvaParams$potObj$parameters$peaks),max(nonStationaryEvaParams$potObj$parameters$peaks))

if (haz=="flood") wr2 <- c(seq(min(ExRange),max(ExRange),length.out=700))
if (haz=="drought") wr2 <- c(seq(1.1*min(ExRange),0.1*max(ExRange),length.out=700))

Plot1= tsEvaPlotGPDImageScFromAnalysisObj(wr2, nonStationaryEvaParams, stationaryTransformData, minYear = '1950',trans="rev")
timeIndex=1
Plot2 = tsEvaPlotReturnLevelsGPDFromAnalysisObj(nonStationaryEvaParams, stationaryTransformData, timeIndex, trans="rev",ylabel="Discharge (m3/s)") # nolint
Plot3 = tsEvaPlotReturnLevelsGEVFromAnalysisObj(nonStationaryEvaParams, stationaryTransformData, timeIndex, trans="rev",ylabel="Discharge (m3/s)")

stationaryTransformData$timeStampsDay=unique(as.Date(stationaryTransformData$timeStamps))
pikos=data.frame(nonStationaryEvaParams$potObj$parameters$peaks,nonStationaryEvaParams$potObj$parameters$peakID,nonStationaryEvaParams$potObj$parameters$peakST, nonStationaryEvaParams$potObj$parameters$peakEN)
names(pikos)=c("value","timeID","tIDstart","tIDend")
pikos$time=timeStamps[pikos$timeID]
pikos$catch=rep(catch,length(pikos[,1]))

plot(timeAndSeries$timestamp,stationaryTransformData$nonStatSeries)
# points(timeAndSeries[pikos$timeID,], col="red")
points(pikos$time,pikos$value, col="blue", pch=16)
#Here I need to convert the timeStamp to a daily one if dt is not 1
dt1=min(diff(timeStamps),na.rm=T)
dt=as.numeric(dt1)
tdim=attributes(dt1)$units
if (tdim=="hours") dt=dt/24
if (dt==1){
  timeDays=stationaryTransformData$timeStamps
}else{
  timeDays=stationaryTransformData$timeStampsDay
}

bounds=c(year(timeDays[1]),year(timeDays[length(timeDays)]))
realbound=10*ceiling(bounds/10)
tbound=c(as.Date(paste0(realbound[1],"-12-31")),as.Date(paste0(realbound[2],"-12-31")))
Impdates=seq(tbound[1],tbound[2],by="1 years")
tindexes=match(Impdates,timeDays)

# Compute return periods and levels
RPgoal=100
timeIndex=tindexes[1]
RLevs100=ComputeReturnLevels(nonStationaryEvaParams, RPgoal, timeIndex)

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
    RLevs100i=ComputeReturnLevels(nonStationaryEvaParams, RPgoal, timeIndex)
    params=c(catch,year(Impdates[t]),timeIndex,RLevs100i$Params)
    names(params)[1:3]=c("catchment","Year","timeIndex")

    Rper=RPcalc(params,RPiGEV=RLevs100$ReturnLevels[2],RPiGPD=RLevs100$ReturnLevels[3])
    nRPgpd=c(nRPgpd,Rper[2])
    nRPgev=c(nRPgev,Rper[1])
    RLgev=cbind(RLgev,RLevs100i$ReturnLevels[2])
    RLgpd=cbind(RLgpd,RLevs100i$ReturnLevels[3])
    ERgev=cbind(ERgev,RLevs100i$ReturnLevels[4])
    ERgpd=cbind(ERgpd,RLevs100i$ReturnLevels[5])
    if (length(parlist)>1) colnames(parlist)=names(params)
    parlist=rbind(parlist,params)
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

  peaklist=rbind(peaklist,pikos)
}
