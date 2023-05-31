######### basic cliamte data ################
# library(imputeTS)
library(lubridate)
# library(ggplot2)

names<-list.files(path = "./Data_Clipick")

for (i in 1:length(names)){
  EC<-read.csv(paste("./Data_Clipick/",names[i],sep = ''))
  EC$date<-as.character(as.Date(paste(EC$Year,EC$Month,EC$Day,sep='.'),format = '%Y.%m.%d'))
  
  EC$PAR<-EC$rss*0.44*4.56
  # MJ / sq meter TO mol/m2 # 4.56 (?mol/J* or ?mol s???1W*???1)
  # 0.44 Percents of the energy of solar radiation is actually in the 400 - 700 nm range  
  EC$TAir<-(EC$tasmax+EC$tasmin)/2
  SVP <-  610.7*10^(7.5*EC$TAir/(237.3+EC$TAir))
  EC$VPD <- SVP*(1-(EC$hursmax+EC$hursmin)/2/100)/1000
  EC$Precip<-EC$pr
  write.csv(EC,file=paste('./Weather_PRELES/',names[i],sep = ''))  
}

######################################################################
####### Global daily CO2  #################

rm(list = ls())
library(imputeTS)
library(lubridate)
library(ggplot2)

names<-list.files(path = "./Weather_PRELES/")

glo<-read.csv('CO2.csv')
glo$Day<-15
glo$date<-as.character(as.Date(paste(glo$Year,glo$Month,glo$Day,sep='.'),format = '%Y.%m.%d'))
glo<-as.data.frame(cbind(glo$date,glo$trend))
colnames(glo)<-c('date','CO2')
glo$CO2<-as.numeric(as.character(glo$CO2))

for (i in 1:length(names)){
  EC<-read.csv(paste("./Weather_PRELES/",names[i],sep = ''))
  EC$date<-as.character(as.Date(paste(EC$Year,EC$Month,EC$Day,sep='.'),format = '%Y.%m.%d'))
  total<-merge(EC,glo,by='date',all.x = T)
  summary(total)
  total$CO2<-na_interpolation(total$CO2,option='linear')
  EC$date<-as.character(EC$date)
  total<-merge(EC,glo,by='date',all.x = T)
  # summary(total)
  # total$CO2.y<-na.seasplit(total$CO2,algorithm = "interpolation")
  
  total$CO2<-na_interpolation(total$CO2,option='linear')
  # total$CO2<-na.kalman(total$CO2,smooth = T)
  
  total$date<-ymd(total$date)
  check<-ggplot(data = total)+
    geom_line(aes(x=date,y=CO2))+
    # geom_line(aes(x=date,y=CO2.x),color='red')+
    theme_bw()
  ggsave(plot=check,paste('./CheckCO2/',names[i],'co2-glo.jpg'),width =15 ,height =15 ,units = 'cm')
  
  data<-data.frame(Date=total$date,
                   Year=total$Year,
                   Month=total$Month,
                   Day=total$Day,
                   PAR = total$PAR,
                   TAir= total$TAir,
                   VPD= total$VPD,
                   Precip= total$Precip,
                   CO2= total$CO2
  )
  summary(data)
  write.csv(data,file=paste('./Weather_PRELES_CO2/',names[i],sep = ''))
}
#########################################################################################

