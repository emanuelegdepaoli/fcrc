EA.Contries <- c("BLR","EST","LTU","LVA","RUS","UKR")

for(i in 1:6){
    if (i==1){
        dati <- read.csv(paste(EA.Contries[i],"_d_interm_idr.csv",sep=""),header=TRUE,sep=",",
                         na.strings = c("", ".", "NA"))}
    if (i>1){
        dati <- rbind(dati,read.csv(paste(EA.Contries[i],"_d_interm_idr.csv",sep=""),header=TRUE,sep=",",
                                    na.strings = c("", ".", "NA")))}
}

colnames(dati)[7:32] <- c("All","Age0","Age1-4","Age5-9","Age10-14","Age15-19","Age20-24","Age25-29","Age30-34","Age35-39","Age40-44","Age45-49","Age50-54","Age55-59","Age60-64","Age65-69","Age70-74","Age75-79","Age80-84","Age85+","Age85-89","Age90+","Age90-94","Age95+","Age95-99","Age100+")


InfectList <- c(1:9)

#From C00 to D48 (excluding C33-C34)
NeopList <- c(10:17,20:34)
LungList <- c(18,19)
EndocList <- c(35,36,37,38)
MentalList <- c(39:42)
NervList <- c(43:47)
CircList <- c(48:63)
RespList <- c(64:72)
DigestList <- c(73:81)
SkinList <- c(82,83)
UroList <- c(84:87)
InfantList <- c(88,89,91)
CongeList <- c(90)
ExternList <- c(92:103)

dati$CauseM <- NA
dati$CauseM <- ifelse(dati$cause%in%InfectList, 1,dati$CauseM)
dati$CauseM <- ifelse(dati$cause%in%NeopList, 2,dati$CauseM)
dati$CauseM <- ifelse(dati$cause%in%LungList, 3,dati$CauseM)
dati$CauseM <- ifelse(dati$cause%in%EndocList, 4,dati$CauseM)
dati$CauseM <- ifelse(dati$cause%in%CircList, 5,dati$CauseM)
dati$CauseM <- ifelse(dati$cause%in%RespList, 6,dati$CauseM)
dati$CauseM <- ifelse(dati$cause%in%DigestList, 7,dati$CauseM)
dati$CauseM <- ifelse(dati$cause%in%ExternList, 8,dati$CauseM)
dati$CauseM <- ifelse(dati$cause%in%MentalList, 9,dati$CauseM)
dati$CauseM <- ifelse(dati$cause%in%NervList, 10,dati$CauseM)
dati$CauseM <- ifelse(dati$cause%in%UroList, 11,dati$CauseM)
dati$CauseM <- ifelse(dati$cause%in%SkinList, 12,dati$CauseM)
dati$CauseM <- ifelse(dati$cause%in%CongeList, 13,dati$CauseM)
dati$CauseM <- ifelse(dati$cause%in%InfantList, 14,dati$CauseM)
dati$CauseM <- ifelse(dati$cause==0, 0,dati$CauseM)

AllDeaths <- dati[dati$CauseM==0,c(1,2,3,6,8:26)]
SpecDeaths <- dati[dati$CauseM>0,]
AggS.Deaths <- aggregate(apply(SpecDeaths[,c(8:26)],2,as.numeric), by=SpecDeaths[,c("country","year","sex","CauseM")], sum, na.rm = T) 

AllDeaths2 <- reshape(AllDeaths,varying = c("Age0","Age1-4","Age5-9","Age10-14", "Age15-19","Age20-24", "Age25-29", "Age30-34", "Age35-39", "Age40-44", "Age45-49", "Age50-54", "Age55-59","Age60-64", "Age65-69", "Age70-74", "Age75-79", "Age80-84", "Age85+"),timevar="Age", times=c("0","1-4","5-9","10-14", "15-19","20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59","60-64", "65-69", "70-74", "75-79", "80-84", "85+"),direction="long",v.names = "Counts")
AllDeaths2$CauseM <- 0


AggS.Deaths2 <- reshape(AggS.Deaths,varying = c("Age0","Age1-4","Age5-9","Age10-14", "Age15-19","Age20-24", "Age25-29", "Age30-34", "Age35-39", "Age40-44", "Age45-49", "Age50-54", "Age55-59","Age60-64", "Age65-69", "Age70-74", "Age75-79", "Age80-84", "Age85+"),timevar="Age", times=c("0","1-4","5-9","10-14", "15-19","20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59","60-64", "65-69", "70-74", "75-79", "80-84", "85+"),direction="long",v.names = "Counts")


final <- rbind(AggS.Deaths2[,c(1,2,3,4,5,6)],AllDeaths2[,c(1,2,3,8,5,6)])
final$Counts <- as.numeric(final$Counts)

complete <- NULL
for (j in 1:length(table(final$country))){
    c <- names(table(final$country))[j]
    for(i in 1:length(table(final$year[final$country==c]))){
        y <- as.numeric(names(table(final$year[final$country==c])))[i]
        All0 <- final[final$country==c&final$year==y&final$CauseM==0,]
        final2 <- final[final$CauseM>0&final$country==c&final$year==y,]
        Agg2 <- aggregate(as.numeric(final2[,6]), by=final2[,c("Age","sex")],sum,na.rm=TRUE)
        rest <- merge(Agg2,All0,by=c("Age","sex"))
        rest$Counts2 <- rest$Counts-rest$x
        rest$CauseM <- 16
        rest <- rest[,c(4,5,2,6,1,8)]
        colnames(rest)[6] <- "Counts"
        final3 <- rbind(final2,rest)
        total <- All0[,-c(4,7)]
        colnames(total)[5] <- "AllCauses"
        comp <- merge(final3[,-7],total, by=c("country","year","sex","Age"))
        complete <- rbind(complete,comp)
    }
}

write.csv(complete,file="Eastern_HCoD2.csv",sep=",")        



