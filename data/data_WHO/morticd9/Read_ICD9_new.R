ICD9 <- read.table("Morticd9",sep=",",header=TRUE)
# 5020,Australia;4010,Austria;4020,Belgium;2090,Canada;4050,Denmark;4070,Finland;4080,France;4140,Greece;4150,Hungary;4160,Iceland;4170,Ireland;4180,Italy;3160,Japan;4190,Luxembourg;4210,Netherlands;4220,Norway;4230,Poland;4240,Portugal;5150,New Zealand;4280,Spain;4290,Sweden;4300,Switzerland;4308,United Kingdom;2450,United States of America

mycountry <- c(5020,4010,4020,2090,4050,4070,4080,4140,4150,4160,4170,4180,3160,4190,4210,4220,4230,4240,5150,4280,4290,4300,4308,2450)

ICD9 <- ICD9[which(ICD9$Country %in% mycountry),]

colnames(ICD9)[10:35] <- c("All","Age0","Age1","Age2","Age3","Age4","Age5-9","Age10-14","Age15-19","Age20-24","Age25-29","Age30-34","Age35-39","Age40-44","Age45-49","Age50-54","Age55-59","Age60-64","Age65-69","Age70-74","Age75-79","Age80-84","Age85-89","Age90-94","Age95+","AgeUNK")

#PROBLEM with HIV (B184): In ICD9 it is classified among Endocrine and nutritional diseases, while in ICD10 is listed among Infective and parasitic diseases, it would be better to hamonise, so we move deaths for cause B184 among Macro-Cause 1 (Infective and parasitic) also for ICD9. HIV related deaths are not recorded in Japan (3160), Finland(4070), Poland(4230)

InfectList <- c("B01", "B02", "B03", "B04", "B05", "B06", "B07","B184")
NeopList <- c("B08","B09", "B100", "B109", "B11", "B12","B13", "B14", "B15", "B16", "B17")
LungList <- c("B101")
EndocList <- c("B18","B19","B20")
MentalList <- c("B21","B465")
NervList <- c("B22")

##B20 Diseases of blood and blood-forming organs -> Endoc
##B21 Mental disorders -> Mental
##B22 Diseases of the nervous system -> Nervous
##B23 Disorders of the eye and adnexa -> Other
##B24 Diseases of the ear and mastoid process -> Other
EndocListB <- c("B180","B181","B182","B183","B185","B189","B19")
CircList <- c("B25", "B26", "B27", "B28", "B29","B30")
RespList <- c("B31", "B32")
DigestList <- c("B33", "B34")
##B35 Diseases of urinary system -> Uro
##B36-B37 Diseases of male genital organs - Diseases of female genital organs -> Uro
##B38 Abortion -> Infant
##B39 Direct obstetric causes -> Infant
##B40 Indirect obstetric causes, B41 Normal Delivery (???) -> Infant
##B42 Diseases of skin and subcutaneous tissue -> Skin
##B43 Diseases of the musculoskeletal system and connective tissue -> Skin
##B44 Congenital anomalies -> Congenital
##B45 Certain conditions originating in the perinatal period - Infant
##
UroList <- c("B35","B36","B37")
SkinList <- c("B42","B43")
CongeList <- c("B44")
InfantList <- c("B38","B39","B40","B41","B45")
ExternList <- c("B47", "B48", "B49","B50","B51","B52","B53","B54","B55","B56","CH17")
OthList <- c("B23","B24","B461","B462","B463","B464","B466","B467","B469")

ICD9$CauseM <- NA
ICD9$CauseM <- ifelse(ICD9$Cause%in%InfectList, 1,ICD9$CauseM)
ICD9$CauseM <- ifelse(ICD9$Cause%in%NeopList, 2,ICD9$CauseM)
ICD9$CauseM <- ifelse(ICD9$Cause%in%LungList, 3,ICD9$CauseM)
ICD9$CauseM <- ifelse(((ICD9$Cause%in%EndocList)&(ICD9$Country%in%c(3160,4070,4230))), 4,ICD9$CauseM)
ICD9$CauseM <- ifelse(((ICD9$Cause%in%EndocListB)&(ICD9$Country%in%mycountry[-c(6,13,17)])), 4,ICD9$CauseM)
ICD9$CauseM <- ifelse(ICD9$Cause%in%CircList, 5,ICD9$CauseM)
ICD9$CauseM <- ifelse(ICD9$Cause%in%RespList, 6,ICD9$CauseM)
ICD9$CauseM <- ifelse(ICD9$Cause%in%DigestList, 7,ICD9$CauseM)
ICD9$CauseM <- ifelse(ICD9$Cause%in%ExternList, 8,ICD9$CauseM)
ICD9$CauseM <- ifelse(ICD9$Cause%in%MentalList, 9,ICD9$CauseM)
ICD9$CauseM <- ifelse(ICD9$Cause%in%NervList, 10,ICD9$CauseM)
ICD9$CauseM <- ifelse(ICD9$Cause%in%UroList, 11,ICD9$CauseM)
ICD9$CauseM <- ifelse(ICD9$Cause%in%SkinList, 12,ICD9$CauseM)
ICD9$CauseM <- ifelse(ICD9$Cause%in%CongeList, 13,ICD9$CauseM)
ICD9$CauseM <- ifelse(ICD9$Cause%in%InfantList, 14,ICD9$CauseM)
ICD9$CauseM <- ifelse(ICD9$Cause%in%OthList, 15,ICD9$CauseM)


AllDeaths <- ICD9[ICD9$Cause=="B00",c(1,4,7,10:35)]
SpecDeaths <- ICD9[!is.na(ICD9$CauseM),]
AggS.Deaths <- aggregate(SpecDeaths[,c(10:35)], by=SpecDeaths[,c("Country","Year","Sex","CauseM")], sum) 



AllDeaths2 <- reshape(AllDeaths,varying = c("All","Age0","Age1","Age2","Age3","Age4","Age5-9","Age10-14", "Age15-19","Age20-24", "Age25-29", "Age30-34", "Age35-39", "Age40-44", "Age45-49", "Age50-54", "Age55-59","Age60-64", "Age65-69", "Age70-74", "Age75-79", "Age80-84", "Age85-89", "Age90-94", "Age95+","AgeUNK"),timevar="Age", times=c("All","0","1","2","3","4","5-9","10-14", "15-19","20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59","60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "90-94", "95+","UNK"),direction="long",v.names = "Counts")
AllDeaths2$CauseM <- 0

AggS.Deaths2 <- reshape(AggS.Deaths,varying = c("All","Age0","Age1","Age2","Age3","Age4","Age5-9","Age10-14", "Age15-19","Age20-24", "Age25-29", "Age30-34", "Age35-39", "Age40-44", "Age45-49", "Age50-54", "Age55-59","Age60-64", "Age65-69", "Age70-74", "Age75-79", "Age80-84", "Age85-89", "Age90-94", "Age95+","AgeUNK"),timevar="Age", times=c("All","0","1","2","3","4","5-9","10-14", "15-19","20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59","60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "90-94", "95+","UNK"),direction="long",v.names = "Counts")
final <- rbind(AggS.Deaths2,AllDeaths2[,c(1,2,3,7,4,5,6)])

complete <- NULL
for (j in 1:length(table(final$Country))){
    c <- as.numeric(names(table(final$Country)))[j]
    for(i in 1:length(table(final$Year[final$Country==c]))){
        y <- as.numeric(names(table(final$Year[final$Country==c])))[i]
        All0 <- final[final$Country==c&final$Year==y&final$CauseM==0,]
        final2 <- final[final$CauseM>0&final$Country==c&final$Year==y,]
        Agg2 <- aggregate(final2[,6], by=final2[,c("Age","Sex")],sum, na.rm=TRUE)
        rest <- merge(Agg2,All0,by=c("Age","Sex"))
        rest$Counts2 <- rest$Counts-rest$x
        rest$CauseM <- 16
        rest <- rest[,c(4,5,2,6,1,9,8)]
        colnames(rest)[6] <- "Counts"
        final3 <- rbind(final2,rest)
        total <- All0[,-c(4,7)]
        colnames(total)[5] <- "AllCauses"
        comp <- merge(final3[,-7],total, by=c("Country","Year","Sex","Age"))
        complete <- rbind(complete,comp)
    }
}

write.csv(complete,file="ICD9_comp_new.csv",sep=",")        
