ICD10 <- read.table("Morticd10_part1",sep=",",header=TRUE)
ICD10 <- rbind(ICD10, read.table("Morticd10_part2",sep=",",header=TRUE))
ICD10 <- rbind(ICD10, read.table("Morticd10_part3",sep=",",header=TRUE))
ICD10 <- rbind(ICD10, read.table("Morticd10_part4",sep=",",header=TRUE))
ICD10 <- rbind(ICD10, read.table("Morticd10_part5",sep=",",header=TRUE))
# 5020,Australia;4010,Austria;4020,Belgium;2090,Canada;4050,Denmark;4070,Finland;
# 4080,France;4140,Greece;4150,Hungary;4160,Iceland;4170,Ireland;4180,Italy;
# 3160,Japan;4190,Luxembourg;4210,Netherlands;4220,Norway;4230,Poland;
# 4240,Portugal;5150,New Zealand;4280,Spain;4290,Sweden;4300,Switzerland;
# 4308,United Kingdom;2450,United States of America

# 4018 Belarus; 4188 Lithuania; 4272 Russian Federation; 4055 Estonia; 4186 Latvia
# 4303 Ukraine; 4260 Republic of Moldova

mycountry <- c(5020,4010,4020,2090,4050,4070,4080,4140,4150,4160,4170,4180,3160,
               4190,4210,4220,4230,4240,5150,4280,4290,4300,4308,2450)
mycountry = c(mycountry, 4018, 4188, 4272, 4055, 4186, 4303, 4260)

ICD10 <- ICD10[which(ICD10$Country %in% mycountry),]

colnames(ICD10)[10:35] <- c("All","Age0","Age1","Age2","Age3","Age4","Age5-9","Age10-14","Age15-19","Age20-24","Age25-29","Age30-34","Age35-39","Age40-44","Age45-49","Age50-54","Age55-59","Age60-64","Age65-69","Age70-74","Age75-79","Age80-84","Age85-89","Age90-94","Age95+","AgeUNK")

#N.B. Codes UExx and CHxx refer to Portugal 2004-2005
##From A00 to B99

##
ICD10$Cause <- as.factor(ICD10$Cause)

InfectList <- c(levels(ICD10$Cause)[(substr(levels(ICD10$Cause),1,1) %in% c("A","B"))],"UE02","UE03","UE04","UE05")
InfectList <- InfectList[-which(InfectList=="AAA")] #AAA includes all deaths

#From C00 to D48 (excluding C33-C34)
NeopList <- c(levels(ICD10$Cause)[(substr(levels(ICD10$Cause),1,1)=="C"&!(as.numeric(substr(levels(ICD10$Cause),2,3))%in% c(33,34)))],levels(ICD10$Cause)[(substr(levels(ICD10$Cause),1,1)=="D"&as.numeric(substr(levels(ICD10$Cause),2,3))<49)],"UE08","UE09","UE10","UE11","UE12","UE13","UE14","UE16","UE17","UE18","UE19","UE20","UE21","UE22","UE23","UE24")

#PROBLEM: UE15 includes also neoplasm of larynx. Beware the change (from two to three digits) between list 103 and 104
LungList <- c("C33","C34","C340","C341","C342","C343","C348","C349")

##D50-D64 Anaemias, D65-D89 remainder of diseases of blood and bllod-forming organs

#From E00 to E88
EndocList <- c(levels(ICD10$Cause)[(substr(levels(ICD10$Cause),1,1)=="E")],"CH04")
MentalList <- c(levels(ICD10$Cause)[(substr(levels(ICD10$Cause),1,1)=="F")], "CH05","R40",  "R400", "R401","R402","R41","R410","R411","R413","R418","R45","R450","R451","R453","R454","R456","R457","R458","R46",  "R460", "R461", "R463", "R468","R54","R95")
NervList <- c(levels(ICD10$Cause)[(substr(levels(ICD10$Cause),1,1)=="G")], "CH06")
SkinList <- c(levels(ICD10$Cause)[(substr(levels(ICD10$Cause),1,1)=="L")],levels(ICD10$Cause)[(substr(levels(ICD10$Cause),1,1)=="M")], "CH12","CH13")
UroList <- c(levels(ICD10$Cause)[(substr(levels(ICD10$Cause),1,1)=="N")], "CH14")
InfantList <- c(levels(ICD10$Cause)[(substr(levels(ICD10$Cause),1,1)=="O")],levels(ICD10$Cause)[(substr(levels(ICD10$Cause),1,1)=="P")],"R95" "CH15","CH16")
CongeList <- c(levels(ICD10$Cause)[(substr(levels(ICD10$Cause),1,1)=="Q")], "CH17")
##F01-F99 Mental and behavioural disorders
##G00-G98 Diseases of the nervous system
##H00-H57 Diseases of the eye and adnexa, H60-H93 Diseases of the ear and mastoid process

#From I00 to I99
CircList <- c(levels(ICD10$Cause)[(substr(levels(ICD10$Cause),1,1)=="I")],"CH09")

#From J00 to J98
RespList <- c(levels(ICD10$Cause)[(substr(levels(ICD10$Cause),1,1)=="J")],"CH10")

#From K00 to K93
DigestList <- c(levels(ICD10$Cause)[(substr(levels(ICD10$Cause),1,1)=="K")],"UE43","UE44")
##L00-L98 Diseases of the skin and subcutaneous tissue
##M00-M99 Diseases of the musculoskeletal system and connective tissue
##N00-N98 Diseases of the genitourinary system
##O00-O99 Pregnancy, childbirth and the puerperium
##P00-P96 Certain conditions originating in the perinatal period
##Q00-Q99 Congenital malformations, deformations and chromosomal abnormalities

#From V01 to Y89
ExternList <- c(levels(ICD10$Cause)[(substr(levels(ICD10$Cause),1,1)%in% c("V", "W", "X", "Y"))],"CH20")

#From R00 to R99
OthList <- c(levels(ICD10$Cause)[(substr(levels(ICD10$Cause),1,1)=="R")],"CH18")[-c(121:129,139:150,174)]

ICD10$CauseM <- NA
ICD10$CauseM <- ifelse(ICD10$Cause%in%InfectList, 1,ICD10$CauseM)
ICD10$CauseM <- ifelse(ICD10$Cause%in%NeopList, 2,ICD10$CauseM)
ICD10$CauseM <- ifelse(ICD10$Cause%in%LungList, 3,ICD10$CauseM)
ICD10$CauseM <- ifelse(ICD10$Cause%in%EndocList, 4,ICD10$CauseM)
ICD10$CauseM <- ifelse(ICD10$Cause%in%CircList, 5,ICD10$CauseM)
ICD10$CauseM <- ifelse(ICD10$Cause%in%RespList, 6,ICD10$CauseM)
ICD10$CauseM <- ifelse(ICD10$Cause%in%DigestList, 7,ICD10$CauseM)
ICD10$CauseM <- ifelse(ICD10$Cause%in%ExternList, 8,ICD10$CauseM)
ICD10$CauseM <- ifelse(ICD10$Cause%in%MentalList, 9,ICD10$CauseM)
ICD10$CauseM <- ifelse(ICD10$Cause%in%NervList, 10,ICD10$CauseM)
ICD10$CauseM <- ifelse(ICD10$Cause%in%UroList, 11,ICD10$CauseM)
ICD10$CauseM <- ifelse(ICD10$Cause%in%SkinList, 12,ICD10$CauseM)
ICD10$CauseM <- ifelse(ICD10$Cause%in%CongeList, 13,ICD10$CauseM)
ICD10$CauseM <- ifelse(ICD10$Cause%in%InfantList, 14,ICD10$CauseM)
ICD10$CauseM <- ifelse(ICD10$Cause%in%OthList, 15,ICD10$CauseM)


#ICD10$CauseM <- ifelse(is.na(ICD10$CauseM), 9,ICD10$CauseM)
#Austria: Move diabete mellitus with circulatory complications to "Circulatory" Diseases, to be consistent with ICD 9 classification
#ICD10$CauseM[ICD10$Country==4010&ICD10$Cause=="E145"] <- 5

AllDeaths <- ICD10[ICD10$Cause%in% c("AAA","CH00"),c(1,4,7,10:35)]
SpecDeaths <- ICD10[!is.na(ICD10$CauseM),]
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
        Agg2 <- aggregate(final2[,6], by=final2[,c("Age","Sex")],sum)
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

write.csv(complete,file="ICD10_comp_new2.csv",sep=",")        
