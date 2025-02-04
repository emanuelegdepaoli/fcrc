ICD8 <- read.table("Morticd8",sep=",",header=TRUE)
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
ICD8 <- ICD8[which(ICD8$Country %in% mycountry),]

colnames(ICD8)[10:35] <- c("All","Age0","Age1","Age2","Age3","Age4","Age5-9","Age10-14","Age15-19","Age20-24","Age25-29","Age30-34","Age35-39","Age40-44","Age45-49","Age50-54","Age55-59","Age60-64","Age65-69","Age70-74","Age75-79","Age80-84","Age85-89","Age90-94","Age95+","AgeUNK")

InfectList <- c("A001", "A002", "A003", "A004", "A005", "A006", "A007", "A008", "A009", "A010", "A011", "A012","A013", "A014", "A015", "A016", "A017", "A018", "A019", "A020", "A021", "A022", "A023", "A024", "A025", "A026", "A027", "A028","A029", "A030", "A031", "A032", "A033", "A034", "A035", "A036", "A037", "A038", "A039", "A040", "A041", "A042", "A043", "A044")
NeopList <- c("A045", "A046", "A047", "A048", "A049","A050", "A052", "A053", "A054", "A055", "A056", "A057", "A058", "A059", "A060","A061")
LungList <- c("A051")
EndocList <- c("A062","A063", "A064", "A065", "A066","A067","A068")
MentalList <- c("A069","A070","A071","A136")
NervList <- c("A072","A073","A074","A078","A079")

##A067 Anaemias, A068 Other diseases of blood and blood-forming organs -> Endoc
##A069 Psychoses, A070 Neuroses -> Mental
##A071 Mental retardation -> Mental
##A072 Meningitis, A073 Multiple sclerosis, A074 Epilepsy -> Nervous
##A075-A077 Inflammatory diseases of eye, Cataract, Glaucoma, -> Other
##A078 Otitis media and mastoiditis -> Other
##A079 Other diseases of nervous system and sense organs -> Nervous
CircList <- c("A080", "A081", "A082", "A083", "A084", "A085", "A086","A087","A088")
RespList <- c("A089", "A090", "A091", "A092", "A093", "A094", "A095", "A096")
DigestList <- c("A097","A098", "A099", "A100", "A101", "A102", "A103", "A104")
UroList <- c("A105","A106","A107","A108","A109","A110","A111","B039")
SkinList <- c("A119","A120","A121","A122","A123","A124","A125")
CongeList <- c("A126","A127","A128","A129","A130")
InfantList <- c("A112","A113","A114","A115","A116","A117","A118","A131","A132","A133","A134","A135")

##A105-A106 Nephritis and nephrosis -> Uro
##A107 Infections of kidney -> Uro
##A108 Calculus of urinary system, A109 Hyperplasia of prostate, A110 Diseases of breast, A111 Other diseases of genito-urinary system -> Uro
##A112 Toxaemias of pregnancy and the puerperium, A113 Haemorrhage of pregnancy and childbirth -> Infant
##A114 Abortion induced for legal indications, A115 Other and unspecified abortion -> Infant
##A116-A117 Other complications of pregnancy, childbirth and the puerperium -> Infant
##A118 Delivery without mention of complication -> Infant
##A119-A120 Infections of skin -> Skin
##A121 Arthritis and spondylitis, A122 Non-articular rheumatism and rheumatism unspecified, A123 Osteomyelitis and periostitis, -> Skin
##A124 -A125 Other diseases of musculoskeletal system and connective tissue -> Skin
##A126- A130 Congenital anomalies -> Congenital
##A131 - A135 causes of perinatal morbidity and mortality -> Infant
##A136 Senility without mention of psychosis -> Mental
ExternList <- c("A138", "A139", "A140","A141", "A142", "A143", "A144", "A145", "A146", "A147", "A148", "A149", "A150","S47")
OthList <- c("A075","A076","A077","A137","B045")

ICD8$CauseM <- NA
ICD8$CauseM <- ifelse(ICD8$Cause%in%InfectList, 1,ICD8$CauseM)
ICD8$CauseM <- ifelse(ICD8$Cause%in%NeopList, 2,ICD8$CauseM)
ICD8$CauseM <- ifelse(ICD8$Cause%in%LungList, 3,ICD8$CauseM)
ICD8$CauseM <- ifelse(ICD8$Cause%in%EndocList, 4,ICD8$CauseM)
ICD8$CauseM <- ifelse(ICD8$Cause%in%CircList, 5,ICD8$CauseM)
ICD8$CauseM <- ifelse(ICD8$Cause%in%RespList, 6,ICD8$CauseM)
ICD8$CauseM <- ifelse(ICD8$Cause%in%DigestList, 7,ICD8$CauseM)
ICD8$CauseM <- ifelse(ICD8$Cause%in%ExternList, 8,ICD8$CauseM)
ICD8$CauseM <- ifelse(ICD8$Cause%in%MentalList, 9,ICD8$CauseM)
ICD8$CauseM <- ifelse(ICD8$Cause%in%NervList, 10,ICD8$CauseM)
ICD8$CauseM <- ifelse(ICD8$Cause%in%UroList, 11,ICD8$CauseM)
ICD8$CauseM <- ifelse(ICD8$Cause%in%SkinList, 12,ICD8$CauseM)
ICD8$CauseM <- ifelse(ICD8$Cause%in%CongeList, 13,ICD8$CauseM)
ICD8$CauseM <- ifelse(ICD8$Cause%in%InfantList, 14,ICD8$CauseM)
ICD8$CauseM <- ifelse(ICD8$Cause%in%OthList, 15,ICD8$CauseM)

AllDeaths <- ICD8[ICD8$Cause=="A000",c(1,4,7,10:35)]
SpecDeaths <- ICD8[!is.na(ICD8$CauseM),]
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

write.csv(complete,file="ICD8_comp_new2.csv",sep=",")        
