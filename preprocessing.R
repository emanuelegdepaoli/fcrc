################################################################################
# preprocessing
# edit 01/02/2022: added ex-URSS countries, Stefano's modifications ICD data and all causes (7, 12, 12)
# edit 22/04/2022: added class 5-39
################################################################################

# country codes
# 5020,Australia;4010,Austria;4020,Belgium;2090,Canada;4050,Denmark;4070,Finland;
# 4080,France;4140,Greece;4150,Hungary;4160,Iceland;4170,Ireland;4180,Italy;
# 3160,Japan;4190,Luxembourg;4210,Netherlands;4220,Norway;4230,Poland;5150,New Zealand; 4240 Portugal;
# 4280,Spain;4290,Sweden;4300,Switzerland;4308,United Kingdom;2450,United States of America
# 4018 Belarus; 4188 Lithuania; 4272 Russian Federation; 4055 Estonia; 4186 Latvia
# 4303 Ukraine; 4260 Republic of Moldova

require(dplyr)
require(tidyr)
require(reshape2)
require(ggplot2)

################################################################################
# causes of death
################################################################################

deaths_who = read.csv("data/data_WHO/morticd07/ICD7_comp_new2.csv", header=TRUE)
deaths_who = rbind(deaths_who, read.csv("data/data_WHO/morticd08/ICD8_comp_new2.csv", header=TRUE))
deaths_who = rbind(deaths_who, read.csv("data/data_WHO/morticd9/ICD9_comp_new2.csv", header=TRUE))
deaths_who = rbind(deaths_who, read.csv("data/data_WHO/morticd10/ICD10_comp_new2.csv", header=TRUE))

deaths_who$Country = as.factor(deaths_who$Country)
countries = c("CAN", "USA", "JPN", "AUT", "BLR", "BEL", "DNK", "EST", "FIN", "FRA", "GRC", "HUN",
               "ISL", "IRL", "ITA", "LVA", "LTU", "LUX", "NLD", "NOR", "POL", "PRT", "MDA", 
               "RUS", "ESP", "SWE", "CHE", "UKR", "UK", "AUS", "NZL")
levels(deaths_who$Country) = countries
deaths_who$CauseM = as.factor(deaths_who$CauseM)
levels(deaths_who$CauseM) = c("INFE", "NEOP", "LUNG", "END", "CIRC",
                         "RESP", "DIG", "EXT", "MENT", "NERV",
                         "UROG", "SKIN", "CONG", "INFA", "OTHER.SIM",
                         "OTHERS")
# sex 9 not specified
deaths_who$Sex = factor(deaths_who$Sex, levels = c(1, 2, 9), labels = c("M", "F", "NS"))

# select data 
# age-groups 0-4, 5-39, 40-65, 65+
# years 1959-2015
deaths_who$AgeGroup[deaths_who$Age %in% c("0", "1", "2", "3", "4")] = 1
deaths_who$AgeGroup[deaths_who$Age %in% c("5-9", "10-14", "15-19", "20-24", "25-29",
                                          "30-34", "35-39")] = 2
deaths_who$AgeGroup[deaths_who$Age %in% c("40-44", "45-49", "50-54", "55-59", "60-64")] = 3
deaths_who$AgeGroup[deaths_who$Age %in% c("65-69", "70-74", "75-79", "80-84",
                              "85-89", "90-94", "95+")] = 4
table(deaths_who$AgeGroup)
deaths2_who = deaths_who %>% filter(!(AgeGroup %in% c("All", "UNK"))) %>%
        group_by(Country, Year, Sex, AgeGroup, CauseM) %>%
        filter(Year >= 1959 & Year <= 2015 & Sex != "NS") %>%
        summarise(Counts = if (sum(is.na(Counts)) == n()) NA else sum(Counts, na.rm = T),
                  AllCauses = if (sum(is.na(AllCauses)) == n()) NA else sum(AllCauses, na.rm = T))
deaths2_who$AgeGroup = as.factor(deaths2_who$AgeGroup)

levels(deaths2_who$AgeGroup) = c("0-4", "5-39", "40-64", "65+")   
deaths2_who$Sex = droplevels(deaths2_who$Sex)

# check sum(Counts) = Allcauses
check = deaths2_who %>% group_by(Country, Year, Sex, AgeGroup) %>% mutate(check = sum(Counts))
check[(check$check-check$AllCauses)!=0,]
check$Country[(check$check-check$AllCauses)!=0]
rm(check)
# POR 1969 0-4 age has some anomalies, but POR will be removed in the next steps 

# check missing data 
nyears_bycountry = deaths2_who %>% group_by(Country) %>% summarise(n = length(unique(Year)))
nyears_bycountry$Country[which(nyears_bycountry$n != 57)]
sapply(nyears_bycountry$Country[which(nyears_bycountry$n != 57)], 
       function(x) paste(x, (1959:2015)[!(1959:2015 %in% deaths2_who$Year[deaths2_who$Country == x])]))

# 1959-1970 LUX
# 1997-1998 POL
# 1971-1979, 2002-2015 POR
# 2000 UK
# 2005 AUS
# BLR, EST, LVA, LTU, MDA, RUS, UKR many years missing from who database


# remove POR, LUX, GRC (missing life expectancy before 1981)
# remove ISL (noisy data)
# remove BLR, EST, LVA, LTU, MDA, RUS, UKR 
deaths2_who = deaths2_who %>% filter(!(Country %in% c("GRC", "PRT", "LUX", "ISL", 
                                              "BLR", "EST", "LVA", "LTU", "MDA",
                                              "RUS", "UKR"))) 
deaths2_who$Country = droplevels(deaths2_who$Country)
countries = levels(deaths2_who$Country)

# use date of BLR, EST, LVA, LTU, RUS, UKR from HCoD database 
deaths_hcod = read.csv("data/data_HcoD/Eastern_HCoD.csv")
deaths_hcod = deaths_hcod[,-1]
colnames(deaths_hcod) = c("Country", "Year", "Sex", "Age", "CauseM", "Counts", "AllCauses")
deaths_hcod$CauseM = factor(deaths_hcod$CauseM, levels = c(1:14, 16), 
                            labels = c(levels(deaths2_who$CauseM)[1:14], "OTHERS"))
deaths_hcod$Sex = factor(deaths_hcod$Sex, levels = 1:3, labels = c("M", "F", "MF"))
deaths_hcod$AgeGroup[deaths_hcod$Age %in% c("0", "1-4")] = 1
deaths_hcod$AgeGroup[deaths_hcod$Age %in% c("5-9", "10-14", "15-19", "20-24", "25-29",
                                            "30-34", "35-39")] = 2
deaths_hcod$AgeGroup[deaths_hcod$Age %in% c("40-44", "45-49", "50-54", "55-59", "60-64")] = 3
deaths_hcod$AgeGroup[deaths_hcod$Age %in% c("65-69", "70-74", "75-79", "80-84",
                                       "85+")] = 4
table(deaths_hcod$AgeGroup, deaths_hcod$Age)

# for HCoD, counts have decimals!
deaths2_hcod = deaths_hcod %>% 
  group_by(Country, Year, Sex, AgeGroup, CauseM) %>%
  filter(Year >= 1959 & Year <= 2015 & Sex != "MF") %>%
  summarise(Counts = if (sum(is.na(Counts)) == n()) NA else sum(round(Counts), na.rm = T),
            AllCauses = if (sum(is.na(AllCauses)) == n()) NA else sum(round(AllCauses), na.rm = T))
deaths2_hcod$AgeGroup = as.factor(deaths2_hcod$AgeGroup)

levels(deaths2_hcod$AgeGroup) = c("0-4", "5-39", "40-64", "65+")   
deaths2_hcod$Sex = droplevels(deaths2_hcod$Sex)

# check NA
table(deaths2_hcod$Country[is.na(deaths2_hcod$Counts)], 
      deaths2_hcod$CauseM[is.na(deaths2_hcod$Counts)],
      deaths2_hcod$AgeGroup[is.na(deaths2_hcod$Counts)],
      deaths2_hcod$Sex[is.na(deaths2_hcod$Counts)])

# check missing years
nyears_bycountry = deaths2_hcod %>% group_by(Country) %>% summarise(min = min(unique(Year)), 
                                                                    max = max(unique(Year)))
nyears_bycountry
sapply(nyears_bycountry$Country, 
       function(x) paste(x, (1965:2012)[!(1965:2012 %in% deaths2_hcod$Year[deaths2_hcod$Country == x])]))

# merge, years from 1965 to 2012
# exclude BLR 
deaths_who$Country = as.character(deaths_who$Country)
deaths2 = rbind(deaths2_who, deaths2_hcod) %>% filter(Year >= 1965 & Year <= 2012 & Country != "BLR")
deaths2$Country = as.factor(deaths2$Country)

# impute missing values POL, UK, AUS
# use linear interpolation since the integrals are approximated by trapezoidal rule
deaths2 = as.data.frame(deaths2)
range_years = 1965:2012
n_years = length(range_years)
n_countries = length(unique(deaths2$Country))
rownames(deaths2) = NULL
for(country in c("POL", "UK", "AUS")){
  x = unique(deaths2$Year[deaths2$Country == country])
  for(sex in levels(deaths2$Sex)){
    for(age in levels(deaths2$AgeGroup)){
      for(cause in levels(deaths2$CauseM)){
        idx = deaths2$Country == country & 
          deaths2$Sex==sex & deaths2$AgeGroup == age & deaths2$CauseM == cause
        years = range_years[!range_years%in%x]
        y = rep(NA, length(range_years))
        y[range_years%in%x] = deaths2$Counts[idx]
        #smooth_obj = smooth.spline(x = range_years[range_years%in%x], y[range_years%in%x])
        #y_est = predict(smooth_obj, range_years[!range_years%in%x])$y
        y_est = imputeTS::na_interpolation(y, option = "linear")[!range_years%in%x] # for linear interpolation
        for(j in 1:length(years)) deaths2[nrow(deaths2)+1,] = list(country, years[j], sex, age, cause, 
                                               round(y_est[j]), NA)
      }
    }
  }
}

deaths2 = deaths2 %>% group_by(Country, Year, Sex, AgeGroup) %>%
  mutate(AllCauses = sum(Counts))

# check sum(Counts) = Allcauses
check = deaths2 %>% group_by(Country, Year, Sex, AgeGroup) %>% mutate(check = sum(Counts))
check[(check$check-check$AllCauses)!=0,] 
rm(check)
# check missing data 
nyears_bycountry = deaths2 %>% group_by(Country) %>% summarise(n = length(unique(Year)))
nyears_bycountry$Country[which(nyears_bycountry$n != length(range_years))]
# check zeros
check = deaths2[which(deaths2$Counts==0),]
table(check$CauseM, check$AgeGroup, check$Sex)
x = table(check$CauseM, check$AgeGroup, check$Sex, check$Country)
for(j in 1:n_countries) {
  print(levels(check$Country)[j])
  print(x[,3,,j])
}
rm(check, x)
# zero counts
# 0-4 Infect., Neoplasm, Respir., External, Nervous
# 40-64, 65+ Mental 

# reasonable number of zeros to impute for 0-4 Infect., Neoplasm, Respir., External,
# Nervous, Congenital, Infant

# Excluding Lung, Urog, Skin, imputing Mental for 5-39

# select causes of death 
# allcauses_sub is the sum of counts for selected causes
inf_causes = c("CONG", "INFA")
other_causes = c("OTHER.SIM", "OTHERS")
# code for old selected causes of death 
# deaths3 = deaths2 %>% 
#   filter((AgeGroup == "0-4" & (CauseM %in% c(inf_causes, "Infect.")) |
#             (AgeGroup == "40-64" & !(CauseM %in% c(inf_causes, other_causes, "Mental", "Skin"))) |
#             (AgeGroup == "65+" & !(CauseM %in% c(inf_causes, "Infect.", other_causes, "Skin"))))) %>%
#   mutate(AllCauses_sub = sum(Counts))
deaths3 = deaths2 %>% 
  filter((AgeGroup == "0-4" & (CauseM %in% c(inf_causes, "INFE", "NEOP", "RESP",
                                             "EXT", "NERV")) |
            (AgeGroup == "5-39" & !(CauseM %in% c("LUNG", "UROG", "SKIN", 
                                    other_causes, inf_causes))) |
            (AgeGroup == "40-64" & !(CauseM %in% c(other_causes, inf_causes))) |
            (AgeGroup == "65+" & !(CauseM %in% c(other_causes, inf_causes))))) %>%
  mutate(AllCauses_sub = sum(Counts))

# check sum(Counts) = Allcauses_sub
check = deaths3 %>% group_by(Country, Year, Sex, AgeGroup) %>% mutate(check = sum(Counts))
check[(check$check-check$AllCauses_sub)!=0,] 
rm(check)

# check zeros
check = deaths3[which(deaths3$Counts==0),]
x = table(check$CauseM, check$AgeGroup, check$Sex, check$Country)
for(j in which(apply(x, 4, sum)!=0)) {
  print(countries[j])
  print(x[,,,j])
}
rm(check, x)
deaths3$Counts[deaths3$Counts == 0] = 0.5 
# the zero count is replaced by the maximum rounding error of 0.5 commonly used 
# in compositional and microbiome data analysis [Aitchison (2003), Kurtz et al. (2015)].

# normalize
deaths3$Prop = deaths3$Counts/deaths3$AllCauses_sub

# datasets in long format for analysis - men
# countries automatically sorted in alphabetical order 
long_deaths_M_props = dcast(deaths3[deaths3$Sex == "M",], 
                    Country + Year ~ AgeGroup + CauseM, value.var = "Prop")
long_deaths_F_props = dcast(deaths3[deaths3$Sex == "F",], 
                    Country + Year ~ AgeGroup + CauseM, value.var = "Prop")
long_deaths_M_counts = dcast(deaths3[deaths3$Sex == "M",], 
                            Country + Year ~ AgeGroup + CauseM, value.var = "Counts")
long_deaths_F_counts = dcast(deaths3[deaths3$Sex == "F",], 
                            Country + Year ~ AgeGroup + CauseM, value.var = "Counts")

Z_M = log(aperm(`dim<-`(t(long_deaths_M_props[,-(1:2)]), c(ncol(long_deaths_M_props)-2, n_years, n_countries)),
             c(2, 1, 3)))
Z_F = log(aperm(`dim<-`(t(long_deaths_F_props[,-(1:2)]), c(ncol(long_deaths_F_props)-2, n_years, n_countries)),
            c(2, 1, 3)))
X_M = aperm(`dim<-`(t(long_deaths_M_counts[,-(1:2)]), c(ncol(long_deaths_M_counts)-2, n_years, n_countries)),
                c(2, 1, 3))
X_F = aperm(`dim<-`(t(long_deaths_F_counts[,-(1:2)]), c(ncol(long_deaths_F_counts)-2, n_years, n_countries)),
                c(2, 1, 3))

################################################################################
# life expectancy
################################################################################

filenames = list.files("data/data_lifexp/lifexp/E0per", pattern = "*.txt", full.names=TRUE)
ldf = lapply(filenames, function(x) read.table(x, skip = 1, header = T))
data = do.call(rbind, ldf)
countries_lifexp = sapply(filenames, function(x) gsub(".*[/]([^.]+)[.].*", "\\1", x))
rows_per_countries = sapply(ldf, nrow)
data[,5] = rep(countries_lifexp, rows_per_countries)
colnames(data)[5] = "Country"
data$Female = as.numeric(data$Female)
data$Male = as.numeric(data$Male)
data$Total = as.numeric(data$Total)

# FRATNP total pop; GBR_NP all UK; NZL_NP total pop
levels(deaths3$Country)
countries_lifexp = c("AUS", "AUT", "BEL", "CAN", "CHE", "DNK", "ESP", "EST", "FIN",
                     "FRATNP", "HUN", "IRL", "ITA", "JPN", "LTU", "LVA", "NLD", "NOR",
                     "NZL_NP", "POL", "RUS", "SWE", "GBR_NP", "UKR", "USA")
lifexp = data %>% filter(Country %in% countries_lifexp & Year >= 1965 & Year <= 2012)
lifexp$Country[lifexp$Country == "FRATNP"] = "FRA"
lifexp$Country[lifexp$Country == "GBR_NP"] = "UK"
lifexp$Country[lifexp$Country == "NZL_NP"] = "NZL"
lifexp = arrange(lifexp, Country)
lifexp$Country = as.factor(lifexp$Country)
levels(lifexp$Country)
# check correspondence between deaths and lifexp datasets
levels(deaths3$Country) == levels(lifexp$Country)
sum(!(long_deaths_M_counts$Country == lifexp$Country))

print(lifexp %>% group_by(Country) %>% summarise(n()), n = 30)

# GRC starts from 1981
# LUX starts from 1960
# NZL ends in 2013 

# impute missing values NZL
# y1 = rev(na_interpolation(c(lifexp$Female[lifexp$Country == "NZL"], rep(NA, 2)), "linear"))[1:2]
# y2 = rev(na_interpolation(c(lifexp$Male[lifexp$Country == "NZL"], rep(NA, 2)), "linear"))[1:2]
# y3 = rev(na_interpolation(c(lifexp$Total[lifexp$Country == "NZL"], rep(NA, 2)), "linear"))[1:2]
# lifexp[nrow(lifexp)+1,] = list(2014, y1[1], y2[1], y3[1], "NZL")
# lifexp[nrow(lifexp)+1,] = list(2015, y1[1], y2[1], y3[1], "NZL")
# lifexp = lifexp  %>% arrange(Country, Year)

lifexp_M = lifexp$Male
lifexp_F = lifexp$Female
Y_M = matrix(lifexp_M, ncol = n_countries)
Y_F = matrix(lifexp_F, ncol = n_countries)
causes = colnames(long_deaths_M_counts[,-(1:2)])
countries = levels(lifexp$Country)

save(Z_M, Z_F, X_M, X_F, Y_M, Y_F, deaths3, deaths2, causes, countries, range_years,
     file = "data_fullCOD.RData")
