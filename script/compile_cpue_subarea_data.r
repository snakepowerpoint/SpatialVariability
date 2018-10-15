wd = "C:\\Users\\b9930\\Google ¶³ºÝµwºÐ\\publication\\SpatialVariability\\"
setwd(paste0(wd, "data\\raw"))

df = read.csv("CPUE_per_length_per_subarea_AtoD.csv", header=T) # A-D
lgth = df
df = read.csv("CPUE_per_length_per_subarea_EtoI.csv", header=T) # E-I
lgth = rbind(lgth, df)
df = read.csv("CPUE_per_length_per_subarea_L.csv", header=T) # L
lgth = rbind(lgth, df)
df = read.csv("CPUE_per_length_per_subarea_M.csv", header=T) # M
lgth = rbind(lgth, df)
df = read.csv("CPUE_per_length_per_subarea_NtoO.csv", header=T) # N-O
lgth = rbind(lgth, df)
df = read.csv("CPUE_per_length_per_subarea_PtoR.csv", header=T) # P-R
lgth = rbind(lgth, df)
df = read.csv("CPUE_per_length_per_subarea_StoZ.csv", header=T) # S-Z
lgth = rbind(lgth, df)

# keep variables we want
df = subset(lgth, select=c(Year, Quarter, Subarea, Species, LngtClas, CPUE_number_per_hour))
rm(lgth)
# use data on quarter 1 and 3
df = droplevels(subset(df, subset=Quarter%in%c(1,3)))

lgth.species = split(df, f=df$Species) # split data by species

# sum up CPUE according to each year, quarter, and subarea
cpue.species = lapply(lgth.species, FUN=function(x){
    data = data.frame(x)
    with(data, aggregate(CPUE_number_per_hour, by=list(Year=Year, Quarter=Quarter, Subarea=Subarea), FUN=sum))
})

# sort data by "Year", "Quarter", and "Subarea"
cpue.species = lapply(cpue.species, FUN=function(x){
    df = data.frame(x)
    df = df[with(df, order(Year, Quarter, Subarea)), ]
})

# the maximal number of subarea among quarters is ...
max(unlist(lapply(cpue.species, FUN=function(x){
    data = data.frame(x)
    max(with(data, table(Year, Quarter)))
}))) # 183

# therefore, we drop out quarters with less than 10 subareas (5% of maximum)
cpue.species <- lapply(cpue.species, FUN=function(x){
    data1 = droplevels(subset(data.frame(x), subset=Quarter%in%1))
    data3 = droplevels(subset(data.frame(x), subset=Quarter%in%3))
    # for Q1
    yr.table1 = table(data1$Year) # summarize how many subareas within a quarter for each year
    yr.keep1 = names(yr.table1)[which(yr.table1>=10)] # keep years with more than 10 subareas
    data1 = droplevels(subset(data1, subset=Year%in%yr.keep1))
    # for Q3
    yr.table3 = table(data3$Year)
    yr.keep3 = names(yr.table3)[which(yr.table3>=10)] 
    data3 = droplevels(subset(data3, subset=Year%in%yr.keep3))
    # combind Q1 and Q3
    result = rbind(data1, data3)
    result[with(result, order(Year, Quarter, Subarea)), ]
})

# restrict CPUE data on 9 standard species during 1991-2015
species = c("Clupea harengus", "Gadus morhua", "Melanogrammus aeglefinus",
            "Merlangius merlangus", "Pleuronectes platessa", "Pollachius virens",
            "Scomber scombrus", "Sprattus sprattus", "Trisopterus esmarkii")
cpue.species = cpue.species[species]
cpue.species = lapply(cpue.species, FUN=function(x){
    data = data.frame(x)
    subset(data, subset = Year%in%c(1991:2015))
})

rm(df, lgth.species)

### delete subareas never large than 0
# totoal CPUE for each subarea
cpue.subarea = lapply(cpue.species, FUN=function(x){
    data = data.frame(x)
    data1 = with(data, aggregate(x, by = list(Subarea), FUN = sum))
    colnames(data1) = c('Subarea', 'TotalCPUE')
    data1
})

# subareas having CPUE larger than 0
cpue.subarea.non0 = lapply(cpue.subarea, FUN=function(x){
    data = data.frame(x)
    droplevels(data$Subarea[data$TotalCPUE > 0])
})

# keep subareas having CPUE larger than 0
cpue.species.drop0 = list()
for (i in 1:length(cpue.species)){
    cpue.species.drop0[[i]] = subset(data.frame(cpue.species[[i]]), 
                                     subset = Subarea %in% cpue.subarea.non0[[i]])
}
names(cpue.species.drop0) = names(cpue.species)

cpue.species.drop0 = lapply(cpue.species.drop0, FUN=function(x){
    data = data.frame(x)
    data$Subarea = droplevels(data$Subarea)
    data
})

for (i in 1:length(species)){
    cpue.species.drop0[[i]]$Species = species[i]
}

# save compiled CPUE data as a csv file
cpue.data.final = do.call(rbind, cpue.species.drop0)
rownames(cpue.data.final) = NULL
cpue.data.final$Subarea = as.character(cpue.data.final$Subarea)
cpue.data.final$Subarea = format(cpue.data.final$Subarea, scientific = FALSE) 
cpue.data.final$Species = as.character(cpue.data.final$Species)
cpue.data.final = cpue.data.final[, c("Year", "Quarter", "Subarea", "Species", "x")]
colnames(cpue.data.final) = c("Year", "Quarter", "Subarea", "Species", "CPUE")

write.csv(cpue.data.final, file=paste0(wd, "data\\compiled_cpue_data.csv"), row.names=FALSE)
