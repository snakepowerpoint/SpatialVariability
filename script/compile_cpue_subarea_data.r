setwd("C:\\Users\\b9930\\Google ¶³ºÝµwºÐ\\publication\\SpatialVariability\\data\\raw")

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

write.csv(lgth, file="CPUE_per_length_per_subarea.csv", row.names=F, sep=",")
