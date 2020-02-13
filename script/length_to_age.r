wd = "C:\\Users\\b9930\\Google ¶³ºÝµwºÐ\\publication\\SpatialVariability\\"
setwd(wd)


### compile data
# age-length key (ALK) data
load("data\\alk\\lemonsole_alks.Rdata")

# CPUE per length per subarea data
cpue = read.csv("data\\raw\\cpue_per_length_per_subarea\\Microstomus kitt.csv")
cpue$LngtClass_cm = cpue$LngtClass / 10
head(cpue)



### fill ALK
## fill all missing age class
fill_all_age = function(data, min_age){
    max_age = max(as.numeric(names(data)))
    
    # empty ALK for all age
    full_alk = data.frame(matrix(0, nrow=nrow(data), ncol=max_age-min_age+1))
    rownames(full_alk) = rownames(data)
    colnames(full_alk) = c(min_age:max_age)
    
    # fill empty ALK with sampled ALK
    full_alk[colnames(data)] = data
    return(full_alk)
}

# min_age = 1 for Q1, and min_age = 0 for other quarters
alk.q1.filled = lapply(alk.q1, function(data){fill_all_age(data, min_age = 1)})
alk.q3.filled = lapply(alk.q3, function(data){fill_all_age(data, min_age = 0)})


## fill ALK for min and max length
fill_min_max_length = function(data){
    # find min and max age class in ALK
    min_age = as.character(min(as.numeric(names(data))))
    max_age = as.character(max(as.numeric(names(data))))
    
    # length class having no ALK
    length_no_catch = apply(data, 1, function(x){sum(x, na.rm=TRUE)})
    
    # fill ALK for length class < min ALK and length class > max ALK
    data[1:(min(which(length_no_catch>0))-1), min_age] = 1
    data[(max(which(length_no_catch>0))+1):nrow(data), max_age] = 1
    return(data)
}

alk.q1.filled = lapply(alk.q1.filled, function(data){fill_min_max_length(data)})
alk.q3.filled = lapply(alk.q3.filled, function(data){fill_min_max_length(data)})


## fill ALK for min_length < length < max_length
fill_mid_length = function(data){
    # check if there is missing ALK for any length class
    length_no_catch = apply(data, 1, function(x){sum(x, na.rm=TRUE)})
    
    # if an age class has no ALK, fill it by the nearest length class
    if (sum(length_no_catch == 0) > 0){
        counter = 0
        for (i in 2:length(length_no_catch)){
            if (length_no_catch[i] == 0){
                counter = counter + 1
                next
            } 
            if (counter > 0){
                next_alk = i
                last_alk = next_alk - counter - 1
                if (counter %% 2 == 0){
                    num_row = counter / 2
                    data[(last_alk+1):(last_alk+num_row), ] = data[last_alk, ]
                    data[(last_alk+num_row+1):(last_alk+counter), ] = data[next_alk, ]
                } else {
                    if (counter == 1){
                        data[last_alk + counter, ] = (data[last_alk, ] + data[next_alk, ]) / 2
                    } else {
                        num_row = (counter - 1) / 2
                        num_row_mid = ceiling(counter / 2)
                        data[(last_alk+1):(last_alk+num_row), ] = data[last_alk, ]
                        data[(last_alk+num_row_mid+1):(last_alk+counter), ] = data[next_alk, ]
                        data[last_alk+num_row_mid, ] = (data[last_alk, ] + data[next_alk, ]) / 2
                    }
                }
            }
            counter = 0
        }
        return(data)
    } else {
        return(data)
    }
}

alk.q1.filled = lapply(alk.q1.filled, function(data){fill_mid_length(data)})
alk.q3.filled = lapply(alk.q3.filled, function(data){fill_mid_length(data)})

# check if all length class has ALK
# change "alk.q1.filled" to "alk.q3.filled" if needed
lapply(alk.q1.filled, function(data){  
    length_no_catch = apply(data, 1, function(x){sum(x)})
    sum(length_no_catch == 0)
})

# convert ALK to proportion
alk.q1.filled = lapply(alk.q1.filled, FUN=function(data){
    data_ratio = t(apply(data, 1, FUN=function(x){x/sum(x)}))
    return(data_ratio)
})

alk.q3.filled = lapply(alk.q3.filled, FUN=function(data){
    data_ratio = t(apply(data, 1, FUN=function(x){x/sum(x)}))
    return(data_ratio)
})



### assign age for each length class by ALK
# install.packages("FSA")
library(FSA)

# delete CPUE data that are not collected within the period of ALK
# and separate it into quarter 1 and quarter 3
cpue.q1 = subset(cpue, subset=cpue$Year %in% names(alk.q1.filled) & cpue$Quarter == 1)
cpue.q3 = subset(cpue, subset=cpue$Year %in% names(alk.q3.filled) & cpue$Quarter == 3)

# assign age for each length class
assign_age = function(alk, cpue_data){
    cpue_data$age = NA
    years = as.character(unique(cpue_data$Year))
    for (year in years){
        row_to_assign = which(cpue_data$Year==year)
        sub_cpue = cpue_data[row_to_assign, ]
        sub_key = alk[[year]]
        
        filled = alkIndivAge(key=sub_key, formula=~LngtClass_cm, data=sub_cpue)
        cpue_data[row_to_assign, ] = filled
    }
    return(cpue_data)
}

cpue.q1 = assign_age(alk.q1.filled, cpue.q1)
cpue.q3 = assign_age(alk.q3.filled, cpue.q3)

# check if all length class has been assigned age
sum(is.na(cpue.q1$age))
sum(is.na(cpue.q3$age))

# merge data on Q1 and Q3
cpue.filled = merge(cpue.q1, cpue.q3, all=TRUE)

# change format of data
library(tidyr)

cpue_agg = cpue.filled %>% spread(age, CPUE_number_per_hour)

# fill missing age class
age_class = as.numeric(na.omit(as.numeric(names(cpue_agg))))
all_age_class = seq(min(age_class), max(age_class), 1)
missing_age_class = all_age_class[!(all_age_class %in% age_class)]

for (a in missing_age_class){
    cpue_agg[, as.character(a)] = NA
}
head(cpue_agg)

# sorted
col_to_sort = !is.na(as.numeric(names(cpue_agg)))
sort_cpue_agg = cpue_agg[, col_to_sort]
sort_cpue_agg = sort_cpue_agg[, order(as.numeric(names(sort_cpue_agg)))]

cpue_agg = cbind(cpue_agg[, !col_to_sort], sort_cpue_agg)
head(cpue_agg)

# change column names
col_to_keep = is.na(as.numeric(names(cpue_agg)))
names(cpue_agg) = c(names(cpue_agg)[col_to_keep], 
                    paste0("Age_", names(cpue_agg)[!col_to_keep]))
head(cpue_agg)

sub_cpue_agg = subset(cpue_agg, select=c("Year", "Quarter", "SubArea", "Species", 
                                         grep(c("Age_"), names(cpue_agg), value=TRUE)))

library(dplyr)
output_cpue = sub_cpue_agg %>%
    group_by(Year, Quarter, SubArea, Species) %>%
    summarise_all(sum, na.rm=TRUE)
    
# save data
write.csv(output_cpue, file="data\\compiled_age_data_mk.csv", row.names=F)



### Appendix
## data mining
test = alk.q1$`2009`

# quarter 1
lapply(alk.q1, FUN=function(x){
    names(x)
})
lapply(alk.q1, FUN=function(x){
    row.names(x)
})
lapply(alk.q1, FUN=function(x){
    rowSums(x)
})

# quarter 3
lapply(alk.q3, FUN=function(x){
    names(x)
})
lapply(alk.q3, FUN=function(x){
    row.names(x)
})

#
lapply(alk.q3.filled, FUN=function(x){
    colnames(x)
})


## testing code
# test on function "fill_mid_length"
data = alk.q3.filled$`2007`
test = alk.q3.filled$`2007`

length_no_catch = apply(data, 1, function(x){sum(x)})
counter = 0
for (i in 2:length(length_no_catch)){
    if (length_no_catch[i] == 0){
        counter = counter + 1
        next
    } 
    if (counter > 0){
        next_alk = i
        last_alk = next_alk - counter - 1
        if (counter %% 2 == 0){
            num_row = counter / 2
            data[(last_alk+1):(last_alk+num_row), ] = data[last_alk, ]
            data[(last_alk+num_row+1):(last_alk+counter), ] = data[next_alk, ]
        } else {
            if (counter == 1){
                data[last_alk + counter, ] = (data[last_alk, ] + data[next_alk, ]) / 2
            } else {
                num_row = (counter - 1) / 2
                num_row_mid = ceiling(counter / 2)
                data[(last_alk+1):(last_alk+num_row), ] = data[last_alk, ]
                data[(last_alk+num_row_mid+1):(last_alk+counter), ] = data[next_alk, ]
                data[last_alk+num_row_mid, ] = (data[last_alk, ] + data[next_alk, ]) / 2
            }
        }
    }
    counter = 0
}

# test on FSA package
alk_q1_2007 = alk.q1.filled$`2007`
alk_q1_2007 = t(apply(alk_q1_2007, 1, FUN=function(x){x/sum(x)}))

cpue_q1_2007 = subset(cpue, subset=cpue$Year==2007 & cpue$Quarter==1)
max_length = max(as.numeric(rownames(alk_q1_2007)))
cpue_q1_2007_plus_group = subset(cpue_q1_2007, subset=cpue_q1_2007$LngtClass_cm > max_length)

cpue_q1_2007_zero = cpue_q1_2007[cpue_q1_2007$CPUE_number_per_hour == 0, ]
cpue_q1_2007 = cpue_q1_2007[cpue_q1_2007$CPUE_number_per_hour > 0, ]


test = alkIndivAge(key=alk_q1_2007, formula=~LngtClass_cm, data=cpue_q1_2007)

