#please load function 'robust' in 'robustness_test_all_variables.r'

#load data
load("D:/Dropbox/Homework/Graduate/MyselfProposal/Data/20171029/20161001combined_compact.RData")

##########Robustness test (all variables)
###C. harengus (dim = 6, need not do robustness test)

###G. morhua (dim = 3)
scale.gm = subset(result.gm, select = c(CV.CPUE.std, Mean.age.std, Shannon.age.std,
                                        Total.CPUE.std, AMO.std, mean.temp.std, cv.temp.std))
colnames(scale.gm) = c('CV.CPUE', 'Mean.age', 'Shannon.age',
                       'Total.CPUE', 'AMO', 'mean.temp', 'cv.temp')
robust(scale.gm, 3, colnames(scale.gm)[c(3,4,5,6,7)], c(4,8,2,0,3))

#save table
john = robust(scale.gm, 3, colnames(scale.gm)[c(3,4,5,6,7)], c(4,8,2,0,3))
write.csv(john, file = 'gm.csv')

###M. aeglefinus (dim = 5)
scale.ma = subset(result.ma, select = c(CV.CPUE.std, Mean.age.std, Shannon.age.std,
                                        Total.CPUE.std, AMO.std, mean.temp.std, cv.temp.std))
colnames(scale.ma) = c('CV.CPUE', 'Mean.age', 'Shannon.age',
                       'Total.CPUE', 'AMO', 'mean.temp', 'cv.temp')
robust(scale.ma, 5, colnames(scale.ma)[c(3,4,5,6,7)], c(7,8,8,8,4))

#save table
john = robust(scale.ma, 5, colnames(scale.ma)[c(3,4,5,6,7)], c(7,8,8,8,4))
write.csv(john, file = 'ma.csv')

###M. merlangus (dim = 2, need not do robustness test)

###P. platessa (dim = 5)
scale.pp = subset(result.pp[19:48, ], select = c(CV.CPUE.std, Mean.age.std, Shannon.age.std,
                                                 Total.CPUE.std, AMO.std, mean.temp.std, cv.temp.std))
colnames(scale.pp) = c('CV.CPUE', 'Mean.age', 'Shannon.age',
                       'Total.CPUE', 'AMO', 'mean.temp', 'cv.temp')
robust(scale.pp, 5, colnames(scale.pp)[c(3,4,5,6,7)], c(3,3,8,0,2))

#save table
john = robust(scale.pp, 5, colnames(scale.pp)[c(3,4,5,6,7)], c(3,3,8,0,2))
write.csv(john, file = 'pp.csv')

###P. virens (dim = 2, need not do robustness test)

###S. scombrus (dim = 2, need not do robustness test)

###S. sprattus (dim = 2, need not do robustness test)

###T. esmarkii (dim = 4)
scale.te = subset(result.te, select = c(CV.CPUE.std, Mean.age.std, Shannon.age.std,
                                        Total.CPUE.std, AMO.std, mean.temp.std, cv.temp.std))
colnames(scale.te) = c('CV.CPUE', 'Mean.age', 'Shannon.age',
                       'Total.CPUE', 'AMO', 'mean.temp', 'cv.temp')
robust(scale.te, 4, colnames(scale.te)[c(3,4,5,6,7)], c(2,5,4,0,7))

#save table
john = robust(scale.te, 4, colnames(scale.te)[c(3,4,5,6,7)], c(2,5,4,0,7))
write.csv(john, file = 'te.csv')


