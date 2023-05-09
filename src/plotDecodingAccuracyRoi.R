rm(list=ls()) #clean console


library(ggplot2)
library(doBy)
library(cowplot)
library(Rmisc)
library(stringr)

# analysis libraries
library(rstatix)
library(dplyr)
library(car)

library(ez) # anova
library(schoRsch)

library(afex) # mixed models
#######################################################

pathResults <- '/Users/battal/Cerens_files/fMRI/Processed/RhythmCateg/RhythmBlock/rhythmBlock_derivatives_cosmoMvpa/'

########
mvpa <- read.csv(paste(pathResults, 'RhythmBlockDecoding_jubrainatlas_s2_ratio150_202305091421.csv', sep ='/'))



mvpa <- mvpa[-c(7:9)]
mvpa$subID <-as.factor(mvpa$subID)

head(mvpa)

voxelSize = '150'

#### order things a bit for easy manipulation for plotting


#name change from mask to roi
mvpa$mask <-as.factor(mvpa$mask)
names(mvpa)[2] <- 'roi'

# let's make a expType to split the no_pitch exp from pitch exp
mvpa$subNb <- as.numeric(mvpa$subID)
mvpa$expType<- ifelse(mvpa$subNb < 23, 'P4', 'P1') # 23 for RhythmBlock, 11 for Nonmetric
mvpa$expType <-as.factor(mvpa$expType)


mvpa$roi_order <- ifelse(mvpa$roi == 'STS1', 1,
                         ifelse(mvpa$roi == 'STS2', 2, 
                                ifelse(mvpa$roi == 'TE10', 3, 
                                       ifelse(mvpa$roi == 'TE11', 4,
                                              ifelse(mvpa$roi == 'TE12', 5,
                                                     ifelse(mvpa$roi == 'TE21', 6,
                                                            ifelse(mvpa$roi == 'TE22', 7,
                                                                   ifelse(mvpa$roi == 'TE30', 8, 
                                                                          ifelse(mvpa$roi == 'TI', 9, 
                                                                                 ifelse(mvpa$roi == 'TPJ', 10, 11))))))))))

# think about other ways of ordering with below function
# mvpa$roi_order <- grepl('nopitch', mvpa$roi, fixed = TRUE)

# currently we  don't have 2 hemispheres for every ROI
# later on consider this separation for "Aud Cx only" analysis
# mvpa$hemis <- ifelse(substr(mvpa$roi,1,1) == 'l', 'left',
#                      ifelse(substr(mvpa$roi,1,1) == 'r', 'right',NA))

##### change the hemisphere header/name here - 15.03/2023
colnames(mvpa)[7] <- 'hemis'

# make everything factor
str(mvpa)
mvpa$image<- as.factor(mvpa$image)
mvpa$subType<-as.factor(mvpa$subType)
mvpa$subTypeLabel<- as.factor(mvpa$subTypeLabel)
mvpa$hemis<-as.factor(mvpa$hemis)

# subset the dataframe for plotting/analysis
img = 't_maps' # or 't_maps' 'beta'
exp = 'P1' 

subsetmvpa = subset(mvpa,expType == exp)
subsetmvpa = subset(subsetmvpa, image == img)

str(subsetmvpa)


# let's subset it again with only 8 rois
#subsetmvpa<- subset(subsetmvpa, roi_order < 9)

df <- summarySE(data = subsetmvpa, 
                groupvars=c('roi_order','roi', 'subType','hemis'),
                measurevar='accuracy', na.rm = TRUE)
df



######### analyze the data ######### 
# let's do t-test 

# separate for ROIs
head(subsetmvpa)

subsetmvpa$roi[11][]
# Levels: STS1 STS2 TE10 TE11 TE12 TE21 TE22 TE30 TI TPJ TeI
roiName = 'TI'
lSTG <-subset(subsetmvpa, roi == roiName & hemis == 'l')
rSTG <-subset(subsetmvpa, roi ==roiName & hemis == 'r')

t.test(lSTG$accuracy, mu = 0.5, alternative = 'greater') # t = 0.19825, df = 9, p-value = 0.4236
t.test(rSTG$accuracy, mu = 0.5, alternative = 'greater')


############# PLOTTING  #############
#is.na(mvpa$accuracy)

min(subsetmvpa$accuracy, na.rm = TRUE)
max(subsetmvpa$accuracy, na.rm = TRUE)

setlimit = c(0.2,0.9) 
setbreak = c(0.25, 0.5, 0.75, 1)


rois = '11ROI'

shapesize = 1
shapetype = 21
shapestroke = 1
transparent = 1 #0.6
jitter  = position_jitterdodge(0.3) # position_jitter(width=0.3)



# colors 
nonmetricGrayBad = "#3d8c55ff" # complex green= 3d8c55ff, nonmetricGrap = 6B6B6B
nonmetricGrayGood = "#3d8c55ff" 
simplePurpleBad = "#c4a4c9"
simplePurpleGood = "#8a4a95"
# conditions
category2 = ''
category1 = ''

cond1= paste0(category2," Bad Tapper")
cond2 = paste0(category2,' Good Tapper')
cond3 = paste0(category1,' Bad Tapper')
cond4 = paste0(category1,' Good Tapper')

cond = 'Complex'

# ##### separate tappers
# fig <- ggplot(data = subsetmvpa, 
#               aes(x = reorder(roi, roi_order),
#                   y = accuracy, 
#                   color = subType,
#                   group = subType)) +
#   geom_point(data=subsetmvpa, 
#              aes(x = reorder(roi, roi_order), 
#                  y = accuracy), 
#              size = shapesize,
#              position = jitter, shape = shapetype, stroke = shapestroke, na.rm = TRUE) + 
#   stat_summary(aes(color = subType), fun=mean, fun.min = mean, fun.max = mean, 
#                geom="crossbar", size=0.6, width=0.6, position = position_dodge(width=.75), 
#                na.rm = TRUE) +
#   theme_classic() +
#   geom_errorbar(data = df, 
#                 aes(ymin = accuracy-se, ymax = accuracy+se, group = subType), 
#                 color = 'black',size=0.5, width=0.15, alpha = transparent, position = position_dodge(width=.75)) +
#   geom_hline(yintercept=c(.5), linetype="dotted", colour="black", size=.5) +
#   
#   ggtitle("") +
#   ylab("") +
#   xlab("") +
#   theme(axis.text.x=element_text(size=8, face = 'bold', angle=0, colour='black')) + # face = 'bold', 
#   theme(axis.text.y=element_text(size=12, angle=0, colour='black')) +
#   theme(axis.title.y=element_text(size=10, angle=90, colour='black')) +
#   scale_y_continuous(limits=setlimit, breaks=setbreak, position="left") +
#   scale_x_discrete(labels = c("lSTG All", "lSTG 10","lSTG 15","lSTG 20","rSTG All","rSTG 10","rSTG 15","rSTG 20"))+
#   scale_color_manual(name = '', labels = c(cond1, cond2), values=c(nonmetricGrayBad, nonmetricGrayGood)) + 
#   # theme(legend.position= "none")
#   theme(legend.position= c(.85, .85)) +
#   theme(legend.text=element_text(size=8)) +
#   theme(legend.title=element_text(size=9))
# fig
# 
# filename <- paste0(pathResults, 'Decoding_Simple_vs_',cond,  'voxelNb-', voxelSize, '_', rois,'tappers.png')
# # ggsave(filename, fig, dpi=300, width=15, height=6, units='cm') 
# ggsave(filename, fig, dpi=300, width=6, height=3) # 1024 x 512


#######################################################
########## ignore tappers #######################################################
#######################################################
df2 <- summarySE(data = subsetmvpa, 
                 groupvars=c('roi_order','roi', 'hemis'),
                 measurevar='accuracy', na.rm = TRUE)
df2


# # make a column for roi x hemis
# roiHemis <- paste(subsetmvpa$roi, subsetmvpa$hemis)
# 
# df2 <- summarySE(data = subsetmvpa, 
#                  groupvars=c('roi_order','roiHemis'),
#                  measurevar='accuracy', na.rm = TRUE)
# df2

# divide the hemispheres into 2 plots
rois = 'Right11Roi'
subsetmvpaL <-subset(subsetmvpa, hemis == 'r')

df2 <- summarySE(data = subsetmvpaL, 
                 groupvars=c('roi_order','roi'),
                 measurevar='accuracy', na.rm = TRUE)
df2


fig <- ggplot(data = subsetmvpaL, 
              aes(x = reorder(roi, roi_order),
                  y = accuracy),
              color = expType) +
  geom_jitter(size = shapesize, shape = shapetype, stroke = shapestroke, width=0.1, color = nonmetricGrayGood, 
              na.rm = TRUE) +
  stat_summary(aes(color = expType), fun=mean, fun.min = mean, fun.max = mean, geom="crossbar", size=0.6, width=0.3,
               na.rm = TRUE) +
  theme_classic() +
  geom_errorbar(data = df2, 
                aes(ymin = accuracy-se, ymax = accuracy+se), 
                color = 'black',size=0.5, width=0.15, alpha = transparent, position = position_dodge(width=.75),
                na.rm = TRUE) +
  geom_hline(yintercept=c(.5), linetype="dotted", colour="black", size=.5) +
  
  ggtitle("") +
  ylab("") +
  xlab("") +
  theme(axis.text.x=element_text(size=8, face = 'bold', angle=0, colour='black')) + # face = 'bold', 
  theme(axis.text.y=element_text(size=12, angle=0, colour='black')) +
  theme(axis.title.y=element_text(size=10, angle=90, colour='black')) +
  scale_y_continuous(limits=setlimit, breaks=setbreak, position="left") +
  #scale_x_discrete(labels = c("lSTG All", "lSTG 10","lSTG 15","lSTG 20","rSTG All","rSTG 10","rSTG 15","rSTG 20"))+
  scale_color_manual(name = '', labels =c('Simple vs. Complex'), values=c(nonmetricGrayGood)) + #
  theme(legend.position= c(.9, 1)) +
  theme(legend.text=element_text(size=9)) +
  theme(legend.title=element_text(size=9))
fig

filename <- paste0(pathResults, 'Decoding_Simple_vs_',cond,  '_', 'voxelNb-',  voxelSize,  '_', rois,  '_', exp,  '_', img, '.png')
ggsave(filename, fig, dpi=300, width=8, height=2) # 1024 x 512











# # figure for roi already cointains hemisphere info
# fig <- ggplot(data = subsetmvpa, 
#               aes(x = reorder(roi, roi_order),
#                   y = accuracy),
#               color = expType) +
#   geom_jitter(size = shapesize, shape = shapetype, stroke = shapestroke, width=0.1, color = nonmetricGrayGood, 
#               na.rm = TRUE) +
#   stat_summary(aes(color = expType), fun=mean, fun.min = mean, fun.max = mean, geom="crossbar", size=0.6, width=0.3,
#                na.rm = TRUE) +
#   theme_classic() +
#   geom_errorbar(data = df2, 
#                 aes(ymin = accuracy-se, ymax = accuracy+se), 
#                 color = 'black',size=0.5, width=0.15, alpha = transparent, position = position_dodge(width=.75),
#                 na.rm = TRUE) +
#   geom_hline(yintercept=c(.5), linetype="dotted", colour="black", size=.5) +
#   
#   ggtitle("") +
#   ylab("") +
#   xlab("") +
#   theme(axis.text.x=element_text(size=8, face = 'bold', angle=0, colour='black')) + # face = 'bold', 
#   theme(axis.text.y=element_text(size=12, angle=0, colour='black')) +
#   theme(axis.title.y=element_text(size=10, angle=90, colour='black')) +
#   scale_y_continuous(limits=setlimit, breaks=setbreak, position="left") +
#   scale_x_discrete(labels = c("lSTG All", "lSTG 10","lSTG 15","lSTG 20","rSTG All","rSTG 10","rSTG 15","rSTG 20"))+
#   scale_color_manual(name = '', labels =c('Simple vs. Nonmetric'), values=c(nonmetricGrayGood)) + #
#   theme(legend.position= c(.85, .9)) +
#   theme(legend.text=element_text(size=9)) +
#   theme(legend.title=element_text(size=9))
# fig
# 
# filename <- paste0(pathResults, 'Decoding_Simple_vs_',cond,  'voxelNb-', voxelSize,  '_', rois, '.png')
# ggsave(filename, fig, dpi=300, width=6, height=3) # 1024 x 512



