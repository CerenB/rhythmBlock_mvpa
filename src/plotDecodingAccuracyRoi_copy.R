library(ggplot2)
library(doBy)
library(cowplot)
library(Rmisc)
library(stringr)



#######################################################
#pathCosmoResults <- '/Users/battal/Cerens_files/fMRI/Processed/RhythmCateg/RhythmBlock/data/derivatives/cosmoMvpa/neurosynth-freesurfer-2beta-diffVoxelNb'
#pathCosmoResults <- '/Users/battal/Cerens_files/fMRI/Processed/RhythmCateg/RhythmBlock/data/derivatives/cosmoMvpa/contrast'

pathCosmoResults <- '/Users/battal/Dropbox/Work/CPPLab/Cerens_files_old/Result_sheets/'

########

mvpa <- NA

#dataNames <- paste(pathCosmoResults, '*.csv', sep ='/')
dataNames = "*20210524.csv"
temp = list.files(path = pathCosmoResults,pattern=dataNames)
csvFileNb <- length(temp)     

resultFiles = list()
for (i in 1:csvFileNb) {
  fileToRead = paste(pathCosmoResults, temp[i], sep ='/')
  x  = read.csv(fileToRead)
  x$FileID <- i
  resultFiles[i] = list(x)
}

# bind txt files using rbind comment
mvpa = do.call(rbind, resultFiles)

#####
### DECODING EB - SC PT-V5
#### order things a bit for easy manipulation for plotting

#make sure subjects are factor
mvpa$sub <-as.factor(mvpa$sub)

#make group order to call them in plot EB first
mvpa$group_order<-ifelse(mvpa$group == 'EB', 1, 2)
mvpa$group <-ifelse(mvpa$group == 'CONT', 'SC', 'EB')

#make roi order to call accordingly
mvpa$roi_order <- ifelse(mvpa$roi == 'lV5_6mm_2.nii', 1, 
                         ifelse(mvpa$roi == 'rV5_6mm_2.nii', 2, 
                                ifelse(mvpa$roi == 'lPT_6mm_mo.nii', 3,4)))

mvpa$roiName <- ifelse(mvpa$roi == 'lV5_6mm_2.nii', 'LeftV5', 
                         ifelse(mvpa$roi == 'rV5_6mm_2.nii', 'RightV5', 
                                ifelse(mvpa$roi == 'lPT_6mm_mo.nii', 'LeftPT','RightPT')))

#add motion or static condition
mvpa$isMotion <- ifelse(tolower(substring(mvpa$conditions, 1, 1)) == 's',0,1) 
mvpa$condition <- ifelse(tolower(substring(mvpa$conditions, 1, 1)) == 's','static','motion') 

# filter out some unnecassary columns
# filter out plane column 
mvpa[,6] <- NULL
mvpa[,6] <- NULL
mvpa[,7] <- NULL
mvpa[,7] <- NULL

# let's multiple accuracy with 100
mvpa$value <- mvpa$value * 100
mvpa.motion <- subset(mvpa, isMotion ==1)
mvpa.static <- subset(mvpa, isMotion == 0)

# take only 4 motion
mvpa.4MS <- subset(mvpa, conditions == 'Motion4' | conditions =='Static4')
# below does not work because it converts the second letter lower case, I want: SLeftvsRight
# mvpa.static$conditions <-str_to_title(mvpa.static$conditions)

# # take only 4 motion
# mvpa.4motion <-subset(mvpa.motion, conditions =='Motion4')
# mvpa.4static <-subset(mvpa.static, conditions =='Static4')


# summary stats
mvpa.4MS$condRoi <- paste(mvpa.4MS[,'condition'], mvpa.4MS[,'roiName'])

# small trial to only see one condition one roi across subjects
# mvpa.4M <- subset(mvpa.4MS, conditions == 'Motion4' & roiName == 'RightV5')
  
df <- summarySE(data = mvpa.4MS, 
                groupvars=c('group', 'group_order', 'condRoi','roiName', 'roi_order','condition'),
                measurevar='value')
df

setlimit = c(10,55) 
setbreak = c(10,20,30,40,50)

shapesize = 2
shapetype = 21
shapestroke = 1
transparent = 1 #0.6
jitter  = position_jitterdodge(0.2) # position_jitter(width=0.3)


fig <- ggplot(data = mvpa.4MS, 
              aes(x = reorder(condRoi,roi_order), 
                  y = value, 
                  color = group,
                  group = group)) +
  geom_point(data=mvpa.4MS,aes(x = reorder(condRoi,roi_order), y = value), size = shapesize,
             position = jitter, shape = shapetype, stroke = shapestroke) + 
  
  stat_summary(aes(color=group), fun=mean, fun.min = mean, fun.max = mean, geom="crossbar", size=0.6, width=0.6,position = position_dodge(width=.75)) +
  
  theme_classic() +
  geom_errorbar(data = df, 
                aes(ymin = value-se, ymax = value+se, group = reorder(group, group_order)), 
                color = 'black',size=0.5, width=0.15, alpha = transparent, position = position_dodge(width=.75)) +


  geom_hline(yintercept=c(25), linetype="dotted", colour="black", size=.5) +
  ggtitle("") +
  ylab("Decoding Accuracy (%)") +
  xlab("") +
  theme(axis.text.x=element_text(size=8, face = 'bold', angle=0, colour='black')) + # face = 'bold', 
  theme(axis.text.y=element_text(size=8, angle=0, colour='black')) +
  theme(axis.title.y=element_text(size=11, angle=90, colour='black')) +
  scale_y_continuous(limits=setlimit, breaks=setbreak, position="left") +
  scale_x_discrete(labels = c("Moving","Static","Moving","Static","Moving","Static","Moving","Static"))+
  # theme(text=element_text(family="Microsoft Sans Serif")) +
  scale_color_manual(values=c("purple","gray")) +
  theme(legend.position="none")
fig

filename <- paste(pathCosmoResults, 'Decoding_4Motion4Static_EBSC.png', sep = '')

# ggsave(filename, fig, dpi=300, width=8, height=2.4)

ggsave(filename, fig, dpi=300, width=4.5, height=3)

#### add fonts for ggpplot. 
filename <- paste(pathCosmoResults, 'Decoding_4Motion4Static_EBSC.pdf', sep = '')
ggsave(filename, fig, dpi=300, width=8, height=2.4)
