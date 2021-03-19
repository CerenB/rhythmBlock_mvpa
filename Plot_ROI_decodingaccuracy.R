library(ggplot2)
library(doBy)
library(cowplot)
library(Rmisc)




#######################################################
#pathCosmoResults <- '/Users/battal/Cerens_files/fMRI/Processed/RhythmCateg/RhythmBlock/derivatives/cosmoMvpa/before16032021'
pathCosmoResults <- '/Users/battal/Cerens_files/fMRI/Processed/RhythmCateg/RhythmBlock/derivatives/cosmoMvpa/'

########

mvpa <- NA

dataNames <- paste(pathCosmoResults, '*.csv', sep ='/')
dataNames = "*.csv"
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




#######

# audBG <- read.csv(paste(pathCosmoResults, 'RhythmBlockDecoding_freesurfer_s2_20210303.csv', sep ='/'))
# 
# motor <- read.csv(paste(pathCosmoResults, 'RhythmBlockDecoding_neurosnyth_s2_20210303.csv', sep ='/'))
# 
# BG <- read.csv(paste(pathCosmoResults, 'RhythmBlockDecoding_BG_s0_vx520_20210308.csv', sep ='/'))
#   
# mvpa <- BG
# 
# mvpa <- rbind(audBG, motor)

# mvpa <- mvpa[-c(7:9)]
mvpa$subID <-as.factor(mvpa$subID)

mvpa$roi_order <- ifelse(mvpa$mask == 'leftAud', 1, 
                         ifelse(mvpa$mask == 'rightAud', 2, 
                                ifelse(mvpa$mask == 'leftBG', 3, 
                                       ifelse(mvpa$mask == 'rightBG',4,
                                              ifelse(mvpa$mask == 'leftPremotor', 5,
                                                     ifelse(mvpa$mask == 'rightPremotor',6,7))))))
mvpa$roi_color_code <- ifelse(mvpa$mask == 'leftAud', 1, 
                         ifelse(mvpa$mask == 'rightAud', 1, 
                                ifelse(mvpa$mask == 'leftBG', 2, 
                                       ifelse(mvpa$mask == 'rightBG',2,
                                              ifelse(mvpa$mask == 'leftPremotor', 3,
                                                     ifelse(mvpa$mask == 'rightPremotor',3,4))))))

mvpa$mask <- ifelse(mvpa$mask == 'leftAud', 'audL', 
                         ifelse(mvpa$mask == 'rightAud', 'audR', 
                                ifelse(mvpa$mask == 'leftBG', 'bgL', 
                                       ifelse(mvpa$mask == 'rightBG','bgR',
                                              ifelse(mvpa$mask == 'leftPremotor', 'preL',

                                                                                                          ifelse(mvpa$mask == 'rightPremotor','preR','SMA'))))))
# sma roi consist of left/right hemispheres but for plotting 
# we can show under the same vx size as the others
# better way would be dividing the masks an re-run the decoding
mvpa$voxelToPlot <- mvpa$choosenVoxNb
mvpa$voxelToPlot<- ifelse(mvpa$mask == 'SMA', mvpa$voxelToPlot/2, mvpa$voxelToPlot)

# ==============================================================================
# summary stats
df <- summarySE(data = mvpa, 
                groupvars=c('mask', 'roi_order', 'image','ffxSmooth','voxelToPlot','roi_color_code'),
                measurevar='accuracy')


#################
pd <- position_dodge(0.1)

filename <- paste(pathCosmoResults, '/plot/', 'Decoding_SimpleVsComplex-AllCondition.png', sep = '')
title <- paste('Simple vs Complex Rhythm Decoding ')


fig <- ggplot(df, aes(x = reorder(mask, roi_order), y = accuracy)) + 
  geom_point(data = mvpa, aes(group=subID), pos=pd, size=2, color=grey(0.8)) + 
  geom_point(size=2,col='black') + 
  geom_hline(yintercept=c(.5), linetype="dotted", colour="red", size=.5) +
  facet_grid(vars(image,ffxSmooth), vars(voxelToPlot)) +
  geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),size=1,width=0.2) + 
  scale_y_continuous(limits=c(0, .90)) +
  xlab('ROIs')+
  ylab('classification acc.')+
  ggtitle(title)+
  theme_cowplot(font_size = 12,
                font_family = "",
                line_size = 0.3,
                rel_small = 6/12,
                rel_tiny = 9/12,
                rel_large = 14/12)


fig

ggsave(filename, device="png", units="in", width=18, height=9.08, dpi=300)  



########
# with for loop
# filterSmoothing <- c('0','2')
# filterImage <- c('beta', 't_maps')
# filterVoxelNb <- c('100','150','250','400','520')

# for (iImage in filterImage) {
#   
#   
#   for (iSmoothing in filterSmoothing) {
#     
#     for (iVoxelNb in filterVoxelNb) {
#       
#       print(paste(iImage, iSmoothing, iVoxelNb))
#       
#       title <- paste('SimplevsComplex - ',iImage,'s',iSmoothing, 'voxel', iVoxelNb)
#       
#       filename <- paste(pathCosmoResults, '/plot/', 'simpleVscomplex-', iImage, '_s', iSmoothing,'_', iVoxelNb,'_Vx.png', sep = '')
#       
#       fig <- ggplot(subset(df, image == iImage & choosenVoxNb == iVoxelNb & ffxSmooth == iSmoothing), aes(x = reorder(mask, roi_order), y = accuracy)) + 
#         geom_point(data = subset(mvpa, image == iImage & iVoxelNb == choosenVoxNb & iSmoothing == ffxSmooth), aes(group=subID), pos=pd, size=3, color=grey(0.8)) + 
#         geom_point(size=4,col='black') + 
#         geom_hline(yintercept=c(.5), linetype="dotted", colour="red", size=.5) +
#         geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),size=1,width=0.2) + 
#         scale_y_continuous(limits=c(0, .90)) +
#         xlab('ROIs')+
#         ylab('classification acc.')+
#         ggtitle(title)+
#         theme_classic() +
#         theme(
#           text=element_text(size=16),
#           legend.position = 'right',
#           legend.text=element_text(size=14),
#           axis.line = element_line(size = 0.6),
#           axis.text.x = element_text(size=14,colour="black"),
#           axis.text.y = element_text(size=14, colour='black'))
#       ggsave(filename, device="png", units="in", width=9, height=4.54, dpi=300)  
#     }
#   }
#   
# }



