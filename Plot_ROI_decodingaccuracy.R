library(ggplot2)
library(doBy)
library(cowplot)
library(Rmisc)




#######################################################
pathCosmoResults <- '/Users/battal/Cerens_files/fMRI/Processed/RhythmCateg/RhythmBlock/derivatives/cosmoMvpa'

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
mvpa$mask <- ifelse(mvpa$mask == 'leftAud', 'L Aud', 
                         ifelse(mvpa$mask == 'rightAud', 'R Aud', 
                                ifelse(mvpa$mask == 'leftBG', 'L BG', 
                                       ifelse(mvpa$mask == 'rightBG','R BG',
                                              ifelse(mvpa$mask == 'leftPremotor', 'L Pre',
                                                     ifelse(mvpa$mask == 'rightPremotor','R Pre','SMA'))))))
filterSmoothing <- c('0','2')
filterImage <- c('beta', 't_maps')
filterVoxelNb <- c('100','150','250','400','520')

pd <- position_dodge(0.1)
# ==============================================================================
# summary stats
df <- summarySE(data = mvpa, 
                groupvars=c('mask', 'roi_order', 'image','ffxSmooth','choosenVoxNb'),
                measurevar='accuracy')


for (iImage in filterImage) {
  
  
  for (iSmoothing in filterSmoothing) {
    
    for (iVoxelNb in filterVoxelNb) {
      
      print(paste(iImage, iSmoothing, iVoxelNb))
      
      title <- paste('SimplevsComplex - ',iImage,'s',iSmoothing, 'voxel', iVoxelNb)
      
      filename <- paste(pathCosmoResults, '/plot/', 'simpleVscomplex-', iImage, '_s', iSmoothing,'_', iVoxelNb,'_Vx.png', sep = '')
      
      fig <- ggplot(subset(df, image == iImage & choosenVoxNb == iVoxelNb & ffxSmooth == iSmoothing), aes(x = reorder(mask, roi_order), y = accuracy)) + 
        geom_point(data = subset(mvpa, image == iImage & iVoxelNb == choosenVoxNb & iSmoothing == ffxSmooth), aes(group=subID), pos=pd, size=3, color=grey(0.8)) + 
        geom_point(size=4,col='black') + 
        geom_hline(yintercept=c(.5), linetype="dotted", colour="red", size=.5) +
        geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),size=1,width=0.2) + 
        scale_y_continuous(limits=c(0, .90)) +
        xlab('ROIs')+
        ylab('classification acc.')+
        ggtitle(title)+
        theme_classic() +
        theme(
          text=element_text(size=16),
          legend.position = 'right',
          legend.text=element_text(size=14),
          axis.line = element_line(size = 0.6),
          axis.text.x = element_text(size=14,colour="black"),
          axis.text.y = element_text(size=14, colour='black'))
      ggsave(filename, device="png", units="in", width=9, height=4.54, dpi=300)  
    }
  }
  
}




##################


ggplot(df, aes(mask, accuracy)) + 
  geom_point(data=mvpa, aes(group=subID), pos=pd, size=3, color=grey(0.8)) + 
  geom_point(size=4,col='black') + 
  geom_errorbar(aes(ymin=accuracy-ci,ymax=accuracy+ci),size=1,width=0.2) + 
  theme_cowplot()


#################

for (iImage in filterImage) {
  
  
  title <- paste('Simple vs Complex Rhythm - ',iImage,'s',smoothing)
  
  filename <- paste(pathCosmoResults, '/plot/', 'simpleVscomplex-BG-', iImage, '_s', smoothing,'.png', sep = '')
  
  fig <- ggplot(data = subset(mvpa, image == iImage), aes(x = reorder(mask, roi_order), y = accuracy, colour = subID)) +
    theme_classic() +
    geom_point(position=position_jitter(width=0.12,h=0), size = 2, stroke= 0.6) +
    geom_hline(yintercept=c(.5), linetype="dotted", colour="red", size=.5) +
    stat_summary(fun=mean,colour = "black", geom="crossbar", size=0.5, width =0.4) +
    
    scale_y_continuous(limits=c(0, .90)) +
    xlab('ROIs')+
    ylab('classification acc.')+
    ggtitle(title)+
    theme(
      text=element_text(size=16),
      legend.position = 'right',
      legend.text=element_text(size=14),
      axis.line = element_line(size = 0.6),
      axis.text.x = element_text(size=14,colour="black"),
      axis.text.y = element_text(size=14, colour='black'))
  
  ggsave(filename, device="png", units="in", width=9, height=4.54, dpi=300)  
  
}

