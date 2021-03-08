library(ggplot2)
library(doBy)

pathCosmoResults <- '/Users/battal/Cerens_files/fMRI/Processed/RhythmCateg/RhythmBlock/derivatives/cosmoMvpa'


audBG <- read.csv(paste(pathCosmoResults, 'RhythmBlockDecoding_freesurfer_s2_20210303.csv', sep ='/'))

motor <- read.csv(paste(pathCosmoResults, 'RhythmBlockDecoding_neurosnyth_s2_20210303.csv', sep ='/'))

mvpa <- rbind(audBG, motor)

mvpa <- mvpa[-c(7:9)]
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
smoothing <- '2'
filterImage <- c('beta', 't_maps')

# ==============================================================================

    
for (iImage in filterImage) {
  

  title <- paste('Simple vs Complex Rhythm - ',iImage,'s',smoothing)
  
  filename <- paste(pathCosmoResults, '/plot/', 'simpleVscomplex-', iImage, '_s', smoothing,'.png', sep = '')
  
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




##################
pd <- position_dodge(0.1)
df <- summarySEwithin(all_summary, 
                      withinvars='rhythmID',
                      idvar='ID',
                      measurevar='zscore.meter.eeg')
ggplot(df, aes(rhythmID, zscore.meter.eeg)) + 
  geom_point(data=all_summary, aes(group=ID), pos=pd, size=3, color=grey(0.8)) + 
  geom_point(size=4,col='black') + 
  geom_errorbar(aes(ymin=zscore.meter.eeg-ci,ymax=zscore.meter.eeg+ci),size=1,width=0.2) + 
  theme_cowplot()

