library(tidyverse)

pathToFunc <- '/Users/battal/Cerens_files/fMRI/Processed/RhythmCateg/RhythmBlock/code/rhythmBlock_fMRI_analysis/lib/bids-R'
source(pathToFunc)

bidsRoot <- '/Users/barilari/Desktop/data_temp/ses-004_cpp_spm/events_bu' 
taskName <- 'rdkMohBimodalMotion' 

taskEventsFiles <- bidsr_queryEvents(bidsRoot = bidsRoot, 
                                     taskName = taskName)
for (i in 1:length(taskEventsFiles)) {
  
  temp <- read.table(paste(bidsRoot, taskEventsFiles[i], sep = '/'), header = TRUE)
  
  temp$trial_type <- ifelse(temp$direction == '0.000000' | temp$direction == '0', 
                            paste(temp$trial_type, 'right', sep = '_'),
                            ifelse(temp$direction == '90.000000' | temp$direction == '90', 
                                   paste(temp$trial_type, 'up', sep = '_'),
                                   ifelse(temp$direction == '180.000000' | temp$direction == '180', 
                                          paste(temp$trial_type, 'left', sep = '_'), 
                                          ifelse(temp$direction == '270.000000' | temp$direction == '270', 
                                                 paste(temp$trial_type, 'down', sep = '_'), 
                                                 ifelse(temp$trial_type == 'response', 
                                                        'response', 
                                                        'trigger')))))
  
  write.table(temp,
              paste(bidsRoot, taskEventsFiles[i], sep = '/'),
              row.names = FALSE,
              sep = '\t',
              quote = FALSE)
  
}
