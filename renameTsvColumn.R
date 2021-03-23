library(tidyverse)

# let's read tsv files and reorganize the trial_type so each repetition of 
# a condition would be labeled differently and thus we can model them in 
# GLM as separated regressors == separated betas

pathToFunc <- '/Users/battal/Cerens_files/fMRI/Processed/RhythmCateg/RhythmBlock/code/rhythmBlock_fMRI_analysis/lib/bids-R/bidsr_queryEvents.R'
source(pathToFunc)


bidsRoot <- '/Users/battal/Cerens_files/fMRI/Processed/RhythmCateg/RhythmBlock/derivatives/cpp_spm/sub-011' 
taskName <- 'RhythmBlock' 

taskEventsFiles <- bidsr_queryEvents(bidsRoot = bidsRoot, 
                                     taskName = taskName)
for (i in 1:length(taskEventsFiles)) {
  
  tsv <- read.table(paste(bidsRoot, taskEventsFiles[i], sep = '/'), header = TRUE)
  
  # if it is simple_block or complex_block, rewrite it with "simple_block_stepNum"
  tsv$trial_type <- ifelse(tsv$trial_type == 'block_simple' | tsv$trial_type == 'block_complex', 
                            paste(tsv$trial_type, tsv$stepNum, sep = '_'), tsv$trial_type)
  
  write.table(tsv,
              paste(bidsRoot, taskEventsFiles[i], sep = '/'),
              row.names = FALSE,
              sep = '\t',
              quote = FALSE)
  
}
