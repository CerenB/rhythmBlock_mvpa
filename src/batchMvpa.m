clear ;
clc;

%% set paths
  % spm
  warning('off');
  addpath(genpath('/Users/battal/Documents/MATLAB/spm12'));
  % cosmo
  cosmo = '~/Documents/MATLAB/CoSMoMVPA';
  addpath(genpath(cosmo));
  cosmo_warning('once');

  % libsvm
  libsvm = '~/Documents/MATLAB/libsvm';
  addpath(genpath(libsvm));
  % verify it worked.
  cosmo_check_external('libsvm'); % should not give an error
  
  % add cpp repo
    % add cpp repo
  run ../lib/bidspm/initCppSpm.m;
%   run ../../rhythmBlock_fMRI_analysis/lib/CPP_BIDS_SPM_pipeline/initCppSpm.m;
  
  % add mini-helper functions
  addpath(genpath(fullfile(pwd, 'subfun')));
  
  % load your options
  opt = getOptionBlockMvpa();

  %% run mvpa 
  
  % tapping info
  opt.tapper.do = 1;
  opt.tapper.good = [1,2,5,7,8,9,10,12,14,15,16,21,25,27,28,30,32,33];

  % use parcels or NS masks?
  roiSource = 'contrastSTGOnly'; % 'freesurfer', 'neurosynth', 'hmat'
  
  % want to omit something? 
  opt.omit.subject.do = 1;
  opt.omit.subject.mask = 'cerebellum';
  opt.omit.subject.ID = [3, 7, 26,32];
  
  accuracy = calculateMvpa(opt, roiSource);
  
  %% run mvpa on auditory cortex - Juelich atlas
  opt.tapper.do = 1;
  opt.tapper.good = [1,2,5,7,8,9,10,12,14,15,16,21,25,27,28,30,32,33];

  % use parcels or NS masks?
  roiSource = 'jubrainatlas'; % 'freesurfer', 'neurosynth', 'hmat'
  opt.probabThreshold = 10;
  opt.omit.subject.do = 0;

  accuracy = calculateMvpa(opt, roiSource);
  
  
  