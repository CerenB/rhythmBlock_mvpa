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
  run ../../rhythmBlock_fMRI_analysis/lib/CPP_BIDS_SPM_pipeline/initCppSpm.m;
  
  % roi path
  parcelPath = fullfile(fileparts(mfilename('fullpath')), '..', ...
                        '..', '..', '..','RhythmCateg_ROI', 'freesurfer');
     
  % load your options
  opt = getOptionBlockMvpa();

  %% make Freesurfer based ROIs
  % which parcels to use?
  useAudParcel = 0;
  
  % make rois
  % action 1 = makes individual/separate parcels and create binary masks
  % action 2 = combines the left hemisphere masks into 1 (e.g. left-basal
  % ganglia, or left-auditory cortex).
  % action 3 = reslice the masks by using spm func to the functional (4D.nii) image
  % action 4 = smooths the masks for auditory parcels
  action = 3;
  info = parcel2mask(action, parcelPath, opt, useAudParcel);
  
  
  useAudParcel = 1;
  action = 4;
  info = parcel2mask(action, parcelPath, opt, useAudParcel);
  %% make HMAT based ROIs
  
  
  
  %% run mvpa 
  
  % use parcels or NS masks?
  opt.mvpa.useParcel = 1;
  accuracy = calculateMvpa(opt);
  
  
  
  
  
  
  
  
  