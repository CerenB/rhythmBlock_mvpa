clear ;
clc;

%% set paths
  % spm
  warning('off');
  addpath(genpath('/Users/battal/Documents/MATLAB/spm12'));
  
  % add cpp repo
%   run ../../rhythmBlock_fMRI_analysis/lib/CPP_BIDS_SPM_pipeline/initCppSpm.m;
  run /Users/battal/Documents/GitHub/CPPLab/CPP_SPM/initCppSpm.m; 
  
  % load your options
  opt = getOptionBlockMvpa();

  %% make Freesurfer based ROIs
  
    % roi path
  parcelPath = fullfile(fileparts(mfilename('fullpath')), '..', ...
                        '..', '..', '..','RhythmCateg_ROI', 'freesurfer');
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
  action = 3;
  info = parcel2mask(action, parcelPath, opt, useAudParcel);

  action = 4;
  opt.maskFWHM = 2;
  info = parcel2mask(action, parcelPath, opt, useAudParcel);
  
  %% make HMAT based ROIs
  opt.dir.roi = fullfile(fileparts(mfilename('fullpath')),  ...
      '..', '..', '..', '..', 'RhythmCateg_ROI', 'hmat');
  
  % dont start from copying the raw masks
  opt.copyRawToDerivatives.do = true;
  
  % reslice the mask into mni  func space
  opt.reslice.do = true;
  
  % reslice mask to normalised func space
  % (only once is enough for MNI space roi)
  opt.resliceFunc.do = false;
  
  % transform space from "mni func image" to "individual func image"
  opt.inversTransform.do = true;
  
  % count the mask voxel
  opt.countVoxel.do = true;
  
  opt = hmat2mask(opt);
  

  %% make a ROI from JulichBrain Cytotechtonic Atlas
  directory = '/Users/battal/Documents/MATLAB/spm12/toolbox/Anatomy_julichbrain_v29-pmaps-4d';
  atlas = fullfile(directory, 'JULICH_BRAIN_CYTOARCHITECTONIC_MAPS_2_9_MNI152_2009C_NONL_ASYM.pmaps.nii');
                

  % get nifti header
  header = spm_vol(atlas);

  % get atlas content
  volume = spm_read_vols(header(207));
  size(volume)
  
  % create binary mask of map
  binary = volume;
  binary(binary >0.0) = 1.0;
  

  % prepare header for the output
  binary_header = header;
  binary_header.fname = 'mask.nii';
 
  % write mask
  spm_write_vol(binary_header, binary);
  
 %%%% IT DOES NOT WORK -----  
  
  
 
 % ending work session?
 run /Users/battal/Documents/GitHub/CPPLab/CPP_SPM/uninitCppSpm.m;
  
  
  
  