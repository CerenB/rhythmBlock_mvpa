% (C) Copyright 2019 CPP BIDS SPM-pipeline developpers

function opt = getOptionBlockMvpa()
  % returns a structure that contains the options chosen by the user to run
  % decoding with cosmo-mvpa.

  if nargin < 1
    opt = [];
  end

  % suject to run in each group
%   opt.subjects = {'001'};
  opt.subjects = {'001', '002', '003', '004', '005', '006',...
                  '007', '008', '009', '010', '011', ...
                  '012', '013', '014', '015', '016', '017', ...
                  '018', '019', '020', '021', '023', ...
                  '024', '025', '026','027', '028', ...
                  '029', '030', '031', '032', '033'};
              
%   opt.subjects = {'012', '013', '014','015', '016', '017', ...
%                   '018', '019', '020', '021', '023'};

  opt = checkOptions(opt);  
  % Uncomment the lines below to run preprocessing
  % - don't use realign and unwarp
  opt.realign.useUnwarp = true;

  % we stay in native space (that of the T1)
  opt.space = 'MNI'; % 'individual', 'MNI'

  % The directory where the data are located
  opt.dataDir = fullfile(fileparts(mfilename('fullpath')), ...
                         '..', '..', '..', 'rhythmBlock_raw');
  opt.derivativesDir = fullfile(opt.dataDir, '..', 'rhythmBlock_derivatives_cpp_spm');
  
  opt.dir.stats = fullfile(opt.dataDir, '..', ...
                                  'derivatives','cpp_spm-stats');                          
  opt.dir.raw = opt.dataDir;

  opt.pathOutput = fullfile(opt.dataDir, '..', 'rhythmBlock_derivatives_cosmoMvpa');

  % multivariate
  opt.model.file = fullfile(fileparts(mfilename('fullpath')), '..', ...
                            'model', 'model-RhythmBlockDecoding1_smdl.json');

  % task to analyze
  opt.taskName = 'RhythmBlock';


  %% DO NOT TOUCH
%   opt = checkOptions(opt);
%   saveOptions(opt);
  % we cannot save opt with opt.mvpa, it crashes

  %% mvpa options
  
  % set cosmo mvpa structure
  opt.mvpa.condLabelNb = [1 2];
  opt.mvpa.condLabelName = {'simple', 'complex'};
  opt.mvpa.decodingCondition = 'simpleVscomplex';
  
  % define the 4D maps to be used
  opt.funcFWHM = 2;

  % take the most responsive xx nb of voxels
  opt.mvpa.ratioToKeep = 150; % 100 150 250 350 420

  % set which type of ffx results you want to use
  opt.mvpa.map4D = {'beta', 't_maps'};

  % design info
  opt.mvpa.nbRun = 9;
  opt.mvpa.nbTrialRepetition = 1;

  % cosmo options
  opt.mvpa.tool = 'cosmo';
  % opt.mvpa.normalization = 'zscore';
  opt.mvpa.child_classifier = @cosmo_classify_libsvm;
  opt.mvpa.feature_selector = @cosmo_anova_feature_selector;

  % permute the accuracies ?
  opt.mvpa.permutate = 1;

  % want to still run mvpa although the mask is smaller than desired
  % vx number? 
  opt.mvpa.useMaskVoxelNumber = 1;
  
end
