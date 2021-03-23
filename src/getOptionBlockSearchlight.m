% (C) Copyright 2019 CPP BIDS SPM-pipeline developpers

function opt = getOptionBlockSearchlight()
  % returns a structure that contains the options chosen by the user to run
  % searchlight.

  if nargin < 1
    opt = [];
  end

  % group of subjects to analyze
  opt.groups = {''};
  % suject to run in each group
  opt.subjects = {'001'};

  %               '001', '002', '003', '004', '005', '006', '007', ...
  %                   '008', '009', '010', '011'

  % Uncomment the lines below to run preprocessing
  % - don't use realign and unwarp
  opt.realign.useUnwarp = true;

  % we stay in native space (that of the T1)
  opt.space = 'MNI'; % 'individual', 'MNI'

  % The directory where the data are located
  opt.dataDir = fullfile(fileparts(mfilename('fullpath')), ...
                         '..', '..', '..', 'raw');
  opt.derivativesDir = fullfile(opt.dataDir, '..', 'derivatives', 'cpp_spm');

  opt.pathOutput = fullfile(opt.dataDir, '..', 'derivatives', 'cosmoMvpa');
  
  % multivariate
  opt.model.file = fullfile(fileparts(mfilename('fullpath')), '..', ...
                            'model', 'model-RhythmBlockDecoding1_smdl.json');

  % task to analyze
  opt.taskName = 'RhythmBlock';

  opt.parallelize.do = true;
  opt.parallelize.nbWorkers = 4;
  opt.parallelize.killOnExit = true;

  %% DO NOT TOUCH
  opt = checkOptions(opt);
  saveOptions(opt);
  % we cannot save opt with opt.mvpa, it crashes

  %% mvpa options

  % define the 4D maps to be used
  opt.funcFWHM = 2;

  % Define a neighborhood with approximately 100 voxels in each searchlight.
  opt.mvpa.searchlightVoxelNb = 100; % 100 150

  % set which type of ffx results you want to use
  opt.mvpa.map4D = {'beta', 't_maps'};

  % design info
  opt.mvpa.nbRun = 9;
  opt.mvpa.nbTrialRepetition = 1;

  % cosmo options
  opt.mvpa.tool = 'cosmo';

  % Use the cosmo_cross_validation_measure and set its parameters
  % (classifier and partitions) in a measure_args struct.
  opt.mvpa.measure = @cosmo_crossvalidation_measure;

  % Define which classifier to use, using a function handle.
  % Alternatives are @cosmo_classify_{svm,matlabsvm,libsvm,nn,naive_bayes}
  opt.mvpa.classifier = @cosmo_classify_libsvm;

end
