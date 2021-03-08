% (C) Copyright 2019 CPP BIDS SPM-pipeline developpers

function opt = getOptionBlockMvpa()
  % opt = getOption()
  % returns a structure that contains the options chosen by the user to run
  % slice timing correction, pre-processing, FFX, RFX.

  if nargin < 1
    opt = [];
  end

  % group of subjects to analyze
  opt.groups = {''};
  % suject to run in each group
  opt.subjects = {'001', '002', '003', '004', '005', '006','007',...
                  '008', '009', '010','011'};

  % Uncomment the lines below to run preprocessing
  % - don't use realign and unwarp
  opt.realign.useUnwarp = true;

  % we stay in native space (that of the T1)
  opt.space = 'individual'; % 'individual', 'MNI'

  % The directory where the data are located
  opt.dataDir = fullfile(fileparts(mfilename('fullpath')), ...
                         '..', '..',  'raw');
  opt.derivativesDir = fullfile(opt.dataDir, '..', 'derivatives', 'cpp_spm');
  
  opt.pathOutput = fullfile(opt.dataDir, '..', 'derivatives','cosmoMvpa');

  % task to analyze
  opt.taskName = 'RhythmBlock';

  opt.parallelize.do = false;
  opt.parallelize.nbWorkers = 1;
  opt.parallelize.killOnExit = true;



  %% DO NOT TOUCH
  opt = checkOptions(opt);
  saveOptions(opt);


end
