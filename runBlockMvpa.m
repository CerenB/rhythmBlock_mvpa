function accu = runBlockMvpa

%cd(fileparts(mfilename('fullpath')));
%% define paths
%spm - for now
warning('off');
addpath(genpath('/Users/battal/Documents/MATLAB/spm12'));
%cosmo
cosmo = '/Users/battal/Documents/MATLAB/CoSMoMVPA';
addpath(genpath(cosmo));
cosmo_warning('once')

%libsvm
libsvm = '/Users/battal/Documents/MATLAB/libsvm';
addpath(genpath(libsvm));
% verify it worked.
cosmo_check_external('libsvm'); % should not give an error

% add cpp-spm
cppSPM = '/Users/battal/Documents/GitHub/CPPLab/CPP_SPM';
addpath(genpath(fullfile(cppSPM, 'src')));
addpath(genpath(fullfile(cppSPM, 'lib')));

% get options
opt = getOptionBlockMvpa();

checkDependencies();

%define the 4D maps to be used
funcFWHM = 0;

%% roi path
% use parcels or NS mask?
useParcel = 1;


maskPath = fullfile(fileparts(mfilename('fullpath')), '..', ...
              '..','..', 'RhythmCateg_ROI','neurosynth',...
               'functional', 'derivatives');
           
% masks to decode           
maskName = {'leftrrthres_7premotor_FDR_0.01.nii', ...
           'rightrrthres_7premotor_FDR_0.01.nii',...
           'rrthres_10sma_FDR_0.01.nii'};
       
% use in output roi name       
maskLabel = {'rightPremotor', 'leftPremotor','SMA'};      

%use in output name
roiSource = 'neurosnyth';

% parcels
% if use parcels, re-writes mask names:
if useParcel == 1
    maskPath = fullfile(fileparts(mfilename('fullpath')), '..', ...
        '..','..', 'RhythmCateg_ROI','freesurfer');
    
    maskName = {'rrlauditorycx.nii', 'rrrauditorycx.nii', ...
                'rrlbasalganglia.nii',...
                'rrrbasalganglia.nii'};
%     maskName = {'thres1_s1_dec_rrlauditorycx.nii', ...
%         'thres1_s1_dec_rrrauditorycx.nii',...
        %         'rrlbasalganglia.nii',...
%         'rrrbasalganglia.nii'};

    maskLabel = {'leftAud','rightAud', 'leftBG','rightBG'};
    
    roiSource = 'Aud';
end

%% set output folder/name
savefileMat = fullfile(opt.pathOutput, ...
            [opt.taskName, ...
            'Decoding_', ...
            roiSource, ...
            '_s', num2str(funcFWHM), ...
            '_', datestr(now, 'yyyymmdd'), '.mat' ]);
        
savefileCsv = fullfile(opt.pathOutput, ...
            [ opt.taskName, ...
            'Decoding_', ...
            roiSource, ...
            '_s', num2str(funcFWHM), ...
            '_', datestr(now, 'yyyymmdd'), '.csv' ]);   

%% MVPA options
opt.mvpa.tool = 'cosmo';

opt.mvpa.normalization = 'zscore';

opt.mvpa.child_classifier = @cosmo_classify_libsvm;

opt.mvpa.feature_selector = @cosmo_anova_feature_selector;

% take the most responsive xx nb of voxels
opt.mvpa.ratioToKeep = [ 110 ];

% set which type of ffx results you want to use
opt.mvpa.map4D = {'beta', 't_maps'};

% design info
opt.mvpa.nbRun = 9;
opt.mvpa.nbTrialRepetition = 1;

% set cosmo mvpa structure
condLabelNb = [1 2];
condName= {'simple', 'complex'};
decodingCondition = 'simpleVscomplex';


%% let's get going!


% set structure array for keeping the results
accu = struct( ...
    'subID', [], ...
    'mask', [], ...
    'accuracy', [], ...
    'prediction', [], ...
    'maskVoxNb', [], ...
    'choosenVoxNb', [], ...
    'image', [], ...
    'ffxSmooth', [], ...
    'roiSource', [], ...
    'decodingCondition', [], ...
    'permutation', [], ...
    'imagePath', []);                

% get dataset info
[group, opt] = getData(opt);

count = 1;
for iGroup = 1:length(group)
    
    %groupName = group(iGroup).name;
    
    for iSub = 1:group(iGroup).numSub
        
        % get FFX path
        subID = group(iGroup).subNumber{iSub};
        ffxDir = getFFXdir(subID, funcFWHM, opt);
        
        % get subject folder name
        subFolder = ['sub-', opt.subjects{iSub}];
        
        for iImage = 1: length(opt.mvpa.map4D)
            
            for iMask = 1: length(maskName)
                
                %choose the mask
                mask = fullfile(maskPath, maskName{iMask});
                if useParcel == 1
                    mask = fullfile(maskPath, subFolder, maskName{iMask});
                end
                
                % 4D image
                imageName = ['4D_', opt.mvpa.map4D{iImage},'_', num2str(funcFWHM), '.nii'];
                image = fullfile(ffxDir, imageName);
                
                % load cosmo input
                ds = cosmo_fmri_dataset(image,'mask',mask);
                
                % Getting rid off zeros
                zeroMask=all(ds.samples==0,1);
                
                ds = cosmo_slice(ds, ~zeroMask, 2);
                
                %calculate the mask size
                maskVoxel = size(ds.samples,2)
                
                % set cosmo structure
                ds = setCosmoStructure(opt, ds, condLabelNb, condName);
                
                % slice the ds according to your targers (choose your
                % train-test conditions
                ds = cosmo_slice(ds, ds.sa.targets == 1 | ds.sa.targets == 2);
                
                % remove constant features
                ds = cosmo_remove_useless_data(ds);
                
                % partitioning, for test and training : cross validation
                partitions = cosmo_nfold_partitioner(ds);
                
                % define the voxel number for feature selection
                % set ratio to keep depending on the ROI dimension
                opt.mvpa.feature_selection_ratio_to_keep = opt.mvpa.ratioToKeep;
                
                % ROI mvpa analysis
                [pred, accuracy] = cosmo_crossvalidate(ds, ...
                    @cosmo_meta_feature_selection_classifier, ...
                    partitions, opt.mvpa);
                
                
                 %% PERMUTATION PART
                 %number of iterations
                 nbIter = 100;
                 
                 % allocate space for permuted accuracies
                 acc0 = zeros(nbIter,1); 
                 
                 % make a copy of the dataset
                 ds0 = ds; 
                 
                 %for _niter_ iterations, reshuffle the labels and compute accuracy
                 % Use the helper function cosmo_randomize_targets
                 for k=1:nbIter
                     ds0.sa.targets=cosmo_randomize_targets(ds);
                     [~, acc0(k)]=cosmo_crossvalidate(ds0, @cosmo_meta_feature_selection_classifier,...
                         partitions,opt.mvpa);
                 end
                 
                 p=sum(accuracy<acc0)/nbIter;
                 fprintf('%d permutations: accuracy=%.3f, p=%.4f\n', nbIter, accuracy, p);
                
                %% store output
                accu(count).subID = subID;
                accu(count).mask = maskLabel{iMask};
                accu(count).maskVoxNb = maskVoxel;
                accu(count).choosenVoxNb = opt.mvpa.feature_selection_ratio_to_keep;
                accu(count).image = opt.mvpa.map4D{iImage};
                accu(count).ffxSmooth = funcFWHM;
                accu(count).accuracy = accuracy;
                accu(count).prediction = pred;
                accu(count).imagePath = image; 
                accu(count).roiSource = roiSource;
                accu(count).decodingCondition = decodingCondition;
                accu(count).permutation = acc0';
                
                
                count = count + 1;
                
                fprintf(['Sub'  subID ' - area: ' maskLabel{iMask} ', accuracy: ' num2str(accuracy) '\n\n\n']);
                
            end
        end        
    end
     %% save output
     
     % mat file
     save(savefileMat, 'accu');
     
     % csv but with important info for plotting
     csvAccu = rmfield(accu,'permutation');
     csvAccu = rmfield(csvAccu, 'prediction');
     csvAccu = rmfield(csvAccu, 'imagePath');
     writetable(struct2table(csvAccu), savefileCsv)
end

end


function ds = setCosmoStructure(opt, ds, condLabelNb, condName)
% sets up the target, chunk, labels by stimuli condition labels, runs,
% number labels.

%chunk (runs), target (condition), labels (condition names)
trialsPerRun = length(condLabelNb) * opt.mvpa.nbTrialRepetition;

chunks = repmat([1:opt.mvpa.nbRun]',1,opt.mvpa.nbTrialRepetition)';
chunks = chunks(:);
chunks = repmat(chunks,trialsPerRun,1);

targets = repmat(condLabelNb',1,opt.mvpa.nbRun)';
targets = targets(:);

condName = repmat(condName',1,opt.mvpa.nbRun)';
condName = condName(:);

ds.sa.targets = targets;
ds.sa.chunks = chunks;
ds.sa.labels = condName;

%figure; imagesc(ds.sa.targets);

end

