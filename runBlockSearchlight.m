function info = runBlockSearchlight

    % cd(fileparts(mfilename('fullpath')));
    %% define paths
    % spm - for now
    warning('off');
    addpath(genpath('/Users/battal/Documents/MATLAB/spm12'));
    % cosmo
    cosmo = '/Users/battal/Documents/MATLAB/CoSMoMVPA';
    addpath(genpath(cosmo));
    cosmo_warning('once');

    % libsvm
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

    % get the smoothing parameter for 4D map
    funcFWHM = opt.funcFWHM;

    %% define the mask
    % use in output name
    roiSource = 'wholeBrain';

    % mask name to decode
    maskName = {'mask.nii'};

    %% set output folder/name
    savefileMat = fullfile(opt.pathOutput, ...
                           [opt.taskName, ...
                                'Decoding_', ...
                                roiSource, ...
                                '_s', num2str(funcFWHM), ...
                                '_vx', num2str(opt.mvpa.searchlightVoxelNb), ...
                                '_', datestr(now, 'yyyymmdd'), '.mat']);

    %% MVPA options

    % set cosmo mvpa structure
    condLabelNb = [1 2];
    condName = {'simple', 'complex'};
    decodingCondition = 'simpleVscomplex';

    %% let's get going!

    % set structure array for keeping the results
    info = struct( ...
                  'subID', [], ...
                  'maskPath', [], ...
                  'maskVoxNb', [], ...
                  'searchlightVoxelNb', [], ...
                  'image', [], ...
                  'ffxSmooth', [], ...
                  'roiSource', [], ...
                  'decodingCondition', [], ...
                  'imagePath', []);

    % get dataset info
    [group, opt] = getData(opt);

    count = 1;
    for iGroup = 1:length(group)

        % groupName = group(iGroup).name;

        for iSub = 1:group(iGroup).numSub

            % get FFX path
            subID = group(iGroup).subNumber{iSub};
            ffxDir = getFFXdir(subID, funcFWHM, opt);

            %         % get subject folder name
            %         subFolder = ['sub-', opt.subjects{iSub}];

            for iImage = 1:length(opt.mvpa.map4D)

                % 4D image
                imageName = ['4D_', opt.mvpa.map4D{iImage}, '_', num2str(funcFWHM), '.nii'];
                image = fullfile(ffxDir, imageName);

                % mask
                mask = fullfile(ffxDir, maskName);

                % load cosmo input
                ds = cosmo_fmri_dataset(image, 'mask', mask);

                % Getting rid off zeros
                zeroMask = all(ds.samples == 0, 1);

                ds = cosmo_slice(ds, ~zeroMask, 2);

                % calculate the mask size
                maskVoxel = size(ds.samples, 2);

                % set cosmo structure
                ds = setCosmoStructure(opt, ds, condLabelNb, condName);

                % slice the ds according to your targers (choose your
                % train-test conditions
                ds = cosmo_slice(ds, ds.sa.targets == 1 | ds.sa.targets == 2);

                % remove constant features
                ds = cosmo_remove_useless_data(ds);

                % partitioning, for test and training : cross validation
                % can be different e.g. cosmo_oddeven_partitioner(ds_per_run)
                opt.mvpa.partitions = cosmo_nfold_partitioner(ds);

                % define a neightborhood
                nbrhood = cosmo_spherical_neighborhood(ds, 'count', ...
                                                       opt.mvpa.searchlightVoxelNb);

                % Run the searchlight
                svm_results = cosmo_searchlight(ds, nbrhood, opt.mvpa.measure, opt.mvpa);

                % store the relevant info
                info(count).subID = subID;
                info(count).maskPath = maskLabel{iMask};
                info(count).maskVoxNb = maskVoxel;
                info(count).searchlightVoxelNb = opt.mvpa.searchlightVoxelNb;
                info(count).image = opt.mvpa.map4D{iImage};
                info(count).ffxSmooth = funcFWHM;
                info(count).roiSource = roiSource;
                info(count).imagePath = image;
                info(count).decodingCondition = decodingCondition;

                count = count + 1;

                % Store results to disc
                savingResultFile = fullfile(opt.pathOutput, ...
                                            [opt.taskName, ...
                                                 'searchlight', ...
                                                 roiSource, ...
                                                 '_s', num2str(funcFWHM), ...
                                                 '_vx', num2str(opt.mvpa.searchlightVoxelNb), ...
                                                 '_', datestr(now, 'yyyymmdd'), 'nii']);
                cosmo_map2fmri(svm_results, savingResultFile);

                fprintf(['Sub'  subID ' is  being processed ....\n\n\n']);

            end
        end
        %% save output
        % mat file
        save(savefileMat, 'accu');

    end

end

function ds = setCosmoStructure(opt, ds, condLabelNb, condName)
    % sets up the target, chunk, labels by stimuli condition labels, runs,
    % number labels.

    % chunk (runs), target (condition), labels (condition names)
    trialsPerRun = length(condLabelNb) * opt.mvpa.nbTrialRepetition;

    chunks = repmat([1:opt.mvpa.nbRun]', 1, opt.mvpa.nbTrialRepetition)';
    chunks = chunks(:);
    chunks = repmat(chunks, trialsPerRun, 1);

    targets = repmat(condLabelNb', 1, opt.mvpa.nbRun)';
    targets = targets(:);

    condName = repmat(condName', 1, opt.mvpa.nbRun)';
    condName = condName(:);

    ds.sa.targets = targets;
    ds.sa.chunks = chunks;
    ds.sa.labels = condName;

    % figure; imagesc(ds.sa.targets);

end
