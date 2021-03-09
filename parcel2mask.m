function info = parcel2mask(action)

% it takes the atlas defined parcels from FS output (in the subject's
% anatomical native space)
% selectively gets the only parcel-of-interests
% then convert these parcels into subject's native functional resolution
% and lastly binarise these masks to be used in further analysis e.g.
% decoding/FFT/ ...

% action 1 = makes individual/separate parcels and create binary masks
% action 2 = combines the left hemisphere masks into 1 (e.g. left-basal
% ganglia, or left-auditory cortex).
% action 3 = reslice the masks by using spm func to the functional (4D.nii) image
% action 4 = smooths the masks 

% be careful, it has two options, with auditory cortex and with basal
% ganglia (parcel numbers, labels, roi names differ)

%% set paths
% spm
warning('off');
addpath(genpath('/Users/battal/Documents/MATLAB/spm12'));

% add cpp-spm
cppSPM = '/Users/battal/Documents/GitHub/CPPLab/CPP_SPM';
addpath(genpath(fullfile(cppSPM, 'src')));
addpath(genpath(fullfile(cppSPM, 'lib')));

% roi path
parcelPath = fullfile(fileparts(mfilename('fullpath')), '..', ...
    '..','..', 'RhythmCateg_ROI','freesurfer');

%% set info
% get subject options
opt = getOptionBlockMvpa();

% smooth info for ffxDir
funcFWHM = 0;

% which parcels to use?
useAudParcel = 1;

if useAudParcel == 1
    % parcel codes in FS LookUp Table
    parcelCodes =[11136, 11175, 11133, 12136,12175, 12133];
    
    % ROIs labels in FS LookUp table
    parcelLabels ={'lPT' , 'ltransverse', 'lsup_transv',...
        'rPT', 'rtransverse', 'rsup_transv'};
    
    % concat parcel name
    concatParcelName = 'auditorycx.nii';
    
    % choose which masks to realign with functional image
    maskToAlign = {'lauditorycx.nii','rauditorycx.nii'};
    
else
    parcelCodes = [11, 12, 13, 26, 50, 51, 52, 58];
    
    parcelLabels ={'lcaudate' , 'lputamen', 'lpallidum', 'lna',...
        'rcaudate', 'rputamen', 'rpallidum', 'rna'};

    concatParcelName = 'basalganglia.nii';

    maskToAlign = {'lbasalganglia.nii','rbasalganglia.nii'};
    
end

% parcel numbers in 1 hemisphere
parcelNb = length(parcelLabels)/2;

wholeParcelName = 'destrieux.nii';

%% smooth the realigned/resliced mask?
% only applied to AudCx, because BG does not need smoothing

% if yes, then define the masks to be smoothed, and FWHM in mm
maskToSmooth = {'rlauditorycx.nii','rrauditorycx.nii'};
maskFWHM = 1;
prefixDec = 'dec_';
prefixSmooth = ['s',num2str(maskFWHM),'_'];
threshold = 0.05;

%% let's shave fun
switch action
    
    case 1
        
        count = 1;
        for iSub = 1:length(opt.subjects)
            
            % get subject folder name
            subID = ['sub-', opt.subjects{iSub}];
            
            % loop for separate parcels
            for iParcel = 1:length(parcelCodes)
                
                % load the whole parcel
                parcel1 = load_nii(fullfile(parcelPath, subID, wholeParcelName));
                
                % assign zero to all the other parcels
                parcel1.img (parcel1.img ~= parcelCodes(iParcel)) = 0;
                
                %assign one to the selected parcel
                parcel1.img (parcel1.img == parcelCodes(iParcel)) = 1;
                
                % count the voxel number
                voxelNb = sum(parcel1.img(parcel1.img == 1));
                
                %save the new binary parcel/mask
                parcelOfInterest = fullfile(parcelPath, subID, ...
                    [parcelLabels{iParcel},'.nii']);
                save_nii(parcel1,parcelOfInterest);
                
                info(count).SUBname = subID;
                info(count).ROIname = parcelLabels{iParcel};
                info(count).ROIcode = parcelCodes(iParcel);
                info(count).ROInumVox = voxelNb;
                count = count +1;
            end
        end
        
        save(fullfile(parcelPath,'info'),'info');
        
        
    case 2
        % now let's combine the masks in 1 hemisphere
        count = 1;
        for iSub = 1:length(opt.subjects)
            
            % get subject folder name
            subID = ['sub-', opt.subjects{iSub}];
            
            
            % loop for separate parcels
            for iParcel = 1:parcelNb:(length(parcelCodes))
                
                parcelImg = [];
                
                % load the parcels in 1 hemisphere
                for i = iParcel:iParcel+(parcelNb-1)
                    
                    parcelName = [parcelLabels{i},'.nii'];
                    parcel = load_nii(fullfile(parcelPath, subID, parcelName));
                    parcelImg = cat(4,parcelImg,parcel.img);
                    
                end
                
                %open a template to save these parcels into 1 image
                temp = parcel ;
                temp.fileprefix = 'parcel';
                temp.img = [];
                
                % sum all the parcels
                temp.img = squeeze(sum(parcelImg,4));
                % make binary mask
                temp.img(temp.img == parcelNb)= 1;
                
                %save the new binary parcel/mask
                parcelOfInterest = fullfile(parcelPath, subID, ...
                    [parcelName(1),concatParcelName]);
                save_nii(temp,parcelOfInterest);
                
                %calculate the voxel num
                voxelNb = sum(temp.img(temp.img==1));
                
                % save voxel info into a struct
                info(count).SUBname = subID;
                info(count).ROIname = [parcelName(1),concatParcelName];
                info(count).ROInumVox = voxelNb;
                count = count +1;
            end
        end
        
        
    case 3
        % let's realign and reslice the masks with functional image
        % assuming the binary masks are ready to go
        
        count = 1;
        for iSub = 1:length(opt.subjects)
            
            % get subject folder name
            %subID = opt.subjects{iSub};
            subID = ['sub-', opt.subjects{iSub}];
            
            
            % get ffx - 4D image folder
            ffxDir = getFFXdir(opt.subjects{iSub}, funcFWHM, opt);
            
            for iMask = 1: length(maskToAlign)
                
                %choose the mask to realign and reslice
                maskName = maskToAlign{iMask};
                mask = fullfile(parcelPath, subID,maskName);
                image = fullfile(ffxDir,'4D_beta_0.nii');
                
                %% realign and reslice the new-roi
                % so that it is in the same space as your 4D images
                prefix = 'r';
                
                matlabbatch =[];
                matlabbatch{1}.spm.spatial.realign.write.data = {
                    [image,',1']
                    [mask,',1']
                    };
                matlabbatch{1}.spm.spatial.realign.write.roptions.which = [2 1];
                matlabbatch{1}.spm.spatial.realign.write.roptions.interp = 4;
                matlabbatch{1}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
                matlabbatch{1}.spm.spatial.realign.write.roptions.mask = 1;
                matlabbatch{1}.spm.spatial.realign.write.roptions.prefix = 'r';
                spm_jobman('run',matlabbatch);
                
                delete(fullfile(ffxDir, ['r4D_beta*', '.nii']));
                delete(fullfile(ffxDir, ['mean4D_beta*', '.nii']));
                
                % call the realigned image
                realignMaskName = [prefix,maskName];
                realignMaskPath = fullfile(parcelPath, subID,realignMaskName);

%                 resliceRealignMask = fullfile(parcelPath, subID,[prefix,realignMaskName]);
%                 resliceRealignMaskName = [prefix, realignMaskName];
%                 
%                 %reslice so that the resolution would be same as 4D image
%                 reslice_nii(realignMask, resliceRealignMask);
                %% after reslicing, turn it again into binary mask
                %load roi
                realignMask = load_untouch_nii(realignMaskPath); % 
                
                % binarise
                realignMask.img(realignMask.img < threshold) = 0.0;
                realignMask.img(realignMask.img > 0) = 1.0;
                voxelNb = sum(realignMask.img(:));
                
                %save
                save_untouch_nii(realignMask,realignMaskPath);
                
                % save voxel info into a struct
                info(count).SUBname = subID;
                info(count).ROIname = realignMaskName; % resliceRealignMaskName
                info(count).ROInumVox = voxelNb;
                info(count).problem = sum(realignMask.img(realignMask.img~=1));
                count = count +1;
                
            end
        end
        
    case 4
        % do we want to smooth the masks a bit for having a bigger ROIs?
        
        %% decimalise it before smoothing
        count = 1;
        
        for iSub = 1:length(opt.subjects)
            
            % get subject folder name
            subID = ['sub-', opt.subjects{iSub}];
            
            for iMask = 1: length(maskToSmooth)
                
                %choose the mask to realign and reslice
                maskName = maskToSmooth{iMask};
                
                % load the mask
                mask = load_untouch_nii(fullfile(parcelPath, subID, maskName));
                
                % decimalise the mask
                mask.hdr.dime.datatype= 16;
                mask.hdr.dime.bitpix= 16;
                
                
                decimalMask{count,1} = fullfile(parcelPath, subID, ...
                                                [prefixDec,maskName]);
                
                save_untouch_nii(mask,decimalMask{count,1});
                
                count = count +1;
            end
        end
        
        %% smooth it
        matlabbatch = [];
        matlabbatch = setBatchSmoothing(matlabbatch, ...
                                  decimalMask, ...
                                  maskFWHM, ...
                                  prefixSmooth);
        spm_jobman('run', matlabbatch);

        %% threshold the image
        %load smoothed decimal masks
        count = 1;
        for iSub = 1:length(opt.subjects)
            
            % get subject folder name
            subID = ['sub-', opt.subjects{iSub}];
            
            for iMask = 1: length(maskToSmooth)
                
                %choose the mask to realign and reslice
                smoothMaskName = [prefixSmooth,prefixDec, maskToSmooth{iMask}];
                
                %define output name
                outputMaskName = ['thres',num2str(threshold*100),'_', smoothMaskName];
                outputMask = fullfile(parcelPath, subID, outputMaskName);
                
                % load the mask
                smoothMask = load_untouch_nii(fullfile(parcelPath, subID, smoothMaskName));
                
                %calculate the voxel num
                voxelNb = sum((smoothMask.img(:)>threshold));
                
                %% binarise it
                smoothMask.img((smoothMask.img(:)> threshold)) = 1;
                smoothMask.img((smoothMask.img(:)<= threshold)) = 0;
                
                %% save 
                save_untouch_nii(smoothMask,outputMask);
                
                info(count).SUBname = subID;
                info(count).ROIname = outputMaskName;
                info(count).ROInumVox = voxelNb;
                count = count +1;
            end
        end
end

% 
% count = 1;
% for i=1:length(info)
%     
%     a(i)=info(i).ROInumVox;
%     
%     if mod(i,2) == 0
%         
%         b(count) = a(i) + a(i-1);
%         count = count + 1;
%         
%     end
% end
% 
% a=sort(a);
% bar(a);
% disp ('Min VoxNb');
% min(a)
% disp ('Max VoxNb');
% max(a)
% 
% disp ('Min both hemisphere voxNb');
% min(b)
% disp ('Max both hemisphere voxNb');
% max(b)



end