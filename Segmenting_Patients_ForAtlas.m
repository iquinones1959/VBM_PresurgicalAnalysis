function Segmenting_Patients_ForAtlas(path_participants,subject,step,corrDTI,corrt2,reor,path_atlas,atlas,matter)
%% Morphometry analyses
% Developed by Sandra Gisbert-Munoz and Ileana Quinones at the Basque Center of Cognition Brain and Language (BCBL)
% Donostia - San Sebastian, Spain
% This matlab function uses SPM12 (https://www.fil.ion.ucl.ac.uk/spm/)

% REQUIREMENTS:
%   - function 'spm_auto_reorient' (https://github.com/CyclotronResearchCentre/spm_auto_reorient)
%   - function 'trim_image'
%   - templates for the SPM batches
%   - atlas folder with ROIs to analyze in .nii format

% INPUTS:
% corrDTI = 1 -> if you want to corregister T1 to DTI
% corrDTI = 0 -> if you want to use the T1 as template for corregistration
% corrt2 = 1 -> if you have T2 images and want to corregister them to the T1
% corrt2 = 0 -> if you don't have T2 images
% reor = 0 -> if you dont want to do the adjustment of the comissure
% reor = 1 -> if you want to do the adjustment of the comissure


% To analyze data using this function:
% 1) Participant data must be saved following the folder structure:
% - Participants
%   - Participant_1_folder
%       - T1_folder
%           - 0T1_image
%       - T2_folder (if T2 image is available)
%           - 0T2_image
% IMPORTANT! the names of the images must start with a '0'
% 2) The atlas(es) that you want to use must be into already separated ROIs
% and saved in NIFTI.

% This function follows the steps:
% step = 1: Preprocesses T1 image by cutting empty space and reorienting the image

%   1.1) For trimming the image, you have to check manually the coordinates
%   where you want to cut each participant and save them in a txt file 'Trim_estructurales.txt' with 2 rows:
%       top_coordinate_participant_1 top_coordinate_participant_2 ... top_coordinate_participant_N
%       bottom_coordinate_participant_1 bottom_coordinate_participant_2 ... bottom_coordinate_participant_N
%   1.2) corrDTI = 1; corregisters T1 to DTI
%         ; corregisters post to pre
%   1.3) Automatically reorients T1 to set (0,0) at the anterior commisure
%   1.4) Corregister T1 and T2 images from the same participant (time = 0 && corr == 1)
% 3) Segment T1 image (we recommend manually checking the segmentations)
% 4) Smooth segmented tissues [8 8 8]
% 5) Estimate volume for Grey Matter, White Matter and CSF and save the values in the Participant_N_folder.
% 6) Take the atlas(es) into the participant's native space
% 7) Estimate the volume of each of the atlas(es)' regions and save the results in a txt file int he Participant_N_folder.
% 8) Save all the results together into a table

%% ----------- Preprocess structural images ------------
tic;
t1folder = dir([path_participants,filesep,subject,filesep,'*T1_sag*']); %
t2folder = dir([path_participants,filesep,subject,filesep,'*T2_sag*']); %
dtifolder = dir([path_participants,filesep,subject,filesep,'*DIFF*']);
if strfind(subject, '_Preop')
    time = 0; %-> pre   
    corr = 0; %-> you are working with pre data or post data, but you don't want to corregister with the pre
else
    time = 1; %-> pre   corr = 0; %-> you are working with pre data or post data, but you don't want to corregister with the pre
    corr = 1; %-> you are working with pre data or post data, but you don't want to corregister with the pre
    patientPre = [path_participants,filesep,subject(1:end-5),'Preop'];
end
%----------Preprocessing Step 1--------------------
if step == 1
    file = fopen([path_participants,filesep,subject,filesep,'Trim_estructurales.txt'],'r');
    trim_all = fscanf(file,'%f');
    fclose(file);
    trim_x = trim_all(1);
    trim_y = trim_all(2);
    x = 256 - trim_x;
    y = 256 - trim_y; % difference between total y and real bottom
    imageT1 = spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^0.*\.nii$');
    trim_img(imageT1,2,x,y);
    if corrt2 == 1
        imageT2 = spm_select('FPList',[path_participants,filesep,subject,filesep,t2folder.name],'^0.*\.nii$');
        trim_img(imageT2,2,x,y);
    end
    
    if corrDTI == 1
        %------------Preprocessing Step 2 - Batch Corregister (DTI-T1)
        matlabbatch = load('Template_Coregister12.mat');
        matlabbatch = matlabbatch.matlabbatch;
        matlabbatch{1}.spm.spatial.coreg.estimate.ref = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,dtifolder.name],'^*.nii'));
        matlabbatch{1}.spm.spatial.coreg.estimate.source = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^tb.*\.nii$'));
        spm_jobman('run', matlabbatch);
        display(['Coregister T1 to DTI Done - ' subject]);
    elseif time ~= 0 && corr == 1
        % corregister post to pre
        matlabbatch = load('Template_Coregister12.mat');
        matlabbatch = matlabbatch.matlabbatch;
        matlabbatch{1}.spm.spatial.coreg.estimate.ref = cellstr(spm_select('FPList',[path_participants,filesep,patientPre.name,filesep,t1folderPre.name],'^tb.*\.nii$'));
        matlabbatch{1}.spm.spatial.coreg.estimate.source = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^tb.*\.nii$'));
        spm_jobman('run', matlabbatch);
        display(['Coregister T1_post to T1_pre Done - ' subject]);
    end
    if corrt2 == 1
        % -------------- Batch Corregister and Reslice T1-T2 ---------------
        matlabbatch = load('Template_Coregister_Reslice.mat');
        matlabbatch = matlabbatch.matlabbatch;
        matlabbatch{1}.spm.spatial.coreg.estwrite.ref = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^tb.*\.nii$'));
        matlabbatch{1}.spm.spatial.coreg.estwrite.source = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t2folder.name],'^tb.*\.nii$'));
        spm_jobman('run', matlabbatch);
        display(['Coregister T2 to T1 Done - ' subject]);
    end  
elseif step == 2
    % -------------- Batch Segment (anatomical) ---------------
    matlabbatch = load ('Template_SegmentT1T2_spm12.mat');
    matlabbatch = matlabbatch.matlabbatch;
    if corrt2 ==  1
        matlabbatch{1}.spm.spatial.preproc.channel(1).vols =  cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^tb.*.nii$'));
        matlabbatch{1}.spm.spatial.preproc.channel(2) = matlabbatch{1}.spm.spatial.preproc.channel(1);
        matlabbatch{1}.spm.spatial.preproc.channel(2).vols =  cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t2folder.name], '^rtb.*.nii$'));
    else
        matlabbatch{1}.spm.spatial.preproc.channel(1).vols =  cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^tb.*.nii$'));
    end
    matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0]; % to save the grey matter output in native space and normalized
    matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [1 1];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0]; % to save the white matter output in native space and normalized
    matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [1 1];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0]; % to save the csf output in native space and normalized
    matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [1 1];
    matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1]; % to save deformations matrixes
    spm_jobman('run', matlabbatch);
    display(['Done segment - Suj_' subject]);
    % -------------- Batch Smooth (anatomical segmentations) ---------------
    matlabbatch = load ('Template_Smooth.mat');
    matlabbatch = matlabbatch.matlabbatch;
    matlabbatch{1}.spm.spatial.smooth.data = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^wc.*.nii$'));
    matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
    spm_jobman('run', matlabbatch);
    display(['Done smooth - Suj_' subject]);
else
    % ------------------ Prepare tumor ------------------
    % ------------------ Binarize tumor mask ------------
    tumor_name = dir([path_participants,filesep,subject,filesep,t1folder.name,filesep,'s*umor*.nii']);
    matlabbatch = load('Template_ImCalc12.mat');
    matlabbatch = matlabbatch.matlabbatch;
    matlabbatch{1}.spm.util.imcalc.input = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^s.*umor.*\.nii$'));
    matlabbatch{1}.spm.util.imcalc.outdir = cellstr([path_participants,filesep,subject,filesep,t1folder.name]);
    matlabbatch{1}.spm.util.imcalc.output = ['bin_',tumor_name.name];
    matlabbatch{1}.spm.util.imcalc.expression = 'i1>eps';
    matlabbatch{1}.spm.util.imcalc.options.interp = 0;
    spm_jobman('run',matlabbatch);
    % -------------- Normalize Tumor to MNI Space --------------
    matlabbatch = load ('Template_NormaliseWriteT1_12.mat');
    matlabbatch = matlabbatch.matlabbatch;
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^y_tb.*.nii$'));
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^bin_s.*umor.*\.nii$'));
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [90,-126,-72;-90,90,108];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'wT_';
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 0;
    spm_jobman('run', matlabbatch);
    % -------------- Reslice tumor to the t1 space --------------
    matlabbatch = load ('Template_Reslice12.mat');
    matlabbatch = matlabbatch.matlabbatch;
    matlabbatch{1}.spm.spatial.coreg.write.ref = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^tb.*.nii$'));
    matlabbatch{1}.spm.spatial.coreg.write.source = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^bin_s.*umor.*\.nii$'));
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0;
    spm_jobman('run', matlabbatch);
    % -------------- Reslice normalized tumor to the normalized-modulated space --------------
    matlabbatch = load ('Template_Reslice12.mat');
    matlabbatch = matlabbatch.matlabbatch;
    matlabbatch{1}.spm.spatial.coreg.write.ref = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^mwc1.*tb.*.nii'));
    matlabbatch{1}.spm.spatial.coreg.write.source = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^wT_bin_s.*umor.*\.nii$'));
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0;
    spm_jobman('run', matlabbatch);
    % -------------- create a mask adding the resection mask and the tumor mask --------------
    if time ~= 0
        matlabbatch = load('Template_ImCalc12.mat');
        matlabbatch = matlabbatch.matlabbatch;
        matlabbatch{1}.spm.util.imcalc.input = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^rbin_s.*umor.*\.nii$'));
        matlabbatch{1}.spm.util.imcalc.input(2,1) = cellstr(spm_select('FPList',[patientPre,filesep,'T1_sag'],'^rbin_s.*umor.*\.nii$'));
        matlabbatch{1}.spm.util.imcalc.outdir = cellstr([path_participants,filesep,subject,filesep,t1folder.name]);
        matlabbatch{1}.spm.util.imcalc.output = 'bin_tumorandres';
        matlabbatch{1}.spm.util.imcalc.expression = '(i1 + i2) > eps';
        matlabbatch{1}.spm.util.imcalc.options.interp = 0;
        spm_jobman('run',matlabbatch);
        %normalize the mask
        matlabbatch = load ('Template_NormaliseWriteT1_12.mat');
        matlabbatch = matlabbatch.matlabbatch;
        matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^y_tb.*.nii$'));
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'bin_tumorandres.nii$'));
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [90,-126,-72;-90,90,108];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'wT_';
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 0;
        spm_jobman('run', matlabbatch);
    end
    % -------------- Estimate tissue volumes --------------
    if time ~= 0
        lesion_nat = spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^bin_tumorandres.nii$');
        resect_nat = spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^rbin_s.*umor.*\.nii$');
    else
        lesion_nat = spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^rbin_s.*umor.*\.nii$');
    end
    seg_nat =  spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^c.*tb.*.nii$');
    ind_seg = strfind(seg_nat(1,:),filesep);
    tissue = seg_nat(:,ind_seg(end)+1:end);
    % -------------- takes out the tumor volume from the tissues - native --------------
    % 1. generate inverse of tumor
    Vt = spm_vol(lesion_nat);
    Vtumor = spm_read_vols(Vt);
    Vtumor_inv = 0;
    for i=1:size(Vtumor,1)
        for j = 1:size(Vtumor,2)
            for k = 1:size(Vtumor,3)
                if Vtumor(i,j,k) > 0
                    Vtumor(i,j,k) = 1;
                    Vtumor_inv(i,j,k) = 0;
                elseif Vtumor(i,j,k) == 0
                    Vtumor_inv(i,j,k) = 1;
                end
            end
        end
    end
    %2. takes out the tumor from the tissues
    Tvol = zeros(1,10);
    for segNum = 1:3
        Vimg = spm_read_vols(spm_vol(seg_nat(segNum,:))) .* Vtumor_inv;
        Vt.fname = [path_participants,filesep,subject,filesep,t1folder.name,filesep,'noTumor_',tissue(segNum,:)];
        Vt.private.dat.fname = Vt.fname;
        spm_write_vol(Vt,Vimg);
    end
    seg_nat_noTumor =  spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^noTumor_.*.nii$');
    % 3. Estimates volume in native space in cm3 for GM, WM, CSF and calculates TIV from them
    for segNum = 1:3
        Tvol(1,segNum) = spm_summarise(spm_vol(seg_nat_noTumor(segNum,:)),'all','litres')*1000;
    end
    % Estimates volume in native space in cm3 for the tumor mask(s)
    if time ~= 0
        for segNum = 1:3
            Tvol(1,3+segNum) = spm_summarise(spm_vol(seg_nat(segNum,:)),spm_vol(lesion_nat),'litres')*1000;
        end
        Tvol(1,7) = sum(Tvol(:,4:6));
        Tvol(1,8) = spm_summarise(spm_vol(resect_nat),'all','litres')*1000;
    else
        for segNum = 1:3
            Tvol(1,3+segNum) = spm_summarise(spm_vol(seg_nat(segNum,:)),spm_vol(lesion_nat),'litres')*1000;
        end
        Tvol(1,7) = sum(Tvol(:,4:6));
    end
    % save table with volumes for GM (1), WM (2), CSF (3), tumor (4), total without tumor (5), total with tumor (6)
    % total without tumor
    Tvol(:,9) = sum(Tvol(:,1:3),2);
    Tvol(:,10) = sum(Tvol(:,1:6),2);
    fid = fopen([path_participants,filesep,subject,filesep,'TIV.txt'],'w');
    labels = {'noTumor_GM','noTumor_WM','noTumor_CSF','lesion_GM','lesion_WM','lesion_CSF','sum_lesion','resection_vol','TIV_notumor','TIV_all'};
    for j = 1:length(labels)
        fprintf(fid,'%s \t', labels{j});
    end
    fprintf(fid,'\n');
    for i = 1:length(Tvol)
        fprintf(fid,'%f \t',Tvol(i));
    end
    fclose(fid);
    save([path_participants,filesep,subject,filesep,'TIV_',date,'.mat'], 'Tvol');
    
    % ATLAS ANALYSES
    % the atlas(es) you want to use for the analyses must be saved in separate folders inside path_atlas.
    % check how many atlas there are in the path
    for num_atlas = 1:length(atlas)
        roi_path = spm_select('CPath',[path_atlas,filesep,atlas(num_atlas).name]);
        %------------ NATIVE SPACE -----------------
        % Convert aal ROIs into native space and create masks
        % rois will be saved into the same path as the MNI atlas
        % ROIs into native space
        matlabbatch = load ('Template_NormaliseWriteT1_12.mat');
        matlabbatch = matlabbatch.matlabbatch;
        matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^iy.*\.nii$'));
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr(spm_select('FPList',roi_path, '^MNI.*\.nii$'));
        matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = [subject, '_'];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1; % Nearest Neighbor = 0; Trilinear = 1
        spm_jobman('run',matlabbatch);
        % Binarize the regions to create masks
        ROIs = dir([roi_path,filesep,subject,'*.nii']);
        for roiNum = 1:length(ROIs)
            matlabbatch = load('Template_ImCalc12.mat');
            matlabbatch = matlabbatch.matlabbatch;
            matlabbatch{1}.spm.util.imcalc.input = cellstr(spm_select('FPList',roi_path,ROIs(roiNum).name));
            matlabbatch{1}.spm.util.imcalc.outdir = cellstr(roi_path);
            matlabbatch{1}.spm.util.imcalc.output = ['Mask_',ROIs(roiNum).name];
            matlabbatch{1}.spm.util.imcalc.expression = 'i1>eps';
            matlabbatch{1}.spm.util.imcalc.options.interp = 0;
            spm_jobman('run',matlabbatch);
        end
        % Reslice rois to T1 space
        matlabbatch = load ('Template_Reslice12.mat');
        matlabbatch = matlabbatch.matlabbatch;
        matlabbatch{1}.spm.spatial.coreg.write.ref = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^tb.*.nii$'));
        matlabbatch{1}.spm.spatial.coreg.write.source = cellstr(spm_select('FPList',roi_path,['Mask_',subject,'*']));
        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0;
        spm_jobman('run', matlabbatch);
        % Estimating volume of masks in native space
        % go through each ROI and calcualte volume
        ROIs = dir([roi_path,filesep,'rMask*',subject,'*.nii']);
        roi_nat =  spm_select('FPList',roi_path,['rMask_',subject,'*']);
        Tvol_roiNat = zeros(1,length(ROIs));
        Tvol_roiNat_all = zeros(1,length(ROIs));
        for roiNum = 1:length(ROIs)
            % voxels of the segmentations WITHOUT tumor
            Tvol_roiNat(1,roiNum) = spm_summarise(spm_vol(seg_nat_noTumor(matter(num_atlas),:)),spm_vol(roi_nat(roiNum,:)),'litres')*1000;
            % all the voxels in the roi mask
            Tvol_roiNat_all(1,roiNum) = spm_summarise(spm_vol(roi_nat(roiNum,:)),'all','litres')*1000;
        end
        fid = fopen([path_participants,filesep,subject,filesep,'Tvol_',atlas(num_atlas).name,'_nat.txt'],'w');
        for i = 1:length(ROIs)
            fprintf(fid,'%s\t',ROIs(i).name(1:end-4));
        end
        fprintf(fid,'\n');
        for i = 1:length(Tvol_roiNat)
            fprintf(fid,'%f \t',Tvol_roiNat(i));
        end
        fclose(fid);
        save([path_participants,filesep,subject,filesep,'Tvol_',atlas(num_atlas).name,'_nat.mat'], 'Tvol_roiNat');
        
        fid = fopen([path_participants,filesep,subject,filesep,'Tvol_',atlas(num_atlas).name,'_nat_all.txt'],'w');
        for i = 1:length(ROIs)
            fprintf(fid,'%s\t',ROIs(i).name(1:end-4));
        end
        fprintf(fid,'\n');
        for i = 1:length(Tvol_roiNat_all)
            fprintf(fid,'%f \t',Tvol_roiNat_all(i));
        end
        fclose(fid);
        save([path_participants,filesep,subject,filesep,'Tvol_',atlas(num_atlas).name,'_nat_all.mat'],'Tvol_roiNat_all');
        % Estimating %ROI affected by the tumor
        rROIs = dir([roi_path,filesep,'rMask*',subject,'*.nii']);
        % multiply ROI.*tumor
        vol_tumorandROI = zeros(1,length(ROIs));
        for roiNum = 1:length(ROIs)
            Vr = spm_vol(roi_nat(roiNum,:));
            tumorandROI = spm_read_vols(spm_vol(lesion_nat)).* spm_read_vols(Vr);
            Vr.fname = [roi_path,filesep,'tumorXroi_',rROIs(roiNum).name];
            Vr.private.dat.fname = Vr.fname;
            spm_write_vol(Vr,tumorandROI);
        end
        
        tumorrROIs = spm_select('FPList',roi_path,['tumorXroi_rMask_',subject,'*']);
        for roiNum = 1:length(ROIs)
            vol_tumorandROI(1,roiNum) = spm_summarise(spm_vol(tumorrROIs(roiNum,:)),'all','litres')*1000;
        end
        fid = fopen([path_participants,filesep,subject,filesep,'Vol_tumorandROI_',atlas(num_atlas).name,'_nat.txt'],'w');
        for i = 1:length(ROIs)
            fprintf(fid,'%s\t',ROIs(i).name(1:end-4));
        end
        fprintf(fid,'\n');
        for i = 1:length(vol_tumorandROI)
            fprintf(fid,'%f \t',vol_tumorandROI(i));
        end
        fclose(fid);
        save([path_participants,filesep,subject,filesep,'Vol_tumorandROI_',atlas(num_atlas).name,'_nat.mat'], 'vol_tumorandROI');
    end
end
% Controlling the time
display('The pipeline ran without problems');
time_out = toc;
display(['The running time is: ' num2str(time_out)]);
end