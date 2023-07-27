function ForCluster_Segmenting(filter,step,corrDTI,corrt2,reor)
%% VBM Morphometry analyses
% Developed by Ileana Quinones at the Basque Center of Cognition Brain and Language (BCBL)
% Donostia - San Sebastian, Spain
% This matlab function uses SPM12 (https://www.fil.ion.ucl.ac.uk/spm/)
%REQUIREMENTS:
    %function 'Segmenting_Patients_ForAtlas'
    %function 'spm_auto_reorient'
    %templates for the SPM batches
    %atlas folder with ROIS to analyze in nii.format
    % Patient folder with "T1_sag", "T2_sag" and "DIFFUSION" folders
    % "Trim_estructurales.txt" file with slices to cut the structural images.   
%EXAMPLE:
%ForCluster_Segmenting('034*Post4',3,1,1,1)
%INPUTS
%filter= participant 
%step= segmentation step
% step = 1 -- cut and corregister T1-T2-DTI
% step = 2 -- segmentation
    % create tumor mask
    % check segmentation
% step = 3 -- volume calculations
%corrDTI= 1 (corregister T1 to DTI image); 0 (not corregister T1 to DTI)
%corrt2= 1 (corregister T2 to T1 image); 0 (not corregister T2 to T1 image)
%reor= 1(adjustment to the comissure); 0 (not doing adjustment)�

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

%% ----------- Start ------------
display Starting
addpath(genpath('/bcbl/home/public/Presurgical_Ileana/Toolbox/spm12'));
addpath('/bcbl/home/public/Presurgical_Ileana/Templates');
path_atlas = '/bcbl/home/public/Presurgical_Ileana/Matlab_Codes_2Share/Atlas';
atlas = dir(path_atlas);
atlas = atlas(3:end);
matter = [1 2]; %c1* and c2* files which refer to grey and white matter tissues
path_participants = '/bcbl/home/public/Presurgical_Ileana/Morphometry';
cd(path_participants);
participants = dir(filter);

%% ------------ Parallel toolbox for Matlab2014B -----------
display(['', num2str(length(participants))]);
if length(participants)>1 && length(participants)<32
     parobj = parpool('local',length(participants));
elseif length(participants)>12
     parobj = parpool('local',12);
end

%% ------------ Loop per participant ------------
for nsubj = 1: length(participants)
    subject = participants(nsubj).name;
    Segmenting_Patients_ForAtlas(path_participants,subject,step,corrDTI,corrt2,reor,path_atlas,atlas,matter);
    display(['Participant ',subject,' finished successfully']);
end

%% ------------ Close nodes -----------
if length(participants)>1
    delete(parobj);
end
end
