# VBM_PresurgicalAnalysis
Codes derived from completing a presurgical Voxel based morphometry study with brain tumor patients

Morphometry analyses
Developed by Sandra Gisbert-Munoz (mariasandra.gisbert@esic.edu) and Ileana Quinones (i.quinones@bcbl.eu)
at the Basque Center of Cognition Brain and Language (BCBL)
Donostia - San Sebastian, Spain
This matlab function uses SPM12 (https://www.fil.ion.ucl.ac.uk/spm/)

REQUIREMENTS:
- function 'spm_auto_reorient' (https://github.com/CyclotronResearchCentre/spm_auto_reorient)
- function 'trim_image'
- templates for the SPM batches
- atlas folder with ROIs to analyze in .nii format

INPUTS:
corrDTI = 1 -> if you want to corregister T1 with DTI b0 images (the output of this step will be used in further analysis)
corrDTI = 0 -> if you want to use the T1 as template for corregistration
corrt2 = 1 -> if you have T2 images and want to corregister them to the T1
corrt2 = 0 -> if you don't have T2 images
 
To analyze data using this function:
1) Participant data must be saved following the folder structure:
 - Participants
 - Participant_1_folder
       - T1_folder
           - 0T1_image
       - T2_folder (if T2 image is available)
           - 0T2_image

IMPORTANT! the names of the images must start with a '0'
2) The atlas(es) that you want to use must be into already separated ROIs and saved in NIFTI.

This function follows the steps:

1) For trimming the image, you have to check manually the coordinates where you want to cut each participant and save them in a txt file 'Trim_estructurales.txt' with 2 rows:
        top_coordinate_participant_1 top_coordinate_participant_2 ... top_coordinate_participant_N
        bottom_coordinate_participant_1 bottom_coordinate_participant_2 ... bottom_coordinate_participant_N
2) corrDTI = 1; corregisters T1 to DTI; corregisters post to pre
    3) Corregister T1 and T2 images from the same participant (time = 0 && corr == 1)
  3) Segment T1 image (we recommend manually checking the segmentations)
  4) Smooth segmented tissues with a kernel of [8 8 8]
  5) Estimate volume for Grey Matter, White Matter and CSF and save the values in the Participant_N_folder.txt.
  6) Take the atlas(es) into the participant's native space
  7) Estimate the volume of each of the atlas(es)' regions and save the results in a txt file int he Participant_N_folder.
  8) Save all the results together into a table
