%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rafael Romero Garcia
%rr480@cam.ac.uk
%3rd July 2017
%University of Cambridge
%%%%%%%%%%%%%%%INPUT%%%%%%%%%%%%%%%%%%

%Each donor include a parcellation folder that must contain the file parcellation_name 
parcellation_name='500.aparc_dilated_cortical2mm.nii.gz';
%parcellation='aparc_dilated_cortical2mm.nii.gz_linear';

%hemiMirror=
%0-Do not mirror samples between hemisphere and consider both hemispheres
%1-Mirror samples and consider only left hemisphere
%2-Mirror samples between hemispheres (left samples are fliped and
%considered in the right hemisphere and viceversa)
hemiMirror=1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display('Download the data and compute the valid probes');
step1_download();

display('Estimate gene expression from probe');
step2_individualProbe_to_geneExpression();

display('Estimate voxels location of each sample in T1freesurfer space of each donor');
step3_createGeneExpressionTable(hemiMirror);

display('Generate gene expression and coexpression matrices');
step4_interpolateRegions(hemiMirror,parcellation_name);

display('Gene expression matrix created successfully');