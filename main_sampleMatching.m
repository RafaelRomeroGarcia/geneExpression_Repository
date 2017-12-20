%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rafael Romero Garcia
%rr480@cam.ac.uk
%3rd July 2017
%University of Cambridge
%%%%%%%%%%%%%%%INPUT%%%%%%%%%%%%%%%%%%

%Each donor include a parcellation folder (AIBS_map/Allen_FS/donorXXXX/parcellation/) that must contain the file parcellation_name 
parcellation_name='500.aparc_dilated_cortical2mm.nii.gz';
%parcellation_name='aparc_dilated_cortical2mm.nii.gz';

%hemiMirror=
%0-Do not mirror samples between hemisphere and consider both hemispheres
%1-Mirror samples and consider only left hemisphere
%2-Mirror samples between hemispheres (left samples are fliped and
%considered in the right hemisphere and viceversa)
hemiMirror=2;1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display('Downloading the data and computing the valid probes');
step1_download();

display('Mapping probe to gene');
step2_probe2gene_mapping();

display('Estimating gene expression from probe');
step3_individualProbe_to_geneExpression();

display('Estimate voxels location of each sample in T1 freesurfer space of each donor');
step4_createGeneExpressionTable(hemiMirror);

display('Generate gene expression and coexpression matrices');
step5_interpolateRegions(hemiMirror,parcellation_name);

display('Gene expression matrix created successfully');