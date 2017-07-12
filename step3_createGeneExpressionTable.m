%%%%%%%%%%%%%%
%Rafael Romero Garcia
%rr480@cam.ac.uk
%6th July 2017
%University of Cambridge

function step3_createGeneExpressionTable(hemiMirror)
normalizeSubject=true;
includeSubcortical=false;

path_allen='AIBS_map/';
path_allen_FS='Allen_FS/';
%path_parcellation='/home/rr480/FSFiles/subjects/Allen/';
path_cortSubcorTable=[path_allen path_allen_FS '/StructureList_RRGedit.csv'];

folder_name_list={'normalized_microarray_donor9861','normalized_microarray_donor10021','normalized_microarray_donor12876','normalized_microarray_donor14380','normalized_microarray_donor15496','normalized_microarray_donor15697'};

%Select mirroring method
if hemiMirror==0
    path_output=[path_allen 'gene_matrices/'];
elseif hemiMirror==1
    path_output=[path_allen 'gene_matrices_hemi_mirror/'];
elseif hemiMirror==2
    path_output=[path_allen 'gene_matrices_mirror/'];
end

path_probe_dir=[path_allen 'downloaded/'];
[s1 s2 s3]=mkdir(path_output);

table_sub_cell=importdata(path_cortSubcorTable);
table_sub_cell=table_sub_cell(2:end);
for it=1:numel(table_sub_cell)
    cell_tab=table_sub_cell{it};
    table_sub{it}=cell_tab(4:end-1);
    table_sub_cx(it)=str2num(cell_tab(1));
end


samples_coor_mni=[];
for ifol=1:numel(folder_name_list)
    samples_coor_vox=[];
    gene_samples=[];
    folder_name=folder_name_list{ifol};
    %Load gene expression values for each sample
    genes_samples=load([path_probe_dir folder_name '/probe2gene/genes_samples.mat']);
    genes_samples=genes_samples.genes_samples;
    if normalizeSubject
        genes_samples=zscore(genes_samples);
        
    end
    
    %Load original T1 image
    vol=load_nifti([path_allen '/downloaded/' folder_name '/T1.nii.gz']); 
    
    %Load T1 in Freesurfer native space
    vol_parc=load_nifti([path_allen path_allen_FS  folder_name(23:end) '/mri/T1.nii']);

    %Transformation from voxel to mni coordinates
    ras2vox=inv(vol.vox2ras);
    
    %Sample location
    table=importdata([path_probe_dir folder_name '/SampleAnnot.csv']);
    
    %Sample location names
    table_names={table.textdata{2:end,6}};
    
    %Sample location voxel coordinates
    vox_coor=table.data(:,end-5:end-3);
    
    %Sample location mni coordinates
    mni_coor=vol.vox2ras*[vox_coor ones(size(vox_coor,1),1)]';
    mni_coor=mni_coor(1:3,:)';
        
    %For each sample, include only those that have a valid cortical label
    samples_coor_mni=[];gene_samples_cx=[];
    ind=1;
    for is=1:size(mni_coor,1)
        
        if includeSubcortical || search_sub_table(table_names{is},table_sub,table_sub_cx)
            samples_coor_mni(ind,:)=mni_coor(is,:);
            gene_samples_cx(ind,:)=genes_samples(:,is);
            ind=ind+1;
        end
    end
    
    gene_samples=gene_samples_cx;
    
    %Mirror the samples if option 1 or 2 is selected
    if hemiMirror==1
        samples_coor_mni=[-abs(samples_coor_mni(:,1)),samples_coor_mni(:,2),samples_coor_mni(:,3)];
    elseif hemiMirror==2
        samples_coor_mni_1=[samples_coor_mni(:,1),samples_coor_mni(:,2),samples_coor_mni(:,3)];
        samples_coor_mni_2=[-samples_coor_mni(:,1),samples_coor_mni(:,2),samples_coor_mni(:,3)];
        samples_coor_mni=[samples_coor_mni_1;samples_coor_mni_2];
        gene_samples=[gene_samples;gene_samples];
    end
    
    %Calculate the voxel coordinate of each valid sample
    for is=1:size(samples_coor_mni,1)
        samples_coor_vox(is,:)=round(inv(vol_parc.vox2ras)*[samples_coor_mni(is,:) 1]');
    end
    
    samples_coor_vox=samples_coor_vox(:,1:3);
    
  
    mkdir([path_output 'raw/' folder_name]);
    %samples_coor_vox contains the voxels location of each sample in T1
    %freesurfer space
    save([path_output 'raw/' folder_name '/samples_coor_vox.mat'],'samples_coor_vox');
    save([path_output 'raw/' folder_name '/gene_samples.mat'],'gene_samples');
    save([path_output 'raw/' folder_name '/samples_coor_mni.mat'],'samples_coor_mni');
    
    display('Done');
end


