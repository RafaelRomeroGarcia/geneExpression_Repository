function step4_interpolateRegions(hemiMirror,parcellation_name)


%hemiMirror=
%0-Do not mirror samples between hemisphere
%1-Mirror samples and consider only left hemisphere
%2-Mirror samples between hemispheres (left samples are fliped and
%considered in the right hemisphere and viceversa)


        
path_Allen='AIBS_map/';

%Select mirroring method
if hemiMirror==0
    path_genes_general=[path_Allen 'gene_matrices/'];
elseif hemiMirror==1
    path_genes_general=[path_Allen 'gene_matrices_hemi_mirror/'];
elseif hemiMirror==2
    path_genes_general=[path_Allen 'gene_matrices_mirror/'];
end

[s1 s2 s3]=mkdir(path_genes_general);

interpMethod=1;
%Interpolation method, in case there are regions without samples
if interpMethod==1
    interpMethodName='nearest';
elseif interpMethod==2
    interpMethodName='linear';
elseif interpMethod==3
    interpMethodName='extreme';
end


dir_img=dir([path_genes_general 'raw/']);
parc_cortex=[];gene_samples=[];samples_coor_mni=[];samples_coor_vox=[];

%Estimate which samples lies within each region for each donor
folder_name_list={'normalized_microarray_donor9861','normalized_microarray_donor10021','normalized_microarray_donor12876','normalized_microarray_donor14380','normalized_microarray_donor15496','normalized_microarray_donor15697'};
for img=1:numel(folder_name_list)
    donor_name=folder_name_list{img};
    %Voxels that include a sample
    path_subject=[path_genes_general 'raw/' donor_name '/'];
    samples_coor_vox_single=load([path_subject 'samples_coor_vox.mat']);
    samples_coor_vox_single=samples_coor_vox_single.samples_coor_vox;
    %Add one to each index because matlab is 1-based
    samples_coor_vox_single=samples_coor_vox_single+ones(size(samples_coor_vox_single));
    samples_coor_vox=[samples_coor_vox;samples_coor_vox_single];
    
    gene_samples_single=load([path_subject 'gene_samples.mat']);
    
    %mni coordinates of the samples, just in case need to be interpolated
    samples_coor_mni_single=load([path_subject 'samples_coor_mni.mat']);
    samples_coor_mni_single=samples_coor_mni_single.samples_coor_mni;
    samples_coor_mni=[samples_coor_mni;samples_coor_mni_single];
    
    numGenes=size(gene_samples_single.gene_samples,2);
    gene_samples=[gene_samples;gene_samples_single.gene_samples];
    donor_name=donor_name(23:end);
    path_parcellation=[path_Allen 'Allen_FS/' donor_name '/parcellation/'];
    vol_parc=load_nifti([path_parcellation parcellation_name]);
    vol_samples=vol_parc;
    vol_parc=vol_parc.vol;
    vol_samples.vol=zeros(size(vol_samples.vol));
    vol_cortical=vol_parc;
    
    %centroids_parcels=centroidFromParcellationFun([path_Allen 'Allen_FS/' donor_name '/parcellation/' parcellation_name]);

    centroids_parcels=load('/home/rr480/FSFiles/fsaverageSubP/parcellation/centroids_dk.mat');centroids_parcels=centroids_parcels.centroids_parcels;
    %centroids_parcels=load('/home/rr480/FSFiles/fsaverageSubP/parcellation/500.centroids_only_cortex.txt');
    

    %Total parcel depend whether we choose one or both hemisphere
    if hemiMirror==1
        npar=sum(centroids_parcels(:,1)<0); 
    else
        npar=size(centroids_parcels,1);
    end
    npar_vol=numel(unique(vol_parc))-1;
    vol_cortical(find(vol_cortical>npar))=0;
    
    %Create a nifty file showing the location of each sample in each donor
    %brain
    for ic=1:size(samples_coor_vox_single,1)
        parc_cortex(end+1)=vol_parc(samples_coor_vox_single(ic,1),samples_coor_vox_single(ic,2),samples_coor_vox_single(ic,3));
        value_par=vol_parc(samples_coor_vox_single(ic,1),samples_coor_vox_single(ic,2),samples_coor_vox_single(ic,3));
        for ind1=-1:1:1
            for ind2=-1:1:1
                for ind3=-1:1:1
                    vol_samples.vol(samples_coor_vox_single(ic,1)+ind1,samples_coor_vox_single(ic,2)+ind2,samples_coor_vox_single(ic,3)+ind3)=value_par;
                end
            end
        end
    end
    parc_cortex(find(parc_cortex>npar))=0;
    %Location of the samples for Quality control. It should overlap with
    %the file donors_name/mri/T1.nii
    save_nifti(vol_samples,[path_Allen 'samples_location_' donor_name '.nii.gz']);
    
end
%Calculate proportion of the samples used
if hemiMirror==2
    numSamples=numel(parc_cortex)/2;
    ratio1=parc_cortex(1:numSamples)>0;
    ratio2=parc_cortex(1+numSamples:end)>0;
    parc_cortex_ratio=sum((ratio1+ratio2)>0)/numSamples;
else
    parc_cortex_ratio=sum(parc_cortex>0)/numel(parc_cortex);
    numSamples=numel(parc_cortex);
end
parc_cortex_samples=zeros(npar,1);
gene_regional_expression=nan(npar,numGenes);
%Extract all the samples associated to each parcel
for ip=1:npar
    samplesRegion=find(parc_cortex==ip);
    samplesRegion(find(samplesRegion>numSamples))=samplesRegion(find(samplesRegion>numSamples))-numSamples;
    if isempty(samplesRegion)
        if interpMethod==1
            gene_regional_expression(ip,:)=allen_interp_nearest(gene_samples,samples_coor_mni,centroids_parcels(ip,:));
        elseif interpMethod==2
            gene_regional_expression(ip,:)=allen_interp_linear_interp(gene_samples,samples_coor_mni,centroids_parcels(ip,:));
        elseif interpMethod==3
            gene_regional_expression(ip,:)=9999*ones(1,size(gene_regional_expression,2));
        else
            error
        end
    else
        %Calculate the gene expression of a region as the median gene
        %expression values
        gene_regional_expression(ip,:)=median(gene_samples(samplesRegion,:),1);
        parc_cortex_samples(ip)=numel(samplesRegion);
    end
end

%Save files
path_output=path_genes_general;
gene_regional_correlations=corrcoef(zscore(gene_regional_expression)');
[s1 s2 s3]=mkdir(path_output);
save([path_output 'gene_regional_expression.mat'],'gene_regional_expression');
save([path_output 'gene_regional_correlations.mat'],'gene_regional_correlations');
save([path_output 'parc_cortex_samples.mat'],'parc_cortex_samples');
save([path_output 'parc_cortex_ratio.mat'],'parc_cortex_ratio');
display(['Expression matrices created:' path_output 'gene_regional_expression.mat']);
