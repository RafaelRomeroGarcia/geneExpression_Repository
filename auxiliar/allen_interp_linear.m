function gene_regional_expression=allen_interp_linear(gene_samples,samples_coor_mni,centroids_parcel)
closest_thr=2; %Number of closet samples taken
numSamples=size(samples_coor_mni,1);
distances_to_centroid=sqrt(sum(power(samples_coor_mni-repmat(centroids_parcel,numSamples,1),2),2));
[s_distance s_pos]=sort(distances_to_centroid);

X_p=centroids_parcel(1);Y_p=centroids_parcel(2);Z_p=centroids_parcel(3);
samplePos=s_pos(1:closest_thr);
X_s=samples_coor_mni(samplePos,1);
Y_s=samples_coor_mni(samplePos,2);
Z_s=samples_coor_mni(samplePos,3);

X_pos=samples_coor_mni(s_pos,1);
posX=find(X_pos>X_p);
pos_X_less=s_pos(posX(1));
X_s(1)=samples_coor_mni(pos_X_less);

posX=find(X_pos<X_p);
pos_X_more=s_pos(posX(1));
X_s(2)=samples_coor_mni(pos_X_more);

gene_samples_x=gene_samples([pos_X_less pos_X_more],:);
gene_regional_expression_X=[];gene_regional_expression_Y=[];gene_regional_expression_Z=[];
if numel(unique(X_s))>1
    gene_regional_expression_X=interp1(X_s,gene_samples(samplePos,:),X_p,'linear','extrap');
end

if numel(unique(Y_s))>1
  gene_regional_expression_Y=interp1(Y_s,gene_samples(samplePos,:),Y_p,'linear','extrap');
end

if numel(unique(Z_s))>1
gene_regional_expression_Z=interp1(Z_s,gene_samples(samplePos,:),Z_p,'linear','extrap');
end

gene_regional_expression=mean([gene_regional_expression_X;gene_regional_expression_Y;gene_regional_expression_Z]);
%,Y_s,Z_s,gene_samples(samplePos,1)',centroids_parcel(1),centroids_parcel(2),centroids_parcel(3),'spline','extrapval','-1');
end
