function gene_regional_expression=allen_interp_linear(gene_samples,samples_coor_mni,centroids_parcel)
closest_thr=2; %Number of closet samples taken
numSamples=size(samples_coor_mni,1);
distances_to_centroid=sqrt(sum(power(samples_coor_mni-repmat(centroids_parcel,numSamples,1),2),2));
[s_distance s_pos]=sort(distances_to_centroid);
gene_regional_expression_X=[];gene_regional_expression_Y=[];gene_regional_expression_Z=[];
X_p=centroids_parcel(1);Y_p=centroids_parcel(2);Z_p=centroids_parcel(3);
%samplePos=s_pos(1:closest_thr);
%X_s=samples_coor_mni(samplePos,1);
%Y_s=samples_coor_mni(samplePos,2);
%Z_s=samples_coor_mni(samplePos,3);
try
X_pos=samples_coor_mni(s_pos,1);
posX=find(X_pos>X_p);
pos_X_less=s_pos(posX(1));
X_s(1)=samples_coor_mni(pos_X_less);

posX=find(X_pos<X_p);
pos_X_more=s_pos(posX(1));
X_s(2)=samples_coor_mni(pos_X_more);

gene_samples_x=gene_samples([pos_X_less pos_X_more],:);
gene_regional_expression_X=interp1(X_s,gene_samples_x,X_p,'linear');
end
try
Y_pos=samples_coor_mni(s_pos,1);
posY=find(Y_pos>Y_p);
pos_Y_less=s_pos(posY(1));
Y_s(1)=samples_coor_mni(pos_Y_less);

posY=find(Y_pos<Y_p);
pos_Y_more=s_pos(posY(1));
Y_s(2)=samples_coor_mni(pos_Y_more);

gene_samples_y=gene_samples([pos_Y_less pos_Y_more],:);
gene_regional_expression_Y=interp1(Y_s,gene_samples_y,Y_p,'linear');
end
try
Z_pos=samples_coor_mni(s_pos,1);
posZ=find(Z_pos>Z_p);
pos_Z_less=s_pos(posZ(1));
Z_s(1)=samples_coor_mni(pos_Z_less);

posZ=find(Z_pos<Z_p);
pos_Z_more=s_pos(posZ(1));
Z_s(2)=samples_coor_mni(pos_Z_more);

gene_samples_z=gene_samples([pos_Z_less pos_Z_more],:);
gene_regional_expression_Z=interp1(Z_s,gene_samples_z,Z_p,'linear');
end
gene_regional_expression=[gene_regional_expression_X;gene_regional_expression_Y;gene_regional_expression_Z];
if size(gene_regional_expression,1)>1
    gene_regional_expression=median(gene_regional_expression);
end
%,Y_s,Z_s,gene_samples(samplePos,1)',centroids_parcel(1),centroids_parcel(2),centroids_parcel(3),'spline','extrapval','-1');
end
