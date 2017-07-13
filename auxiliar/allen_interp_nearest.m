function gene_regional_expression=allen_interp_nearest(gene_samples,samples_coor_mni,centroids_parcel)
numSamples=size(samples_coor_mni,1);
distances_to_centroid=sqrt(sum(power(samples_coor_mni-repmat(centroids_parcel,numSamples,1),2),2));
[dummy posfind]=min(distances_to_centroid);
gene_regional_expression=gene_samples(posfind,:);
end
