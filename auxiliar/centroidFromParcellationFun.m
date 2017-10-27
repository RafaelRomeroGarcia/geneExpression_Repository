function centroids=centroidFromParcellationFun(path_par)
    vol=load_nifti(path_par);
    uvol=unique(vol.vol);
    pos=strfind(path_par,'/');
    path_dir=path_par(1:pos(end));
    if exist([path_dir 'centroids.txt'])
        centroids=load([path_dir 'centroids.txt']);
    else
        system('rm -f centroids.txt');

        for iv=2:numel(uvol)
            val=uvol(iv);
            cmd=['fslstats ' path_par ' -l ' num2str(val-1)  ' -u ' num2str(val+1) ' -c >> centroids.txt'  ];
            system(cmd);        
        end
        centroids=importdata('centroids.txt');

        system(['cp centroids.txt ' path_dir 'centroids.txt']);
        system('rm -f centroids.txt');
    end
    
end
