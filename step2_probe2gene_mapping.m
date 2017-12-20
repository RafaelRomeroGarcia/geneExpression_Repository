%%%%%%%%%%%%%%
%Rafael Romero Garcia
%rr480@cam.ac.uk
%6th July 2017
%University of Cambridge

%Use probe data to calculate gene expression of each subject
function step2_individualProbe_to_geneExpression()
path_probe_dir='AIBS_map/downloaded/';

donors_name={'normalized_microarray_donor9861',...
    'normalized_microarray_donor10021',...
    'normalized_microarray_donor12876',...
    'normalized_microarray_donor14380',...
    'normalized_microarray_donor15496',...
    'normalized_microarray_donor15697'};

for ifol=1:numel(donors_name)
    donor_name=donors_name{ifol};
    
    %load Probe data
    path_prob=[path_probe_dir donor_name '/MicroarrayExpression.csv'];
    load([path_probe_dir donor_name '/probe2gene/res_gene_symbol_tonum.mat']);
    load([path_probe_dir donor_name '/probe2gene/probe_id.mat']);
    res_gene_symbol_tonum_all{ifol}=res_gene_symbol_tonum;
    fileID{ifol}=fopen(path_prob);
end
    
    probe_to_gene=zeros(size(res_gene_symbol_tonum))';
    gene_to_probe=zeros(1,max(res_gene_symbol_tonum))';
    gene_to_probe_mean=zeros(1,max(res_gene_symbol_tonum))';
    
    %Total number of genes
    numGenes=numel(unique(res_gene_symbol_tonum))-1;
    
    tline1=fgets(fileID{1});
    tline2=fgets(fileID{2});
    tline3=fgets(fileID{3});
    tline4=fgets(fileID{4});
    tline5=fgets(fileID{5});
    tline6=fgets(fileID{6});

    posProbe=1;
    genes_iter=ones(numGenes,1);
    %For each probe
    while ischar(tline1)
        display(num2str(posProbe));
        %splitted=regexp(tline,'\,','split');
        tline_num1=str2num(tline1);
        tline_num1=tline_num1(2:end);
        
        tline_num2=str2num(tline2);
        tline_num2=tline_num2(2:end);
        
        tline_num3=str2num(tline3);
        tline_num3=tline_num3(2:end);
        
        tline_num4=str2num(tline4);
        tline_num4=tline_num4(2:end);
        
        tline_num5=str2num(tline5);
        tline_num5=tline_num5(2:end);
        
        tline_num6=str2num(tline6);
        tline_num6=tline_num6(2:end);
        
        tline=[tline_num1 tline_num2 tline_num3 tline_num4 tline_num5 tline_num6];
        
        %Gene associated to that Probe
        posGene=res_gene_symbol_tonum(posProbe);
        if posGene>0
            %sum of the expression of all gene
            %associated to that probe
            if mean(tline)>gene_to_probe_mean(posGene)
                gene_to_probe(posGene)=posProbe;
                gene_to_probe_mean(posGene)=mean(tline);
            end
        end
        posProbe=posProbe+1;
        
        tline1=fgets(fileID{1});
        tline2=fgets(fileID{2});
        tline3=fgets(fileID{3});
        tline4=fgets(fileID{4});
        tline5=fgets(fileID{5});
        tline6=fgets(fileID{6});
    
        
    end
    
    for ig=1:numel(gene_to_probe)
        if gene_to_probe(ig)>0
            probe_to_gene(gene_to_probe(ig))=ig;
        end
    end
    
    save([path_probe_dir '/gene_to_probe.mat'],'gene_to_probe');
    save([path_probe_dir '/probe_to_gene.mat'],'probe_to_gene');
end