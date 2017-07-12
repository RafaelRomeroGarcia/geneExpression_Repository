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
    
    %Total number of genes
    numGenes=numel(unique(res_gene_symbol_tonum))-1;
    fileID=fopen(path_prob);
    tline=fgets(fileID);
    numSamples=numel(strfind(tline,','));
    expression_samples=zeros(numGenes,numSamples);
    posProbe=1;
    genes_iter=ones(numGenes,1);
    %For each probe
    while ischar(tline)
        
        display(num2str(posProbe));
        splitted=regexp(tline,'\,','split');
        tline_num=str2num(tline);
        tline_num=tline_num(2:end);
        %Gene associated to that Probe
        posGene=res_gene_symbol_tonum(posProbe);
        if posGene>0
            %sum of the expression of all gene
            %associated to that probe
            expression_samples(posGene,:)=expression_samples(posGene,:)+tline_num;
            %Number of genes added to that probe
            genes_iter(posGene)=genes_iter(posGene)+1;
        end
        posProbe=posProbe+1;
        tline=fgets(fileID);
    end
    %Expression sample is the average of all the expression probe associated to each gene    
    genes_samples=expression_samples./repmat(genes_iter-1,1,numSamples);
    save([path_probe_dir donor_name '/probe2gene/genes_samples.mat'],'genes_samples');
end