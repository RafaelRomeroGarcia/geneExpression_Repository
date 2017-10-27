%%%%%%%%%%%%%%
%Rafael Romero Garcia
%rr480@cam.ac.uk
%13th Sept 2017
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
Batch=[];
genes_samples_gen=[];colums_name=[];
numSamples_in=1;
display('Creating files for batch correction');
for ifol=1:numel(donors_name)  
    donor_name=donors_name{ifol};  
    display(donor_name);
    genes_samples=load([path_probe_dir donor_name '/probe2gene/genes_samples.mat'],'genes_samples');
    genes_samples=genes_samples.genes_samples;
    genes_samples_gen=[genes_samples_gen genes_samples];
    numSamples=size(genes_samples,2);
    Batch_numSamples=(numSamples_in:numSamples_in+numSamples-1)';
    Batch_numBatch=(ifol*ones(numSamples,1));
    Batch=[Batch;Batch_numSamples Batch_numBatch];
    numSamples_in=numSamples_in+numSamples;
end
genes_samples_gen=[1:1:size(genes_samples_gen,2);genes_samples_gen];
csvwrite([path_probe_dir 'batch.csv'],Batch);
csvwrite([path_probe_dir 'genes_samples_gen.csv'],genes_samples_gen);
system(['sed -i ''1isamples, batch'' ' path_probe_dir 'batch.csv']);

display('Batch correction');

system(['cd ' pwd ' && R CMD BATCH auxiliar/batch_correction.R']);

display('Split back into subjects');
genes_samples_gen_corr=csvread([path_probe_dir 'genes_samples_gen.csv']);
for ifol=1:numel(donors_name)
    donor_name=donors_name{ifol};   
    genes_samples=genes_samples_gen_corr(:,Batch(:,2)==ifol);
    save([path_probe_dir donor_name '/probe2gene/genes_samples_corr.mat'],'genes_samples');
end