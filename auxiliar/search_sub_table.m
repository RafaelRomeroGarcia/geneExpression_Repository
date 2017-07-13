function out=search_sub_table(table_names,table_sub,table_sub_cx)
error=true;
for it=1:numel(table_sub)
    if  not(isempty(strfind(table_names,'nucleus')))      ||  ...
        not(isempty((strfind(table_names,'nuclei'))))     ||  ...
        not(isempty((strfind(table_names,'paravermis')))) ||  ...
        not(isempty((strfind(table_names,'medullary'))))  ||  ...
        not(isempty((strfind(table_names,'olivary'))))    ||  ...
        not(isempty((strfind(table_names,'gigantocellular'))))    ||  ...
        not(isempty((strfind(table_names,'VII'))))      ||  ...
        not(isempty((strfind(table_names,'Crus'))))     ||  ...
        not(isempty((strfind(table_names,'pons'))))    ||  ...
        not(isempty((strfind(table_names,'colliculus'))))      
    
    
    
        %display(table_names);
        out=0;
        error=false;
    end
    
    if strcmp(table_names,table_sub{it})
        out=table_sub_cx(it);
        error=false;
    end
end
if error
   display(table_names)
   %display('area not found');
   FAIL
end
end