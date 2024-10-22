%Supporting function for conversion to PHREEQC format. Needed for sorting.
function [indices2, solution_master_species]=sortmasterspeciesbychargeandOH(indices,solution_master_species,defined_aqueous_species,nrline,mssize2,mssize1)                        
indices2=zeros(mssize1,1);
    for i = 1:nrline
    for j=indices(i,1):indices(i,2)
        totalcharge=str2num(solution_master_species{j,2});
        [elnames,elnumbers]=extractspeciesfromname(solution_master_species{j,3});
        if length(elnames)>1
            for k=1:length(elnames);
                if strcmp(elnames{k},'H')
                    totalcharge=totalcharge-str2num(elnumbers{k});
                end
            end
            if strcmp(solution_master_species{j,1},'O')<1
                if strcmp(elnames{k},'O')
                    totalcharge=totalcharge+2*str2num(elnumbers{k});
                end
            end
        end
        solution_master_species{j,2}=num2str(totalcharge);
    end
    [chargeorder,chargeindex]=sort(solution_master_species(indices(i,1):indices(i,2),2));
    for j=2:mssize2
        solution_master_species(indices(i,1):indices(i,2),j)=solution_master_species(indices(i,1):indices(i,2),j)(chargeindex);
    end
    for j=indices(i,1):indices(i,2)
        for k=1:length(defined_aqueous_species);
            if strcmp(solution_master_species{j,3},defined_aqueous_species{k})
                indices2(j)=1;
            end
        end
    end
    indices2(indices(i,1):indices(i,2))=indices2(indices(i,1):indices(i,2))<1;
end