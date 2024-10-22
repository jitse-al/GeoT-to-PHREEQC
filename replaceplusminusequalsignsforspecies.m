%Supporting function for conversion to PHREEQC format. Makes for easier reading/conversion
function speciesname=replaceplusminusequalsignsforspecies(speciesname);
    speciesname=deblank(speciesname);%The following replaces -2, -4, +4 etc to --, ----, ++++
    lengthspname=length(speciesname);
    if lengthspname>3
        if isdigit(speciesname(lengthspname-1))
            if speciesname(lengthspname-2) == '+'
                replacementstring='';
                for j=1:str2num(speciesname(lengthspname-1))
                    replacementstring=[replacementstring '+'];
                end
                speciesname=[speciesname(1:lengthspname-3) replacementstring ''''];
            end
            if speciesname(lengthspname-2) == '-'
                replacementstring='';
                for j=1:str2num(speciesname(lengthspname-1))
                    replacementstring=[replacementstring '-'];
                end
                speciesname=[speciesname(1:lengthspname-3) replacementstring ''''];
            end
        end
    end
    equalsigns = speciesname == '=';
    if sum(equalsigns) > 0
        equalsigns2 = 1:length(equalsigns);
        equalsigns2 = equalsigns2(equalsigns>0);
        m=0;
        for k=1:length(equalsigns2)
            if equalsigns2<length(speciesname)
                speciesname=[speciesname(1:equalsigns2(k)-1+m) '-' '-' speciesname(equalsigns2(k)+1+m:length(speciesname)+m)];
            elseif equalsigns2(k) == length(speciesname)
                speciesname=[speciesname(1:equalsigns2(k)-1+m) '-' '-'];
            end
        end
    end %End of replacements of signs