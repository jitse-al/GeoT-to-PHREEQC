%Supporting function. Removes a part of a string, like spaces or brackets.
function newstring=removeparts(string,parts)
lengths1=length(parts);
lengths2=[];
for i=1:lengths1;
    lengths2=[lengths2 length(parts{i})];
end
for j = 1:lengths1;
    if length(string)>=lengths2(j);
        removables=strfind(string,parts{j});
        speciesones=ones(size(string));
        for k=1:length(removables);                    
            speciesones(removables(k):removables(k)+lengths2(j)-1)=zeros(1,lengths2(j));    
        end
        speciesones=speciesones>0;
        string = string(speciesones);
    end
end
newstring=string;