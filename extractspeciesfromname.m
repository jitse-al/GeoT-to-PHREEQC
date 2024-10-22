%Supporting function for conversion to PHREEQC. Extract individaul elements from a string containing a chemical formula.

function [elements,numbers]=extractspeciesfromnames(speciesname);                        
    strtest1=isdigit(speciesname);
    strtest2=speciesname=='.';
    strtest3=strtest1+strtest2;
    strtest4=lower(speciesname);
    strtest4= strtest4 == speciesname;
    elementposition=[];
    for i = 1:length(strtest4);
        if strtest4(i) == 0
            elementposition = [elementposition; i];
        end
    end
    
    strtest5=speciesname == '(';
    strtest6=speciesname == ')';
    brackets=[];
    for i = 1:length(strtest4);
        if strtest5(i) == 1
            brackets = [brackets; i 0];
        end
        if strtest6(i) == 1
            [sizea,sizeb]=size(brackets);
            brackets(sizea,2)=i;
        end
    end
    [sizea,sizeb]=size(brackets);
    for i = 1:length(elementposition);
        for j=1:sizea
            if elementposition(i)>brackets(j,1);
                if elementposition(i)<brackets(j,2);
                    elementposition(i)=0;
                end
            end
        end
    end
    elementposition=elementposition(elementposition>0);
    strtest5=upper(speciesname) == speciesname;
    elements=cell(size(elementposition));
    for i=1:length(elementposition)
        k=1;
        if elementposition(i)<length(strtest5);
            try
            while strtest5(elementposition(i)+k)<1
                k=k+1;
            end
            end
        end
        elements{i}=speciesname(elementposition(i):elementposition(i)+k-1);
    end
    for i=1:sizea;
        elements=[elements;{speciesname(brackets(i,1):brackets(i,2))}];
        strtest3(brackets(i,1):brackets(i,2))=zeros(1,brackets(i,2)+1-brackets(i,1));
        elementposition=[elementposition; brackets(i,1)];
    end
    numbers=cell(size(elements));
    for i=1:length(strtest3);
        if strtest3(i)==1;
            el=i-elementposition;
            el=el(el>0);
            el=i-min(el);
            for j=1:length(elementposition);
                if elementposition(j)==el;
                    numbers{j}=[numbers{j} speciesname(i)];
                end
            end
        end
    end
    for i=1:length(numbers);
        if isempty(numbers{i})
            numbers{i}='1';
        end
    end
    

                