%Supporting function for conversion to PHREEQC. For mineral/specie strings, extracts amount of O and H and puts it in brackets at end of string. For easier reading of mineral formulas.

function [elementstotalnew, numberstotalnew]=correctOH(elementstotal, numberstotal)
    elementstotalnew=unique(elementstotal);
    numberstotalnew=zeros(size(elementstotalnew));
    OHpresent=0;
    Hpresent=0;
    Opresent=0;
    for i=1:length(elementstotal);
        for j=1:length(elementstotalnew);
            if strcmp(elementstotal{i},elementstotalnew{j})
                numberstotalnew(j)=numberstotalnew(j)+numberstotal(i);
            end
        end
    end
    for i=1:length(numberstotalnew)
        if strcmp(elementstotalnew{i},'(OH)')
            OHpresent=1;
        end
        if strcmp(elementstotalnew{i},'H')
            Hpresent=1;
        end
        if strcmp(elementstotalnew{i},'O')
            Opresent=1;
        end  
    end
    if OHpresent==0
        elementstotalnew=[elementstotalnew;'(OH)'];
        numberstotalnew=[numberstotalnew;0];
    end
    if Hpresent==0
        elementstotalnew=[elementstotalnew;'H'];
        numberstotalnew=[numberstotalnew;0];
    end
    if Opresent==0
        elementstotalnew=[elementstotalnew;'O'];
        numberstotalnew=[numberstotalnew;0];
    end
    for i=1:length(elementstotalnew);
        switch elementstotalnew{i}
            case '(OH)'
                OHvalue=i;
            case 'H'
                Hvalue=i;
            case 'O'
                Ovalue=i;
            end
        end
numberstotalnew(Ovalue)=numberstotalnew(Ovalue)+numberstotalnew(OHvalue);
numberstotalnew(Hvalue)=numberstotalnew(Hvalue)+numberstotalnew(OHvalue);
numberstotalnew(OHvalue)=0;
if numberstotalnew(Ovalue)>0 & numberstotalnew(Hvalue) > 0
    maxOHvalue=min(numberstotalnew(Hvalue),numberstotalnew(Ovalue));
    numberstotalnew(OHvalue)=numberstotalnew(OHvalue)+maxOHvalue;
    numberstotalnew(Ovalue)=numberstotalnew(Ovalue)-maxOHvalue;
    numberstotalnew(Hvalue)=numberstotalnew(Hvalue)-maxOHvalue;
end
elementstotalnew=elementstotalnew((numberstotalnew==0)<1);
numberstotalnew=numberstotalnew((numberstotalnew==0)<1);
        
