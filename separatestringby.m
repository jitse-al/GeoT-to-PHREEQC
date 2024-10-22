%Supporting function. Separates strings by a delimiter, and puts all separate strings into cells.
function a=separatestringby(string, separationstring)
a={};
startnr=0;
apos=0;
    for i1=1:length(string)
        if string(i1) == separationstring
            startnr=0;
        else
            if startnr==0
                apos=apos+1;
                startnr=1;
                a=[a {string(i1)}];
            else
                a{apos}=[a{apos} string(i1)];
            end
        end
    end
