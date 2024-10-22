%Supporting function for conversion to PHREEQC. Selects redox reactions that use the preferred species determined in the beginning of the PHREEQC Conversion script

function samenames=determine_most_important_of_double_reactions(samenames,reactionmatrix,preferredspecies, preferredspecies2);
    samenames2=zeros(size(samenames));
    for i=2:length(samenames);
        if samenames(i) == 1
            if samenames(i-1) == 0
                k=i;
                try
                    while samenames(k) == 1
                        k=k+1;
                    end
                end
                k2=k-1;
                k=i-1;
                samenames2(k:k2)=ones(k2-k+1,1);
                reactionmatrixtemp=reactionmatrix(k:k2,1);
                speciestemp=cell(size(reactionmatrixtemp));
                speciestempvalue=zeros(size(reactionmatrixtemp));
                speciestempcomparison=extractspeciesfromname(reactionmatrixtemp{1}{length(reactionmatrixtemp{1})});
                for j=1:length(reactionmatrixtemp);
                    for l=1:(length(reactionmatrixtemp{j})-1)
                        for h=1:length(preferredspecies2);
                            if strcmp(reactionmatrixtemp{j}{l},preferredspecies2{h})
                                speciestempvalue(j)=speciestempvalue(j)+1;
                            end
                        end
                        speciestemp{j}=[speciestemp{j};extractspeciesfromname(reactionmatrixtemp{j}{l})];
                    end
                    speciestemp{j}=unique(speciestemp{j});
                    speciestempsorter=ones(size(speciestemp{j}));
                    for l=1:length(speciestempsorter)
                        for h=1:length(speciestempcomparison);
                            if strcmp(speciestemp{j}{l},speciestempcomparison{h})
                                speciestempsorter(l)=0;
                            end
                        end
                    end
                    speciestemp{j}=speciestemp{j}(speciestempsorter > 0 );
                    for l=1:length(speciestemp{j});
                        for h=1:length(preferredspecies)
                            if strcmp(speciestemp{j}{l},preferredspecies{h})
                                speciestempvalue(j)=speciestempvalue(j)+1;
                            end
                        end
                    end
                    speciestempvalue(j)=speciestempvalue(j)/l;
                end
                [minvalue,minindex]=max(speciestempvalue);
                samenames2(k+minindex-1)=0;
            end
        end
    end
    samenames=samenames2;