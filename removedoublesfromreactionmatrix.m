%Supporting function for conversion to PHREEQC format. 
function sortingmatrix=removedoublesfromreactionmatrix(reactionmatrix);
    sortingmatrix=zeros(length(reactionmatrix),1);
    for i=2:length(sortingmatrix);
        equaltest=0;
        if length(reactionmatrix{i,1})==length(reactionmatrix{i-1,1});
            equaltest=1;
            for j=1:length(reactionmatrix{i,1})
                equaltest=equaltest.*strcmp(reactionmatrix{i,1}{j},reactionmatrix{i-1,1}{j});
            end
            if equaltest == 1
                equaltest2 = reactionmatrix{i,2}./reactionmatrix{i-1,2};
                equaltest2= equaltest2-max(equaltest2);
                equaltest2=sum(equaltest2);
                if equaltest2==0;
                    sortingmatrix(i)=1;
                end
            end
        end
    end