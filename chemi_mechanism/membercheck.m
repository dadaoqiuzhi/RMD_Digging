%scrit file membercheck
%purpose:
%This program is used to check member
function [indexcol,indexovlp,indexapp]=membercheck(cellin,cellraw)
indexcol=[];indexovlp=[];indexapp=[];
a=length(cellin);
for i=1:a
    id=ismember(cellraw,cellin{i});
    idsum=sum(id);
    if idsum==0%
        d=length(indexapp);
        indexapp(d+1)=i;
    else
        for j=1:length(id)
            if id(j)==1
                c=length(indexcol);
                indexcol(c+1)=j;
                indexovlp(c+1)=i;
            end
        end
    end
end


