%scrit file name molecuweight
%purpose:
%This function is used to calculate molecular weight
function moleweigh=molecuweight(classmatch)
chemele={'C','H','O','N'};
eleweigh=[12;1;16;14];%eleweigh=[12.011;1.008;15.999;14.007]
[a,~]=size(classmatch);
matcheleweigh=[];
for i=1:a
    for j=1:length(chemele)
        if classmatch{i,1}==chemele{j}
          matcheleweigh(i)=j;
        end
    end
end
classmatchseccol=classmatch(:,2);
classmatchseccol=cell2mat(classmatchseccol);
moleweigh=0;
for i=1:length(classmatchseccol)
    moleweigh=moleweigh+eleweigh(matcheleweigh(i))*classmatchseccol(i);
end