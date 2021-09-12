%scrit file name molecuweight
%purpose:
%calculate molecular weight by molecular formula
function moleweigh=molecuweight(classmatch)
chemele={'C','H','O','N','He','Li','Be','B','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Pd','Ag','Cd','In','Sn','Sb','I','Xe','Cs','Ba','Pt','Au','Hg','Pb'};
eleweigh=[12,1,16,14,4,7,9,11,19,20,23,24,27,28,31,32,35,40,39,40,45,48,51,52,55,56,59,59,64,65,70,73,75,79,80,84,106,108,112,115,119,122,127,131,133,137,195,197,201,207];%eleweigh=[12.011;1.008;15.999;14.007]
[a,b]=size(classmatch);
matcheleweigh=[];
for i=1:a
    for j=1:length(chemele)
        if classmatch{i,1}==chemele{j}
          matcheleweigh(i)=j;
        end
    end
end
classmatchseccol=classmatch(:,2);
for i=1:length(classmatchseccol)
    classmatchseccol{i}=str2num(classmatchseccol{i});
end
classmatchseccol=cell2mat(classmatchseccol);
moleweigh=0;
for i=1:length(classmatchseccol)
    moleweigh=moleweigh+eleweigh(matcheleweigh(i))*classmatchseccol(i);
end