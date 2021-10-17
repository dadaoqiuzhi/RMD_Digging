%scrit file name eleme_molecule
%purpose:
%This program is used to add atom element to the molecule it bonded to.
%elenummatch finally gives the integral molecular components.
function elenummatch=eleme_molecule(elenummatch,elementname)
[row,~]=size(elenummatch);
for i=1:row
   if strcmp(elementname,elenummatch{i,1}) 
       elenummatch{i,2}=elenummatch{i,2}+1;
       break
   end
end