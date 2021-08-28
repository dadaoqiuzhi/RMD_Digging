%scrit file name atomnummolecule_strcat
%purpose:
%This program is used to splice atom and number in generated file of bondorder_deepmining
%to form a molecular formula, only for single molecule with two columns.
function atomnummolecule=atomnummolecule_strcat(tarelenummatch)
atomnummolecule={};%initialization
[tarelenumrow,tarelenumcol]=size(tarelenummatch);
strcatbin={};%temporary data
for j=1:tarelenumrow
    if tarelenummatch{j,2}==1%1 is not shown in the molecular formula
        tarelenummatch{j,2}=[];
    else
        tarelenummatch{j,2}=num2str(tarelenummatch{j,2});
    end
end
for k=1:tarelenumrow%molecular composition
    if tarelenummatch{k,2}=='0'%the element with 0 is not shown in the molecular formula
        strcatbin{1,k}=[];
    else
        strcatbin{1,k}=strcat(tarelenummatch{k,1},tarelenummatch{k,2});
    end
end
if length(strcatbin)>=2
    for ii=2:tarelenumrow
        strcatbin{1,ii}=strcat(strcatbin{1,ii-1},strcatbin{1,ii});
    end
    atomnummolecule=strcatbin{1,ii};
else%one element,e.g. radical and carbon nanotube
    atomnummolecule=strcatbin{1,1};
end
