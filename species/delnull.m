%scrit file delnull
%purpose:
%This program is used del null
function datadelnull=delnull(C)
datadelnull={};
k=1;
for i=1:length(C)
    if ~isempty(C{i})
        datadelnull{k}=C{i};
        k=k+1;
    end
end