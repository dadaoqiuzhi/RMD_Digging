%scrit file BO_match
%purpose:
%make sure the target BO is found
function BOmatchdata=BO_match(tarBOinform,i,elepriority)
idset=[];temporary={};
temporary=tarBOinform(i,:);temporary=BOid_element(temporary,tarBOinform);
id=BOele_rank(temporary,tarBOinform,elepriority);
for i=4:7
    if ~ischar(BOmatchdata{i})
        lenidset=length(idset);
        idset(lenidset+1)=BOmatchdata{i};
    end
end

