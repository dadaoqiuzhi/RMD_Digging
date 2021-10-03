%scrit file BOid_element
%purpose:
%id2elementnum
function BOcmp=BOid_element(BOcmp,tarBOinform)
[row,~]=size(tarBOinform);
for i=1:4
    if ~ischar(BOcmp{i+1})
        for j=1:row
            if BOcmp{i+1}==tarBOinform{j,1}
                BOcmp{i+1}=tarBOinform{j,2};
            end
        end
    end
end
