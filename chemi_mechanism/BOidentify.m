%scrit file BOidentify
%purpose:
%This program is used to identify the existence of the same row in
%BOtablerowsetcopy compared with that in tarBOinform
%instantiation of BOtable
function BOidans=BOidentify(BOcmp,BOtablerowsetcopy)
[row,~]=size(BOtablerowsetcopy);BOidans=0;
for i=1:row-1
    if BOtablerowsetcopy{i,1}==BOcmp{1} && BOtablerowsetcopy{i,6}==BOcmp{6} && BOtablerowsetcopy{i,7}==BOcmp{7} && BOtablerowsetcopy{i,8}==BOcmp{8}
        BOidans=1;
        break;
    end
end