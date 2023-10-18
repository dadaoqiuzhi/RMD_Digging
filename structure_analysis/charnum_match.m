%scrit file name charnum_match
%purpose:
%This program is used to matche atom element with its numeric atom type.
function elementname=charnum_match(element,numseq,atomtype)
if atomtype<1 || atomtype>4
    disp('Error! only C H O N are involved in the program for now£¬namely atomtype should be within 1-4£¡')
    elementname='Error!';
else
    for i=1:length(numseq)
        if atomtype==numseq{i}
            elementname=element{i};
            break
        end
    end
end