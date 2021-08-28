%scrit file name charnum_match
%purpose:
%This program is used to matche atom element with its numeric atom type.
function elementname=charnum_match(element,numseq,atomtype)
if atomtype<1 || atomtype>length(numseq)
    disp('Error!!!atomtype shoule be between 1 and %d!',length(numseq))
    elementname='Error!!!';
    error('illegal atomtype parameters!');
else
    for i=1:length(numseq)
        if atomtype==numseq{i}
            elementname=element{i};
            break
        end
    end
end