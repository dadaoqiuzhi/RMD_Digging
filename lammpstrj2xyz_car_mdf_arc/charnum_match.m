%scrit file name charnum_match
%purpose:
%This program is used to matche atom element with its numeric atom type.
function elementname=charnum_match(element,numseq,atomtype)
if atomtype<1 || atomtype>length(numseq)
    fprintf('\nError!!!atomtype shoule be between 1 to %d! Please check atomtype:%d!!!\n',length(numseq),atomtype)
    atomtypestr=num2str(atomtype);Erroratomtype=strcat('Error atomtype:',atomtypestr);
    elementname=Erroratomtype;
    error('Illegal atomtype parameters!');
    return;
else
    for i=1:length(numseq)
        if atomtype==numseq{i}
            elementname=element{i};
            break
        end
    end
end