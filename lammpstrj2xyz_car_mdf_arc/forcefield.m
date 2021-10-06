%scrit file name forcefiled
%purpose:
%This program is used to matche atom element with its forcefield
function ff=forcefield(element,numseq,atomtype)
    ffelement=lower(element);
    for i=1:length(numseq)
        if atomtype==numseq{i}
            ff=ffelement{i};
            break
        end
    end
end