%scrit file name BondLink
%purpose:
%This function is used to update bond topology according to connected
%relation of chemical bond.
%version 1;2021.10.24

function [atomPBC,atomPBCmember]=BondLink(atomPBC,atomPBCmember,tarBOinform_copyi)
if atomPBC(1,3)>0
    for i=1:atomPBC(1,3)
        count3=size(atomPBC,1);
        if ~ismember(atomPBC(1,i+3),atomPBCmember)
            atomPBC(count3+1,1)=atomPBC(1,i+3);
            atomPBCmember(length(atomPBCmember)+1,1)=atomPBC(count3+1,1);
            for j=1:size(tarBOinform_copyi,1)%
                if tarBOinform_copyi{j,1}==atomPBC(count3+1,1)
                    atomPBC(count3+1,2)=tarBOinform_copyi{j,2};
                    if tarBOinform_copyi{j,3}>1
                        atomPBC(count3+1,3)=tarBOinform_copyi{j,3}-1;
                        if atomPBC(count3+1,3)>=1
                            ii=1;
                            for k=1:tarBOinform_copyi{j,3}
                                if tarBOinform_copyi{j,k+3}~=atomPBC(1,1)%
                                    atomPBC(count3+1,ii+3)=tarBOinform_copyi{j,k+3};
                                    ii=ii+1;
                                end
                            end
                        end
                    else
                        atomPBC(count3+1,:)=[];
                        break
                    end
                end
            end
        end
    end
end
end