%scrit file name BondForm
%purpose:
%This program is used to determine the bond formation according to the van
%der Waals Radii (1.15*(r1+r2)/2)
%version 1;2021.10.22

function xyzdata=BondForm(xyzdata,RadiiData,boxsize)
xyzdata(1,6)=0;
xyzdata_copy=xyzdata;

for i=2:size(xyzdata,1)
    xyzdata(i,6)=((xyzdata(i,3)-xyzdata(1,3))^2+(xyzdata(i,4)-xyzdata(1,4))^2+(xyzdata(i,5)-xyzdata(1,5))^2)^0.5;
    %
    breakans=0;
    for j=1:size(RadiiData,1)
        if xyzdata(i,2)+xyzdata(1,2)==RadiiData(j,1)+RadiiData(j,2) && xyzdata(i,2)*xyzdata(1,2)==RadiiData(j,1)*RadiiData(j,2)
            breakans=1;
            if xyzdata(i,6)>RadiiData(j,5)
                for k=1:3
                    if abs(xyzdata(i,k+2)-xyzdata(1,k+2))>=RadiiData(j,5)
                        if abs(xyzdata(i,k+2)+(boxsize(k,2)-boxsize(k,1))-xyzdata(1,k+2))<abs(xyzdata(i,k+2)-xyzdata(1,k+2))
                            xyzdata_copy(i,k+2)=xyzdata(i,k+2)+(boxsize(k,2)-boxsize(k,1));
                        elseif abs(xyzdata(i,k+2)-(boxsize(k,2)-boxsize(k,1))-xyzdata(1,k+2))<abs(xyzdata(i,k+2)-xyzdata(1,k+2))
                            xyzdata_copy(i,k+2)=xyzdata(i,k+2)-(boxsize(k,2)-boxsize(k,1));
                        end
                    end
                end                
            end
        end
        if breakans
            break
        end
    end
    
    
    xyzdata_copy(i,6)=((xyzdata_copy(i,3)-xyzdata_copy(1,3))^2+(xyzdata_copy(i,4)-xyzdata_copy(1,4))^2+(xyzdata_copy(i,5)-xyzdata_copy(1,5))^2)^0.5;
    if xyzdata_copy(i,6)<=xyzdata(i,6)
        if xyzdata_copy(i,6)<=RadiiData(j,5)
            xyzdata(i,:)=xyzdata_copy(i,:);
        elseif xyzdata_copy(i,6)>RadiiData(j,5)*1.5 && xyzdata(i,6)>xyzdata_copy(i,6)
            fprintf('\n')
            warning('Bond length becomes smaller after Unwrap but larger than the criterion!')
            fprintf('\nWarning: bond length after unwrap %.4f is smaller than the original value%.4f, but larger than the 1.5 times of criterion %.4f.\nThe involved atom ids are %d and %d',xyzdata_copy(i,6),xyzdata(i,6),RadiiData(j,5),xyzdata_copy(1,1),xyzdata_copy(i,1))
            xyzdata(i,:)=xyzdata_copy(i,:);
            xyzdata_copy(i,:)
        end
    else
        error('The bond length becomes larger after Unwrap, please check it!')
    end
end
end