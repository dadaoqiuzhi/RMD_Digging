%scrit file name Posi_map
%purpose:
%This program is used to generate the coordinate position of PC by mapping
%the seconsary carbon in benzene ring with the other atoms on backbone
function xyz_Posi=Posi_map(BOrderid,trjdata)
xyz_Posi=[];row=size(BOrderid,1);row2=size(trjdata,1);
k=1;
while row
    if BOrderid(k,3)==0
        for j=1:row2
            if BOrderid(k,2)==trjdata(j,1)
                len=size(xyz_Posi,1);
                xyz_Posi(len+1,1)=len+1;
                xyz_Posi(len+1,2)=BOrderid(k,2);
                xyz_Posi(len+1,3:5)=trjdata(j,3:5);
            end
        end
        k=k+1;
        row=row-1;
    else
        if k+1<=row
            if BOrderid(k+1,3)~=0
                BObenzene=[];
                BObenzene(1:2,1)=BOrderid(k,2:3)';
                BObenzene(3:4,1)=BOrderid(k+1,2:3)';
            else
                BObenzene=[];
                BObenzene(1:2,1)=BOrderid(k,2:3)';
                BObenzene(3,1)=BOrderid(k+1,2)';
                BObenzene(4,1)=BOrderid(k+1,2)';
            end
                for j=1:row2
                    for i=1:4
                        if BObenzene(i,1)==trjdata(j,1)
                            BObenzene(i,2:4)=trjdata(j,3:5);
                        end
                    end
                end
                len=size(xyz_Posi,1);
                xyz_Posi(len+1,1)=len+1;
                xyz_Posi(len+1,3)=mean(BObenzene(:,2));
                xyz_Posi(len+1,4)=mean(BObenzene(:,3));
                xyz_Posi(len+1,5)=mean(BObenzene(:,4));
%             else
%                 error('\nNot a secondary C on successive benze ring, see the %dth line in BOrderid',k+1);
%             end
        end
        k=k+2;
        row=row-2;
    end
end

