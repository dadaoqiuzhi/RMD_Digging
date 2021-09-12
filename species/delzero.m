%%scrit file name delzero
%purpose:
%delete molecular weight with number of 0
%version 1;2019.12.14
fprintf('\nPlease copy the molecualr weight data to MWDdata2\n');
% MWDdata2=[];
% MWDdata2(1,:)=MWDdata(1,:);
% MWDdata2(2,:)=MWDdata(n,:);
i=size(MWDdata2,2);j=1;
while j<=i
    if MWDdata2(2,j)==0
        MWDdata2(:,j)=[];
        i=i-1;
    else
        j=j+1;
    end
end
MWDdata2=MWDdata2';
fprintf('\nstatistical result of molecualr weight in MWDdata2. molecular weight vs number\n');