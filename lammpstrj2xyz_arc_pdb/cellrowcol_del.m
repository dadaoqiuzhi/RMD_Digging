%scrit file name cellrowcol_del
%purpose:
%This program is used to delet the row or column which meets certain conditions.
%bondoutdata and datapython are aimed to be operated.
%keyword:delrow or delcol;contition:'NaN' or something else according to
%formatting operation.
function bondoutdata=cellrowcol_del(bondoutdata,keyword,condition)
if strcmp(keyword,'delrow')
    index=strcmp(bondoutdata(:,1),condition);
    bondoutdata(find(index(:,1)==1),:)=[];
end

if strcmp(keyword,'delcol')
    [~,col]=size(bondoutdata);k=0;
    for i=1:col
        k=k+1;
        if strcmp(bondoutdata{1,k},condition)
            bondoutdata(:,k)=[];
            k=k-1; 
        end
    end
end