%scrit file name lpnumfind
%purpose:
%This program is used to find the row containing lpnum of specified atomid in bondoutdata produced by
%bonds_analysis program
function tarlpnum=lpnumfind(lpnum,tarabondnum)
tarlpnum={};k=1;[row,col]=size(tarabondnum);
for i=1:row
    if lpnum==tarabondnum{i,14}
        tarlpnum(k,:)=tarabondnum(i,:);
        k=k+1;
    end
end