%scrit file name bondnumfind
%purpose:
%This program is used to find the row containing covalent bond numbers of specified atom type in bondoutdata produced by
%bonds_analysis program
function tarabondnum=bondnumfind(bondnum,bonddatacapture)
tarabondnum={};k=1;[row,~]=size(bonddatacapture);
for i=1:row
    if bondnum==bonddatacapture{i,3}
        tarabondnum(k,:)=bonddatacapture(i,:);
        k=k+1;
    end
end