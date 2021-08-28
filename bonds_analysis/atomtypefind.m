%scrit file name atomtypefind
%purpose:
%This program is used to find the row containing specified atom type in bondoutdata produced by
%bonds_analysis program
function taratomtyperow=atomtypefind(atomtype,tartrjnum,tarbondnum,bondoutdata)
taratomtyperow=[];k=1;
for i=tartrjnum+1:tartrjnum+tarbondnum
    if atomtype==bondoutdata{i,2}
        taratomtyperow(k)=i;
        k=k+1;
    end
end