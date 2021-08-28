%scrit file name atomidfind
%purpose:
%This program is used to find the row containing specified atomid in bondoutdata produced by
%bonds_analysis program
function taratomidrow=atomidfind(atomid,tartrjnum,tarbondnum,bondoutdata)
for i=tartrjnum+1:tartrjnum+tarbondnum
    if atomid==bondoutdata{i,1}
        taratomidrow=i;
        break
    end
end