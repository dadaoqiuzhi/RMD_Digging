%scrit file name SysConvert
%purpose:
%This function is used to convert a given number into an expected number
%system. ASCII is used¡£
function target_num=SysConvert(given_num,base)
if given_num==round(given_num) && given_num>0 && base==round(base) && base>0
    control=1;remaindermat=[];
else
    control=0;
    error('\nConversion of number system should be with an integer!!!');
end

while control
    remainder=rem(given_num,base);
    remaindermat(control,1)=control;
    remaindermat(control,2)=remainder;
    control=control+1;
    quotient=(given_num-remainder)/base;
    if quotient==0
        control=0;
    end
    given_num=quotient;
end

atomid='';i=size(remaindermat,1);
while i
    if remaindermat(i,2)<=9
        atomid=strcat(atomid,num2str(remaindermat(i,2)));
    else
        atomid=strcat(atomid,char(remaindermat(i,2)+87));
    end
    i=i-1;
end
target_num=atomid;

