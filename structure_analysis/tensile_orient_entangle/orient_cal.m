%scrit file name orient_cal
%purpose:
%This program is used to calculate the orientation degree of PC
function Porient=orient_cal(xyz_Posi,eunitvect)
Pvector=[];
len=size(xyz_Posi,1);
if len<3
    error('\nThe atom number is less than 3 after eliminating some nonsignificant atom, please check it!');
else
    for i=1:len-2%
        Pvector(i,1)=i;
        Pvector(i,2)=xyz_Posi(i+2,3)-xyz_Posi(i,3);
        Pvector(i,3)=xyz_Posi(i+2,4)-xyz_Posi(i,4);
        Pvector(i,4)=xyz_Posi(i+2,5)-xyz_Posi(i,5);
    end
    for i=1:len-2
        vector_product=Pvector(i,2)*eunitvect(1)+Pvector(i,3)*eunitvect(2)+Pvector(i,4)*eunitvect(3);
        vector_modxyz=(Pvector(i,2)*Pvector(i,2)+Pvector(i,3)*Pvector(i,3)+Pvector(i,4)*Pvector(i,4))^0.5;
        vector_modeunitvect=(eunitvect(1)^2+eunitvect(2)^2+eunitvect(3)^2)^0.5;
        Pvector(i,5)=(vector_product/vector_modxyz/vector_modeunitvect)^2;
    end
    Porient=1.5*mean(Pvector(:,5))-0.5;
end
