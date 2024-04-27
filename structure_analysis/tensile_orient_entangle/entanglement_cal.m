%scrit file name entanglement_cal
%purpose:
%This program is used to calculate the entanglement degree of PC
function [entangle_deg,degnum]=entanglement_cal(xyz_Posi,entanglength)
entangle_deg=[];degnum=0;
len=size(xyz_Posi,1);
if len<2*entanglength+1
    error('\nToo large interval between entangled monomers, please check it!');
else
    vector=[];
    for i=1:len-2*entanglength
        vector(1,1:3)=xyz_Posi(i,3:5)-xyz_Posi(entanglength+i,3:5);
        vector(2,1:3)=xyz_Posi(2*entanglength+i,3:5)-xyz_Posi(entanglength+i,3:5);
        vector_product=vector(1,1)*vector(2,1)+vector(1,2)*vector(2,2)+vector(1,3)*vector(2,3);
        vector_modA=(vector(1,1)^2+vector(1,2)^2+vector(1,3)^2)^0.5;
        vector_modB=(vector(2,1)^2+vector(2,2)^2+vector(2,3)^2)^0.5;
        entangle_deg(i,1)=acosd(vector_product/vector_modA/vector_modB);
        vector=[];
    end
end

entanglenum=0;
for i=1:size(entangle_deg,1)
    if entangle_deg(i,1)<90
        entanglenum=entanglenum+1;    
    end
end
degnum=entanglenum/length(entangle_deg);
