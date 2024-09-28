%scrit file bond_idtopo2typeid
%purpose: This function is used to determine the type and atom id mapping for any line in bonds.* file
%version 1.0， 2024.04.07

function [type_bondnum_type,id_bondnum_id] = bond_idtopo2typeid(ii,bondoutdata,id_type)
type_bondnum_type(1,1) = bondoutdata(ii,2);
type_bondnum_type(1,2) = bondoutdata(ii,3);
id_bondnum_id(1,1) = bondoutdata(ii,1);
id_bondnum_id(1,2) = bondoutdata(ii,2);
id_bondnum_id(1,3) = bondoutdata(ii,3);
for i = 1:bondoutdata(ii,3)
    type_line = find(id_type(:,1) == bondoutdata(ii,i+3));
    type = id_type(type_line,2);
    type_bondnum_type(1,2+i) = type;
    id_bondnum_id(1,3+i) = bondoutdata(ii,i+3);
end
end