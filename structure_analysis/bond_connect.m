%scrit file bond_connect
%purpose: This function is used to determine the connected atom id mapping for any line in bonds.* file
%version 1.0, 2024.04.11

function TopoBond = bond_connect(bondoutdata,ijk)
TopoBond = [];
TopoBond(1,1) = bondoutdata(ijk,1);
for i = 1:bondoutdata(ijk,3)
    TopoBond(length(TopoBond)+1) = bondoutdata(ijk,3+i);
end
end