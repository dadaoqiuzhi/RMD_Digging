%scrit file IdType_Bond
%purpose: This function is used to find out the atom groups for different
%bonds via bonds.* file
%version 1.0， 2024.04.05

function bond_target=IdType_Bond(bond_target,bondoutdata,id_type)
bond_target(:,3) = bond_target(:,1)+bond_target(:,2);
bond_target(:,4) = bond_target(:,1).*bond_target(:,2);
bond_target(:,5) = 0;
for i =1:size(bondoutdata,1)
    atom1_type = find(id_type(:,1) == bondoutdata(i,1));
    atom1_type = id_type(atom1_type,2);
    for j = 1:bondoutdata(i,3)
        atom2_type = find(id_type(:,1) == bondoutdata(i,j+3));
        atom2_type = id_type(atom2_type,2);
        RowIdx = find(ismember(bond_target(:,3:4),[atom1_type+atom2_type atom1_type*atom2_type],'rows'));
        if ~isempty(RowIdx)
            bond_target(RowIdx,5) = bond_target(RowIdx,5) + 1;
        end
    end
end
end