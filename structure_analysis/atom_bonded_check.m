%scrit file atom_bonded_check
%purpose: This function is used to find out whether bond exists for atoms
%in two group by restrained 444n or 555n, or the same atom exists by 555n
%restriction.
%version 1.0， 2024.04.16

function atom_bonded_YesNo = atom_bonded_check(atom_centre_444n1,atom_centre_444n2,bondoutdata)
atom_bonded_YesNo = 'n';
for i = 1:length(atom_centre_444n1)
    line = find(atom_centre_444n1(i) == bondoutdata(:,1));
    bonded_atom = bondoutdata(line,4:3+bondoutdata(line,3));
    for j = 1:length(bonded_atom)
        if ismember(bonded_atom(j),atom_centre_444n2)
            atom_bonded_YesNo = 'y';
            break
        end
    end
    if strcmpi(atom_bonded_YesNo,'y')
        break
    end
end
if sum(ismember(atom_centre_444n1,atom_centre_444n2)) >= 1 %
    atom_bonded_YesNo = 'y';
end
end