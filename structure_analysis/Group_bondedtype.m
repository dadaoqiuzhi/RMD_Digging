%scrit file Group_bondedtype
%purpose: This function is used to generate matched atom type file via the
%group file
%version 1.0， 2024.04.10

function group_target_copy_bondedtype_temp = Group_bondedtype(group_target_copy,ii)
group_target_copy_bondedtype_temp = [];
group_target_copy_line = group_target_copy(ii,:);
num_nonatom = 0;
for j = 4:size(group_target_copy_line,2)
    if group_target_copy_line(1,j) >= 888
        num_nonatom = num_nonatom + 1;
    end
end

if num_nonatom > 0
    if group_target_copy_line(1,3) > num_nonatom  %存在显式的键接原子，用于后续确认考察的bondoutdata数据键接原子type是否满足限定
        for i = 4:size(group_target_copy_line,2)
            if group_target_copy_line(1,i) < 888 && group_target_copy_line(1,i) > 0
                group_target_copy_bondedtype_temp(length(group_target_copy_bondedtype_temp)+1) = group_target_copy_line(1,i);
            else %6666标记后面的限制type未记录
                break
            end
        end
    end
elseif num_nonatom == 0
    for i = 4:size(group_target_copy_line,2)
        if group_target_copy_line(1,i) > 0
            group_target_copy_bondedtype_temp(length(group_target_copy_bondedtype_temp)+1) = group_target_copy_line(1,i);
        elseif group_target_copy_line(1,i) == 0
            break
        end
    end
end
end