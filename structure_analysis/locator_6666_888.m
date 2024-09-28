%scrit file locator_6666_888
%purpose: This function is used to locate the position of 8888 and 666
%keywords in one line of the investigated group file
%version 1.0， 2024.04.12

function table_6666_888_posi = locator_6666_888(group_target_copy,group_target_copy_block,ii,table_line_mark,block_num)
table_6666_888_posi(1,1) = 6666;
table_6666_888_posi(2,1) = 888;
for i = table_line_mark(block_num,1):table_line_mark(block_num,2)
    if ismember(group_target_copy_block(ii,:),group_target_copy(i,:))
        line = i;
    end
end
num_6666 = 0;
num_888 = 0;
table_6666_posi = [];
table_888_posi = [];
for i = 1:size(group_target_copy,1)
    if i < line
        for j = 1:size(group_target_copy(i,:),2)
            if group_target_copy(i,j) == 6666
                num_6666 = num_6666 + 1;
            elseif group_target_copy(i,j) == 888
                num_888 = num_888 + 1;
            end
        end
    elseif i == line %储存6666和888的序号
        for j = 1:size(group_target_copy(i,:),2)
            if group_target_copy(i,j) == 6666
                num_6666_ok = num_6666 + 1;
                table_6666_posi(length(table_6666_posi)+1) = num_6666_ok;
            elseif group_target_copy(i,j) == 888
                num_888_ok = num_888 + 1;
                table_888_posi(length(table_888_posi)+1) = num_888_ok;
            end
        end
    else
        break
    end
end
table_6666_888_posi(1,2:length(table_6666_posi)+1) = table_6666_posi;
table_6666_888_posi(2,2:length(table_888_posi)+1) = table_888_posi;
end