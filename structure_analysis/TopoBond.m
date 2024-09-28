%scrit file TopoBond
%purpose: This function is used to generate bond information with topo
%depth required by search_group_inform file.
%version 1.0， 2024.04.11

function [table_TopoBond,topo_deep_OK] = TopoBond(bondoutdata,ijk,topo_deep,group_target_copy,table_line_mark,n_bond,table_444n)
topo_deep_OK = 'n';
block_num = 1;
atom_centre = [];
atom_bonded = [];
table_TopoBond = [];
table_TopoBond(1,1) = 0;%首行标记,深度行标记一律用group_target_copy的标记
table_TopoBond(1,2:bondoutdata(ijk,3)+2) = bond_connect(bondoutdata,ijk);%中心原子+键接原子id
atom_centre = table_TopoBond(1,2); %中心原子，避免后面重复检索
atom_bonded(1,1) = 0;%首行标记
atom_bonded(1,2:1+length(table_TopoBond(3:end))) = table_TopoBond(3:end); %中心原子为单键的不计入
block_num = block_num + 1;
if block_num > topo_deep
    block_num = 0; %达到指定深度，为1
    topo_deep_OK = 'y';
    block_num_copy = block_num;
end
    
while block_num
    atom_bonded_temp = [];
    num_table_TopoBond_check1 = size(table_TopoBond,1);
    block_num_copy = block_num-1;
    for j = 1:size(atom_bonded,1) %一行一个block
        for i =2:length(atom_bonded(j,:))
            if atom_bonded(j,i) ~= 0
                line = find(bondoutdata(:,1) == atom_bonded(1,i));
                if bondoutdata(line,3) > n_bond && ~ismember(bondoutdata(line,1),table_TopoBond(:,2)) %中心原子不为单键的才计入
                    table_TopoBond(size(table_TopoBond,1)+1,1) = group_target_copy(table_line_mark(block_num,1),1); %深度行标记一律用group_target_copy的标记
                    table_TopoBond(size(table_TopoBond,1),2:bondoutdata(line,3)+2) = bond_connect(bondoutdata,line);%新的中心原子+键接原子id
                    atom_centre(length(atom_centre)+1) = table_TopoBond(size(table_TopoBond,1),2); %更新atom_centre
                    atom_bonded_temp(length(atom_bonded_temp)+1:length(atom_bonded_temp)+bondoutdata(line,3)) = table_TopoBond(size(table_TopoBond,1),3:bondoutdata(line,3)+2); %未处理atom_bonded临时储存，可能有单键原子
                elseif ismember(bondoutdata(line,1),table_TopoBond(:,2)) %检查和记录444n标记原子id
                    if ismember(group_target_copy(table_line_mark(block_num,1),1),table_444n(:,1)) %block标记中存在444n
                        && bondoutdata(line,2) == table_444n(:,4) %
                        
                        
                        
                    end
                end
            else
                break
            end
        end
        num_table_TopoBond_check2 = size(table_TopoBond,1);%检查是否增加了block行对应的中心原子区块
        if num_table_TopoBond_check2 > num_table_TopoBond_check1
            block_num = block_num + 1;
            block_num_copy = block_num;
        else %键接原子均无法作为新的中心原子增加block行，比如不考虑单键
            block_num = 0;
            break
        end
    end
    if block_num > topo_deep
        block_num = 0; %达到指定深度
    else %更新atom_bonded
        
        atom_bonded = [];
        for i = 1:length(atom_bonded_temp)
            if ~ismember(atom_bonded_temp(i),atom_centre)
                line = find(bondoutdata(:,1) == atom_bonded_temp(i));
                if bondoutdata(line,3) > n_bond && isempty(atom_bonded)%未处理键接原子不为单键，重复键接原子不重复录入
                    atom_bonded(size(atom_bonded,1)+1,1) = group_target_copy(table_line_mark(block_num,1),1);
                    atom_bonded(size(atom_bonded,1),2) = atom_bonded_temp(i);
                elseif bondoutdata(line,3) > n_bond && ~ismember(atom_bonded_temp(i),atom_bonded(:,2)) %%未处理键接原子不为单键，重复键接原子不重复录入
                    atom_bonded(size(atom_bonded,1)+1,1) = group_target_copy(table_line_mark(block_num,1),1);
                    atom_bonded(size(atom_bonded,1),2) = atom_bonded_temp(i);
                end
            end
        end
    end
    
    if isempty(atom_bonded)
        block_num = 0; %无新的键接原子可供分析
    end
end

block_num_copy = block_num_copy - 1;
if block_num_copy == topo_deep
    topo_deep_OK = 'y';
end

end
