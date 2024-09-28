%scrit file Type_Combination
%purpose: This function is used to generate all combination of the bonded
%atoms by group restrain
%version 1.0， 2024.04.08

function BondedAtom_Combination = Type_Combination(group_target_copy_bondedtype,table_6666,table_888,table_6666_888_posi,id_type,n66i)
%确定是否含6666或888关键词，以及所对应位置
table_6666_posi = [];
table_888_posi = [];
for i = 1:size(table_6666_888_posi,1)
    table_6666_888_posi_temp = table_6666_888_posi(i,:);
    if length(table_6666_888_posi_temp) >= 2
        for j = 2:length(table_6666_888_posi_temp)
            if table_6666_888_posi_temp(j) ~= 0 && i == 1
                table_6666_posi(length(table_6666_posi)+1) = table_6666_888_posi_temp(j);
            elseif table_6666_888_posi_temp(j) ~= 0 && i == 2
                table_888_posi(length(table_888_posi)+1) = table_6666_888_posi_temp(j);
            end
        end
    end
end

temp_6666 = [];
if ~isempty(table_6666_posi) %6666一行最多只有一个关键字
    table_6666_current = table_6666(table_6666_posi(1),:);
    for i = 3:length(table_6666_current)%记录6666指定原子类型
        if table_6666_current(i) > 0
            temp_6666(length(temp_6666)+1) = table_6666_current(i);
        end
    end  
end

table_888_current = [];
if ~isempty(table_888_posi)
    for i = 1:length(table_888_posi)
        table_888_current(size(table_888_current,1)+1,:) = table_888(table_888_posi(i),:);
    end
end
temp_888 = [1:max(id_type(:,2))]; %原子type类型
num_888 = size(table_888_current,1); %888限制原子个数，可能不止一个

if ~isempty(temp_6666) && ~isempty(table_888_posi)
    total_num = length(temp_6666)*length(temp_888)^num_888; %总的可能个数
    type_target = []; %可能的指定type矩阵
    for i = 1:total_num
        if ~isempty(group_target_copy_bondedtype)
            type_target(i,1:length(group_target_copy_bondedtype)) = group_target_copy_bondedtype;
        end
    end
    if num_888 == 1
        x1 = temp_6666; 
        x2 = temp_888;  
        [x2,x1] = ndgrid(x2,x1);
        combinatorics = [x1(:) x2(:)]; %所有排列组合情况
    elseif num_888 == 2
        x1 = temp_6666; 
        x2 = temp_888; 
        x3 = temp_888; 
        [x3,x2,x1] = ndgrid(x3,x2,x1);
        combinatorics = [x1(:) x2(:) x3(:)]; %所有排列组合情况
    elseif num_888 == 3
        x1 = temp_6666; 
        x2 = temp_888; 
        x3 = temp_888; 
        x4 = temp_888; 
        [x4,x3,x2,x1] = ndgrid(x4,x3,x2,x1);
        combinatorics=[x1(:) x2(:) x3(:) x4(:)]; %所有排列组合情况
    elseif num_888 == 4
        x1 = temp_6666;
        x2 = temp_888;
        x3 = temp_888;
        x4 = temp_888;
        x5 = temp_888;
        [x5,x4,x3,x2,x1] = ndgrid(x5,x4,x3,x2,x1);
        combinatorics = [x1(:) x2(:) x3(:) x4(:) x5(:)]; %所有排列组合情况
    else
        error('group文件一行中888字符太多，超过4个，请检查')
    end
    if ~isempty(group_target_copy_bondedtype)
        type_target(1:size(type_target,1),size(type_target,2)+1:size(type_target,2)+size(combinatorics,2)) = combinatorics;
    else
        type_target = combinatorics;
    end
    
elseif ~isempty(temp_6666) && isempty(table_888_posi)
    total_num = length(temp_6666); %总的可能个数
    type_target = []; %可能的指定type矩阵
    for i = 1:total_num
        if ~isempty(group_target_copy_bondedtype)
            type_target(i,1:length(group_target_copy_bondedtype)) = group_target_copy_bondedtype;
        end
    end
    combinatorics= temp_6666'; %所有排列组合情况
    if ~isempty(group_target_copy_bondedtype)
        type_target(1:size(type_target,1),size(type_target,2)+1) = combinatorics;
    else
        type_target = combinatorics;
    end
    
    
elseif isempty(temp_6666) && ~isempty(table_888_posi)
    total_num = length(temp_888)^num_888; %总的可能个数
    type_target = []; %可能的指定type矩阵
    for i = 1:total_num
        if ~isempty(group_target_copy_bondedtype)
            type_target(i,1:length(group_target_copy_bondedtype)) = group_target_copy_bondedtype;
        end
    end
    if num_888 == 1
        x2 = temp_888;  
        combinatorics=x2'; %所有排列组合情况
    elseif num_888 == 2
        x2 = temp_888; 
        x3 = temp_888; 
        [x3,x2] = ndgrid(x3,x2);
        combinatorics=[x2(:) x3(:)]; %所有排列组合情况
    elseif num_888 == 3
        x2 = temp_888; 
        x3 = temp_888; 
        x4 = temp_888; 
        [x4,x3,x2] = ndgrid(x4,x3,x2);
        combinatorics = [x2(:) x3(:) x4(:)]; %所有排列组合情况
    elseif num_888 == 4
        x2 = temp_888;
        x3 = temp_888;
        x4 = temp_888;
        x5 = temp_888;
        [x5,x4,x3,x2] = ndgrid(x5,x4,x3,x2);
        combinatorics = [x2(:) x3(:) x4(:) x5(:)]; %所有排列组合情况
    else
        error('group文件一行中888字符太多，超过4个，请检查')
    end
    if num_888 == 1
        if ~isempty(group_target_copy_bondedtype)
            type_target(1:size(type_target,1),size(type_target,2)+1) = combinatorics;
        else
            type_target = combinatorics;
        end
    elseif num_888 > 1
        if ~isempty(group_target_copy_bondedtype)
            type_target(1:size(type_target,1),size(type_target,2)+1:size(type_target,2)+size(combinatorics,2)) = combinatorics;
        else
            type_target = combinatorics;
        end
    end
    
elseif isempty(temp_6666) && isempty(table_888_posi)
    type_target = group_target_copy_bondedtype;
end

if n66i > 66 %66n标记存在
    datadelimiter = {'66'};
    [C,~] = strsplit(num2str(n66i),datadelimiter,'CollapseDelimiters',false);
    num_bond = str2double(C{2});
    if  num_bond >= size(type_target,2)
        type_target(:,size(type_target,2)+1) = n66i; %66n标记可能具有更多的键接原子，且原子type是不限定的
    end
end

BondedAtom_Combination = type_target;
end