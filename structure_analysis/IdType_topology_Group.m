%scrit file IdType_topology_Group
%purpose: This function is used to find out the number of functional groups via bonds.* file
%version 1.0， 2024.04.07

function group_output = IdType_topology_Group(group_target_copy,bondoutdata,id_type,n_bond)
group_output = 0;

%逐行考察group文件每个基团的bond信息是否被bonds.*文件每一行为伊始的挖掘拓扑信息满足与否
if size(group_target_copy,1) < 2 %group文件每个基团的信息至少有两行
    fprintf('\n')
    error('search_group_inform.txt文件输入格式有误，某基团的信息少于2行，请检查！')
end

%解析group文件，建立相应的关键字解析文件
group_target_copy_bondedtype = [];
table_444n = [];
table_555n = [];
table_6666 = [];
table_888 = [];
for ii = 1:size(group_target_copy,1)
    if group_target_copy(ii,1) ~= 9999
        group_target_copy_bondedtype_temp = Group_bondedtype(group_target_copy,ii);%显式键接原子类型
        [table_444n_temp,table_555n_temp,table_6666_temp,table_888_temp] = Group_KeywordsAnalysis(group_target_copy,ii,id_type); %关键字解析表
    else
        break
    end
    group_target_copy_bondedtype(size(group_target_copy_bondedtype,1)+1:size(group_target_copy_bondedtype,1)+size(group_target_copy_bondedtype_temp,1),1:size(group_target_copy_bondedtype_temp,2)) = group_target_copy_bondedtype_temp;
    table_444n = [table_444n;table_444n_temp];
    table_555n = [table_555n;table_555n_temp];
    table_6666(size(table_6666,1)+1:size(table_6666,1)+size(table_6666_temp,1),1:size(table_6666_temp,2)) = table_6666_temp; %维度可能不一致
    table_888 = [table_888;table_888_temp];
end

%group文件键接原子依据行标记分块（同一中心原子）
line_mark = 1;
loop_control = 1;
table_line_mark = [];
table_line_mark(1,1) = 1;
for i = 2:size(group_target_copy,1)
    if line_mark == group_target_copy(i,1) && loop_control ==1 && group_target_copy(i,1) ~= 9999
        table_line_mark(size(table_line_mark,1)+1,1) = i;
        loop_control = loop_control + 1;
    elseif line_mark ~= group_target_copy(i,1) && group_target_copy(i,1) ~= 9999
        table_line_mark(size(table_line_mark,1)+1,1) = i;
        line_mark = group_target_copy(i,1);
        loop_control = loop_control + 1;
    end
end
for i = 1:size(table_line_mark)-1
    table_line_mark(i,2) = table_line_mark(i+1,1)-1;
end
table_line_mark(size(table_line_mark,1),2) = size(group_target_copy,1)-1;

%分析group文件，建立所有可能的原子type拓扑信息
topo_deep = size(table_line_mark,1); 
block_num = 1;
Atom_Combination_collect = [];
while block_num <= topo_deep
    group_target_copy_block = group_target_copy(table_line_mark(block_num,1):table_line_mark(block_num,2),:);
    n66 = [];
    for i = 1:size(group_target_copy_block,1)
        n66(size(n66,1)+1,1) = group_target_copy_block(i,3);
    end
    for i = 1:size(group_target_copy_block,1)
        group_target_copy_bondedtype = Group_bondedtype(group_target_copy_block,i);
        table_6666_888_posi = locator_6666_888(group_target_copy,group_target_copy_block,i,table_line_mark,block_num);%获得6666和888在group_target_copy_block考察行的位置列表
        BondedAtom_Combination = Type_Combination(group_target_copy_bondedtype,table_6666,table_888,table_6666_888_posi,id_type,n66(i,1));
        Centre_BondedAtom_Combination = [group_target_copy_block(1,2)*ones(size(BondedAtom_Combination,1),1),BondedAtom_Combination]; %首列加上中心原子typy
        Atom_Combination_collect(size(Atom_Combination_collect,1)+1:size(Atom_Combination_collect,1)+size(Centre_BondedAtom_Combination,1),1:size(Centre_BondedAtom_Combination,2)) = Centre_BondedAtom_Combination;
        Atom_Combination_collect(size(Atom_Combination_collect,1)+1,:) = 0; %加0数字行区分每个block中的每一行的可能type组合
    end
    Atom_Combination_collect(size(Atom_Combination_collect,1)+1,:) = -1; %加-1数字行区分每个block
    block_num = block_num + 1;
end
%获得所有block区分线所在行，即-1所在行
line_nege1 = find(Atom_Combination_collect(:,1) == -1);
block_line = [];
block_line(:,2) = line_nege1';
for i = 1:size(block_line,1)
    if i == 1
        block_line(i,1) = 1;
    else 
        block_line(i,1) = block_line(i-1,2)+1;
    end
end
%bondoutdata中逐行依照group文件限制进行拓扑结构检索辨认
for ijk = 1:size(bondoutdata,1)
    [type_bondnum_type,~] = bond_idtopo2typeid(ijk,bondoutdata,id_type); %bonds.*文件中键接type表和atom id表
    if type_bondnum_type(1,1) == group_target_copy(1,2) %与group文件首行中心原子类型相同
        if group_target_copy(1,3) > 66 %检查是否限制的成键数，而非确定的成键数
            datadelimiter = {'66'};
            [C,~] = strsplit(num2str(group_target_copy(1,3)),datadelimiter,'CollapseDelimiters',false);
            if type_bondnum_type(1,2) >= str2num(C{2}) %满足66n限定
                first_line_OK = 'y';
            else
                continue %继续检索bondoutdata下一行
            end
        elseif type_bondnum_type(1,2) == group_target_copy(1,3) %是确定的成键数,满足相等
            first_line_OK = 'y';
        else
            continue %继续检索bondoutdata下一行
        end
    else
        continue %继续检索bondoutdata下一行
    end
      
    if strcmpi(first_line_OK,'y') %开始搜索建立等深度拓扑结构
        topo_deep = size(table_line_mark,1); %group区块数决定了拓扑结构分析深度，每个区块的行数决定了分析宽度（键接原子结构）
        [table_TopoBond,topo_deep_OK] = TopoBond(bondoutdata,ijk,topo_deep,group_target_copy,table_line_mark,n_bond,table_444n); %根据group文件深度从bonds.*文件考察每行开始建立等深的拓扑文件，深度不够则不满足要求
        
        if strcmpi(topo_deep_OK,'y') %满足拓扑深度
            block_num = 1;
            Block_Yes_Num = 0;
            type_bondnum_type_all = [];
            match_block_line_all = {}; %记录bonds.*所有block以及其分行匹配情况，方便处理444n和555n关键词
            while block_num <= topo_deep %按深度循环分析table_TopoBond是否满足group文件要求
                group_target_copy_block_instance =  Atom_Combination_collect(block_line(block_num,1):block_line(block_num,2),:); %group指定的block的所有可能实例
                line_zero = find(group_target_copy_block_instance(:,1) == 0); %建立同一个block中group文件每行指定可能类型分区位置
                block_member_line = [];
                block_member_line(:,2) = line_zero';
                for i = 1:size(block_member_line,1)
                    if i == 1
                        block_member_line(i,1) = 1;
                    else
                        block_member_line(i,1) = block_member_line(i-1,2) + 1;
                    end
                end
                block_member_line(:,2) = block_member_line(:,2) - 1;%0行不包含
                table_TopoBond_block = table_TopoBond(table_line_mark(block_num,1):table_line_mark(block_num,2),:);%假设满足要求，取相同的行组成block进行比对
                if sum(table_TopoBond_block(:,1)) ~= size(table_TopoBond_block,1)*table_TopoBond_block(1,1) %行标记有不同值，表明拓扑不完全匹配
                    break
                end
                
                %逐行和逐块分析group文件和bonds.*文件提取的block是否匹配，考虑行不匹配，同一行元素位置不匹配，关键字带来的多种可能性
                %以bonds.*文件提取的block信息为基础，匹配group文件指定的所有可能性
                type_bondnum_type_copy = [];
                for ii = 1:size(table_TopoBond_block,1)
                    table_TopoBond_block_line = table_TopoBond_block(ii,:);
                    for i = 3:size(table_TopoBond_block_line,2) %去掉尾部填充0，添加成键数
                        if table_TopoBond_block_line(i) == 0
                            if i == size(table_TopoBond_block_line,2)
                                table_TopoBond_block_line(i) = [];
                                break
                            else
                                table_TopoBond_block_line(i:end) = [];
                                break
                            end
                        end
                    end
                    table_TopoBond_block_line_temp = table_TopoBond_block_line(1,3:end); %可能缺少单键键接原子，注意不能用来计算成键数
                    table_TopoBond_block_line(1,3) = bondoutdata(find(table_TopoBond_block_line(1,2) == bondoutdata(:,1)),3); %成键原子数
                    table_TopoBond_block_line(1,4:length(table_TopoBond_block_line_temp)+3) = table_TopoBond_block_line_temp;
                    
                    [type_bondnum_type,~] = bond_idtopo2typeid(find(table_TopoBond_block_line(1,2) == bondoutdata(:,1)),bondoutdata,id_type);%bonds.*提取type拓扑关系
                    type_bondnum_type_copy(size(type_bondnum_type_copy,1)+1,1) = type_bondnum_type(1,1);%bonds.*文件提取的纯type表，用于后续比对
                    type_bondnum_type_copy(size(type_bondnum_type_copy,1),2:length(type_bondnum_type(3:end))+1) = type_bondnum_type(3:end);
                end
                type_bondnum_type_all(size(type_bondnum_type_all,1)+1:size(type_bondnum_type_all,1)+size(type_bondnum_type_copy,1),1:size(type_bondnum_type_copy,2)) = type_bondnum_type_copy;%用于后面辅助处理444n和555n标记
                
                %开始对比bonds.*文件提取的同一block中各行是否均满足group文件指定
                match_block = {}; %记录匹配的block编号
                for i = 1:size(type_bondnum_type_copy,1)
                    type_bondnum_type_copy_i = nonzeros(type_bondnum_type_copy(i,:))'; %去掉0
                    match_block{i,1} = i;
                    match_block{i,2} = [];
                    for j = 1:size(block_member_line,1)
                        group_target_copy_block_instance_member = group_target_copy_block_instance(block_member_line(j,1):block_member_line(j,2),:);
                        group_target_copy_block_instance_member1 = nonzeros(group_target_copy_block_instance_member(1,:))';%去掉0
                        datadelimiter = {'66'};
                        [~,matches] = strsplit(num2str(group_target_copy_block_instance_member1(1,end)),datadelimiter,'CollapseDelimiters',false);
                        if ~isempty(matches)%最后一列检查是否存在66n标记，若有，同一个分行block尾部均一致含66n
                            match_66n = 1;
                        else
                            match_66n = 0;
                        end
                        for k = 1:size(group_target_copy_block_instance_member,1) %处理每一个分行block
                            group_target_copy_block_instance_member_k = nonzeros(group_target_copy_block_instance_member(k,:))';%去掉0
                            if type_bondnum_type_copy_i(1) == group_target_copy_block_instance_member_k(1) %中心原子匹配
                                if match_66n %含有66n标记
                                    if length(group_target_copy_block_instance_member_k) >= 3
                                        if length(type_bondnum_type_copy_i) >= length(group_target_copy_block_instance_member_k) - 1
                                            if prod(ismember(group_target_copy_block_instance_member_k(2:length(group_target_copy_block_instance_member_k)-1),type_bondnum_type_copy_i(2:end)))
                                                match_block{i,2}(length(match_block{i,2})+1) = j;
                                                match_block_line_all{size(match_block_line_all,1)+1,1} = block_num; %记录nonds.*文件提取的block数
                                                match_block_line_all{size(match_block_line_all,1),2} = match_block{i,1}; %记录nonds.*文件提取的拓扑关系该block中行编号
                                                match_block_line_all{size(match_block_line_all,1),3} = match_block{i,2}; %记录nonds.*文件提取的匹配的block分行在group文件同一个block中哪些分行匹配，所有可能分行中心原子和已知键接原子type不变
                                                break
                                            end
                                        end
                                    elseif length(group_target_copy_block_instance_member_k) <= 2 %group文件没写任何原子type，自由度最大，可以匹配任何键接类型
                                        match_block{i,2}(length(match_block{i,2})+1) = j;
                                        match_block_line_all{size(match_block_line_all,1)+1,1} = block_num; %记录nonds.*文件提取的block数
                                        match_block_line_all{size(match_block_line_all,1),2} = match_block{i,1}; %记录nonds.*文件提取的拓扑关系该block中行编号
                                        match_block_line_all{size(match_block_line_all,1),3} = match_block{i,2}; %记录nonds.*文件提取的匹配的block分行在group文件同一个block中哪些分行匹配，所有可能分行中心原子和已知键接原子type不变
                                        break
                                    end
                                else %不存在66n标记
                                    if length(type_bondnum_type_copy_i) == length(group_target_copy_block_instance_member_k)
                                        if length(type_bondnum_type_copy_i) == 1
                                            if ismember(sortrows(type_bondnum_type_copy_i')',sortrows(group_target_copy_block_instance_member_k')','rows') %满足匹配，记录匹配的block编号
                                                match_block{i,2}(length(match_block{i,2})+1) = j;
                                                match_block_line_all{size(match_block_line_all,1)+1,1} = block_num; %记录nonds.*文件提取的block数
                                                match_block_line_all{size(match_block_line_all,1),2} = match_block{i,1};%记录nonds.*文件提取的拓扑关系该block中行编号
                                                match_block_line_all{size(match_block_line_all,1),3} = match_block{i,2}; %记录nonds.*文件提取的匹配的block分行在group文件同一个block中哪些分行匹配，所有可能分行中心原子和已知键接原子type不变
                                                break
                                            end
                                        elseif length(type_bondnum_type_copy_i) >= 2
                                            if ismember(sortrows(type_bondnum_type_copy_i(2:end)')',sortrows(group_target_copy_block_instance_member_k(2:end)')','rows') %满足匹配(中心原子匹配且键接原子匹配)，记录匹配的block编号
                                                match_block{i,2}(length(match_block{i,2})+1) = j;
                                                match_block_line_all{size(match_block_line_all,1)+1,1} = block_num; %记录nonds.*文件提取的block数
                                                match_block_line_all{size(match_block_line_all,1),2} = match_block{i,1}; %记录nonds.*文件提取的拓扑关系该block中行编号
                                                match_block_line_all{size(match_block_line_all,1),3} = match_block{i,2}; %记录nonds.*文件提取的匹配的block分行在group文件同一个block中哪些分行匹配，所有可能分行中心原子和已知键接原子type不变
                                                break
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end 
                end
                
                %检查梳理是否全部匹配
                match_block_all = [1:size(match_block,1)];
                match_matrix_current = zeros(1,length(match_block_all));
                for i = 1:size(match_block,1)
                    match_matrix = ismember(match_block_all,match_block{i,2});
                    match_matrix_current = match_matrix_current + match_matrix; %记录匹配block分行，未匹配的分行数字始终为0，匹配过的由于加和>=1
                    if isempty(match_block{i,2}) %同一block中存在没有匹配的行
                        block_num = topo_deep + 2; %退出while循环
                        break 
                    end
                end
                if prod(match_matrix_current) == 0 %存在未被匹配的group文件指定的block分行编号
                    block_num = topo_deep + 2; %退出while循环
                end
                if block_num <= topo_deep %尚未出现以上不满足的情况，进一步考察是否全匹配
                    Match_YesNo = Match_GroupBlock(match_block);
                    if strcmpi(Match_YesNo,'y')%存在完全匹配
                        if ~isempty(table_444n) || ~isempty(table_555n)%进一步处理444n和555n标记是否满足
                           Match_45n = Match_444n_555n(table_TopoBond,type_bondnum_type_all,Atom_Combination_collect,table_444n,table_555n);
                           if strcmpi(Match_45n,'y') %满足444n和555n限定
                               block_num = block_num + 1;
                               Block_Yes_Num = Block_Yes_Num + 1;
                           else
                               block_num = topo_deep + 2; %退出while循环
                           end
                        else
                            block_num = block_num + 1;
                            Block_Yes_Num = Block_Yes_Num + 1;
                        end
                    else
                        block_num = topo_deep + 2; %退出while循环
                    end
                end  
            end
            if Block_Yes_Num == topo_deep 
                group_output = group_output + 1; %完全匹配，增加数目
            end
        end
    else
        continue
    end
end
    
    
end


