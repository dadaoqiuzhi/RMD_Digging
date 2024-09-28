%scrit file Match_444n_555n
%purpose: This function is used to decide whether the 444n and 555n
%restrains are all satisfied. 
%version 1.0， 2024.04.10

function Match_45n = Match_444n_555n(table_TopoBond,table_line_mark,table_444n,table_555n,match_block_line_all,bondoutdata)
Match_45n = 'n';

if ~isempty(table_444n)
    line_444n = []; 
    for i = 1:size(table_444n,1)
        if ~ismember(table_444n(i,3),line_444n(:,1))
            line_444n(size(line_444n,1)+1,1) = table_444n(i,3); %444n标记
            line_444n(size(line_444n,1),2) = table_444n(i,4); %中心原子type
            line_444n(size(line_444n,1),3) = table_444n(i,2); %444n在group文件所在行
        else
            exist_line = find(table_444n(i,3) == line_444n(:,1));
            if ismember(0,line_444n(exist_line,:))
                for j = 4:size(line_444n(exist_line,:),2)
                    if line_444n(exist_line,j) == 0
                        line_444n(exist_line,j) = table_444n(i,2);
                        break
                    end
                end
            else
                line_444n(exist_line,size(line_444n(exist_line,:),2)+1) = table_444n(i,2);
            end
        end
    end
    %开始444n匹配性检查
    line_444n = sortrows(line_444n,[1 4]); %按照444n标记和所在行列升序排列，每对444n只能有一对中心原子
    Match_444n_each = [];
    for i = 1:size(line_444n,1)/2
        line1_444n = line_444n(2*i-1,4);
        line2_444n = line_444n(2*i,4);
        for j = 1:size(table_line_mark,1)
            if table_line_mark(j,1) <= line1_444n && line1_444n <= table_line_mark(j,2)
                line_444n(2*i-1,5) = j; %记录444n第一个标记匹配的block_num
                line_444n(2*i-1,6) = line1_444n-j+1; %记录相应block中属于第几行
                break
            end
        end
        for j = 1:size(table_line_mark,1)
            if table_line_mark(j,1) <= line2_444n && line2_444n <= table_line_mark(j,2)
                line_444n(2*i,5) = j; %记录444n第二个标记匹配的block_num
                line_444n(2*i,6) = line2_444n-j+1; %记录相应block中属于第几行
                break
            end
        end
    end
    bonds_match_line1 = [];
    bonds_match_line2 = [];
    for j = 1:size(match_block_line_all,1) %查找bonds.*匹配信息与group文件指定匹配信息自洽
        if match_block_line_all{j,1} == line_444n(2*i-1,5)
            if ismember(line_444n(2*i-1,6),match_block_line_all{j,3}) %匹配行在之前的匹配操作中找到成员
                bonds_match_line1(size(bonds_match_line1)+1) = match_block_line_all{j,2}; %匹配bonds在该block中的行编号
            end
        end
    end
    for j = 1:size(match_block_line_all,1) %查找bonds.*匹配信息与group文件指定匹配信息自洽
        if match_block_line_all{j,1} == line_444n(2*i,5)
            if ismember(line_444n(2*i,6),match_block_line_all{j,3}) %匹配行在之前的匹配操作中找到成员
                bonds_match_line2(size(bonds_match_line2)+1) = match_block_line_all{j,2}; %匹配bonds在该block中的行编号
            end
        end
    end
    atom_centre_444n1 = [];
    atom_centre_444n2 = [];
    table_TopoBond_block1 = table_TopoBond(table_line_mark(line_444n(2*i-1,5),1):table_line_mark(line_444n(2*i-1,5),2),:); %获得bonds444n指定所在block 1
    table_TopoBond_block2 = table_TopoBond(table_line_mark(line_444n(2*i,5),1):table_line_mark(line_444n(2*i,5),2),:);%获得bonds444n指定所在block 2
    for j = 1:length(bonds_match_line1)
        atom_centre_444n1(length(atom_centre_444n1)+1) = table_TopoBond_block1(bonds_match_line1(j),2);
    end
    for j = 1:length(bonds_match_line2)
        atom_centre_444n2(length(atom_centre_444n2)+1) = table_TopoBond_block1(bonds_match_line2(j),2);
    end
    if ~isempty(atom_centre_444n1) && ~isempty(atom_centre_444n2)
        Match_444n_each(length(Match_444n_each)+1) = 0;
        atom_bonded_YesNo = atom_bonded_check(atom_centre_444n1,atom_centre_444n2,bondoutdata);
        if strcmpi(atom_bonded_YesNo,'y')
            Match_444n_each(length(Match_444n_each)) = 1;
        end
        if ismember(0,Match_444n_each)
            warning('444n限定未找到拓扑匹配的中心原子')
        end
    else
        warning('444n限定未找到拓扑匹配的中心原子')
        Match_444n_each(length(Match_444n_each)+1) = 0;
    end
    if prod(Match_444n_each) == 1 %444n全部找到键接原子
        Match_444n = 'y';
    else
        Match_444n = 'n';
    end
    
else
    Match_444n = 'y';
end


if ~isempty(table_555n)
    line_555n = []; 
    for i = 1:size(table_555n,1)
        if ~ismember(table_555n(i,3),line_555n(:,1))
            line_555n(size(line_555n,1)+1,1) = table_555n(i,3); %555n标记
            line_555n(size(line_555n,1),2) = table_555n(i,4); %指定键接原子type
            line_555n(size(line_555n,1),3) = table_555n(i,2); %555n在group文件所在行
        else
            exist_line = find(table_555n(i,3) == line_555n(:,1));
            if ismember(0,line_555n(exist_line,:))
                for j = 4:size(line_555n(exist_line,:),2)
                    if line_555n(exist_line,j) == 0
                        line_555n(exist_line,j) = table_555n(i,2);
                        break
                    end
                end
            else
                line_555n(exist_line,size(line_555n(exist_line,:),2)+1) = table_555n(i,2);
            end
        end
    end
    %开始555n匹配性检查
    
    line_555n = sortrows(line_555n,[1 4]); %按照555n标记和所在行列升序排列，每对555n只能有一对原子相连或者相同
    Match_555n_each = [];
    for i = 1:size(line_555n,1)/2
        line1_555n = line_555n(2*i-1,4);
        line2_555n = line_555n(2*i,4);
        for j = 1:size(table_line_mark,1)
            if table_line_mark(j,1) <= line1_555n && line1_555n <= table_line_mark(j,2)
                line_555n(2*i-1,5) = j; %记录444n第一个标记匹配的block_num
                line_555n(2*i-1,6) = line1_555n-j+1; %记录相应block中属于第几行
                break
            end
        end
        for j = 1:size(table_line_mark,1)
            if table_line_mark(j,1) <= line2_555n && line2_555n <= table_line_mark(j,2)
                line_555n(2*i,5) = j; %记录444n第二个标记匹配的block_num
                line_555n(2*i,6) = line2_444n-j+1; %记录相应block中属于第几行
                break
            end
        end
    end
    bonds_match_line1 = [];
    bonds_match_line2 = [];
    for j = 1:size(match_block_line_all,1) %查找bonds.*匹配信息与group文件指定匹配信息自洽
        if match_block_line_all{j,1} == line_555n(2*i-1,5)
            if ismember(line_555n(2*i-1,6),match_block_line_all{j,3}) %匹配行在之前的匹配操作中找到成员
                bonds_match_line1(size(bonds_match_line1)+1) = match_block_line_all{j,2}; %匹配bonds在该block中的行编号
            end
        end
    end
    for j = 1:size(match_block_line_all,1) %查找bonds.*匹配信息与group文件指定匹配信息自洽
        if match_block_line_all{j,1} == line_555n(2*i,5)
            if ismember(line_555n(2*i,6),match_block_line_all{j,3}) %匹配行在之前的匹配操作中找到成员
                bonds_match_line2(size(bonds_match_line2)+1) = match_block_line_all{j,2}; %匹配bonds在该block中的行编号
            end
        end
    end
    atom_bonded_555n1 = [];
    atom_bonded_555n2 = [];
    table_TopoBond_block1 = table_TopoBond(table_line_mark(line_555n(2*i-1,5),1):table_line_mark(line_555n(2*i-1,5),2),:); %获得bonds444n指定所在block 1
    table_TopoBond_block2 = table_TopoBond(table_line_mark(line_555n(2*i,5),1):table_line_mark(line_555n(2*i,5),2),:);%获得bonds444n指定所在block 2
    for j = 1:length(bonds_match_line1)
        atom_bonded_555n1(length(atom_bonded_555n1)+length(table_TopoBond_block1(bonds_match_line1(j),3:end))) = table_TopoBond_block1(bonds_match_line1(j),3:end);
    end
    for j = 1:length(bonds_match_line2)
        atom_bonded_555n2(length(atom_bonded_555n2)+length(table_TopoBond_block2(bonds_match_line2(j),3:end))) = table_TopoBond_block2(bonds_match_line2(j),3:end);
    end
    if ~isempty(atom_centre_555n1) && ~isempty(atom_centre_555n2)
        Match_555n_each(length(Match_555n_each)+1) = 0;
        atom_bonded_YesNo = atom_bonded_check(atom_bonded_555n1,atom_bonded_555n2,bondoutdata);
        if strcmpi(atom_bonded_YesNo,'y')
            Match_555n_each(length(Match_555n_each)) = 1;
        end
        if ismember(0,Match_555n_each)
            warning('555n限定未找到拓扑匹配的中心原子')
        end
    else
        warning('555n限定未找到拓扑匹配的中心原子')
        Match_555n_each(length(Match_555n_each)+1) = 0;
    end
    if prod(Match_555n_each) == 1 %444n全部找到键接原子
        Match_555n = 'y';
    else
        Match_555n = 'n';
    end

else
    Match_555n = 'y';
end




if strcmpi(Match_444n,'y') && strcmpi(Match_555n,'y')
    Match_45n = 'y';
else
    Match_45n = 'n';
end

end