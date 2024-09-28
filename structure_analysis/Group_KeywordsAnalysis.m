%scrit file Group_KeywordsAnalysis
%purpose: This function is used to find out the number of functional groups via bonds.* file
%version 1.0， 2024.04.10

function [table_444n_temp,table_555n_temp,table_6666_temp,table_888_temp] = Group_KeywordsAnalysis(group_target_copy,ii,id_type)
table_444n_temp = [];
table_555n_temp = [];
table_6666_temp = [];
table_888_temp = [];
group_target_copy_line = group_target_copy(ii,:);

if max(group_target_copy_line(1,4:end)) >= 888 %检查是否存在444n、555n、6666和888限定词，并进行预处理
    datadelimiter = {'444','555','6666','888'};
    for j = 1:size(group_target_copy_line,2)-3
        if group_target_copy_line(1,3+j) == 0
            break
        end
        [C,matches]=strsplit(num2str(group_target_copy_line(1,3+j)),datadelimiter,'CollapseDelimiters',false);
        if ~isempty(matches)
            if strcmpi(matches{1},'444') %只有唯一一个中心原子
                table_444n_temp(size(table_444n_temp,1)+1,1) = group_target_copy_line(1,1);%行标记
                table_444n_temp(size(table_444n_temp,1),2) = ii; %记录所在行
                table_444n_temp(size(table_444n_temp,1),3) = group_target_copy_line(1,3+j);%记录444n标记
                table_444n_temp(size(table_444n_temp,1),4) =  group_target_copy_line(1,2);%444n标记对应的中心原子type
            elseif strcmpi(matches{1},'555')
                table_555n_temp(size(table_555n_temp,1)+1,1) = group_target_copy_line(1,1);%行标记
                table_555n_temp(size(table_555n_temp,1),2) = ii; %记录所在行
                table_555n_temp(size(table_555n_temp,1),3) = group_target_copy_line(1,3+j);%记录555n标记,可能多个
                table_555n_temp(size(table_555n_temp,1),4) = group_target_copy_line(1,str2num(C{2})+3);%555n标记对应的键接原子type
            elseif strcmpi(matches{1},'6666') %一行最多只能有1个6666标记
                id_postpone = 0;
                for ia = 4:size(group_target_copy_line,2)
                    if group_target_copy_line(1,ia) == 0
                        break
                    end
                    if group_target_copy_line(1,ia) >= 888
                        id_postpone = id_postpone + 1;
                    end
                    if id_postpone >= 1 && group_target_copy_line(1,ia) < 888 && group_target_copy_line(1,ia) > 0
                        if ~ismember(group_target_copy_line(1,ia),id_type(:,2))
                            error('6666指定的type有误，请检查')
                        end
                        table_6666_temp(1,1) = group_target_copy_line(1,1);%行标记
                        table_6666_temp(1,2) = 6666; %记录6666标记
                        table_6666_temp(1,size(table_6666_temp,2)+1) = group_target_copy_line(1,ia); %记录6666指定的type
                    end
                end
            elseif strcmpi(matches{1},'888') %一行可有多个888标记
                table_888_temp(size(table_888_temp,1)+1,1) = group_target_copy_line(1,1);%行标记
                table_888_temp(size(table_888_temp,1),2) = 888;%记录888标记
                ia = 1;
                table = [];%记录888标记相应的type所有可能值
                while ia
                    if ismember(ia,id_type(:,2))
                        table(1,size(table,2)+1) = ia;
                        ia = ia +1;
                    else
                        ia = 0;
                    end
                end
                table_888_temp(size(table_888_temp,1),3:size(table,2)+2) = table;
            end
        end
    end
end
end