%scrit file Match_GroupBlock
%purpose: This function is used to generate all combination of the bonded
%atoms by group restrain
%version 1.0， 2024.04.14

function Match_YesNo = Match_GroupBlock(match_block)
Match_YesNo = 'n';

match_block_num = size(match_block,1);
match_block_all = [1:match_block_num];

if match_block_num == 1
    if match_block_all == match_block{match_block_num,2}
        Match_YesNo = 'y';
    end
elseif match_block_num == 2
    x1 = match_block{1,2};
    x2 = match_block{2,2};
    [x2,x1] = ndgrid(x2,x1);
    combinatorics = [x1(:) x2(:)]; %所有排列组合情况
    for i = 1:size(combinatorics,1)
        if ismember(sortrows(match_block_all')',sortrows(combinatorics(i,:)')','rows')
            Match_YesNo = 'y';
            break
        end
    end
elseif match_block_num == 3
    x1 = match_block{1,2};
    x2 = match_block{2,2};
    x3 = match_block{3,2};
    [x3,x2,x1] = ndgrid(x3,x2,x1);
    combinatorics = [x1(:) x2(:) x3(:)]; %所有排列组合情况
    for i = 1:size(combinatorics,1)
        if ismember(sortrows(match_block_all')',sortrows(combinatorics(i,:)')','rows')
            Match_YesNo = 'y';
            break
        end
    end
elseif match_block_num == 4
    x1 = match_block{1,2};
    x2 = match_block{2,2};
    x3 = match_block{3,2};
    x4 = match_block{4,2};
    [x4,x3,x2,x1] = ndgrid(x4,x3,x2,x1);
    combinatorics = [x1(:) x2(:) x3(:) x4(:)]; %所有排列组合情况
    for i = 1:size(combinatorics,1)
        if ismember(sortrows(match_block_all')',sortrows(combinatorics(i,:)')','rows')
            Match_YesNo = 'y';
            break
        end
    end
elseif match_block_num == 5
    x1 = match_block{1,2};
    x2 = match_block{2,2};
    x3 = match_block{3,2};
    x4 = match_block{4,2};
    x5 = match_block{5,2};
    [x5,x4,x3,x2,x1] = ndgrid(x5,x4,x3,x2,x1);
    combinatorics = [x1(:) x2(:) x3(:) x4(:) x5(:)]; %所有排列组合情况
    for i = 1:size(combinatorics,1)
        if ismember(sortrows(match_block_all')',sortrows(combinatorics(i,:)')','rows')
            Match_YesNo = 'y';
            break
        end
    end
elseif match_block_num == 6
    x1 = match_block{1,2};
    x2 = match_block{2,2};
    x3 = match_block{3,2};
    x4 = match_block{4,2};
    x5 = match_block{5,2};
    x6 = match_block{6,2};
    [x6,x5,x4,x3,x2,x1] = ndgrid(x6,x5,x4,x3,x2,x1);
    combinatorics = [x1(:) x2(:) x3(:) x4(:) x5(:) x6(:)]; %所有排列组合情况
    for i = 1:size(combinatorics,1)
        if ismember(sortrows(match_block_all')',sortrows(combinatorics(i,:)')','rows')
            Match_YesNo = 'y';
            break
        end
    end
elseif match_block_num == 7
    x1 = match_block{1,2};
    x2 = match_block{2,2};
    x3 = match_block{3,2};
    x4 = match_block{4,2};
    x5 = match_block{5,2};
    x6 = match_block{6,2};
    x7 = match_block{7,2};
    [x7,x6,x5,x4,x3,x2,x1] = ndgrid(x7,x6,x5,x4,x3,x2,x1);
    combinatorics = [x1(:) x2(:) x3(:) x4(:) x5(:) x6(:) x7(:)]; %所有排列组合情况
    for i = 1:size(combinatorics,1)
        if ismember(sortrows(match_block_all')',sortrows(combinatorics(i,:)')','rows')
            Match_YesNo = 'y';
            break
        end
    end
elseif match_block_num == 8
    x1 = match_block{1,2};
    x2 = match_block{2,2};
    x3 = match_block{3,2};
    x4 = match_block{4,2};
    x5 = match_block{5,2};
    x6 = match_block{6,2};
    x7 = match_block{7,2};
    x8 = match_block{8,2};
    [x8,x7,x6,x5,x4,x3,x2,x1] = ndgrid(x8,x7,x6,x5,x4,x3,x2,x1);
    combinatorics = [x1(:) x2(:) x3(:) x4(:) x5(:) x6(:) x7(:) x8(:)]; %所有排列组合情况
    for i = 1:size(combinatorics,1)
        if ismember(sortrows(match_block_all')',sortrows(combinatorics(i,:)')','rows')
            Match_YesNo = 'y';
            break
        end
    end
elseif match_block_num == 9
    x1 = match_block{1,2};
    x2 = match_block{2,2};
    x3 = match_block{3,2};
    x4 = match_block{4,2};
    x5 = match_block{5,2};
    x6 = match_block{6,2};
    x7 = match_block{7,2};
    x8 = match_block{8,2};
    x9 = match_block{9,2};
    [x9,x8,x7,x6,x5,x4,x3,x2,x1] = ndgrid(x9,x8,x7,x6,x5,x4,x3,x2,x1);
    combinatorics = [x1(:) x2(:) x3(:) x4(:) x5(:) x6(:) x7(:) x8(:) x9(:)]; %所有排列组合情况
    for i = 1:size(combinatorics,1)
        if ismember(sortrows(match_block_all')',sortrows(combinatorics(i,:)')','rows')
            Match_YesNo = 'y';
            break
        end
    end
elseif match_block_num == 10
    x1 = match_block{1,2};
    x2 = match_block{2,2};
    x3 = match_block{3,2};
    x4 = match_block{4,2};
    x5 = match_block{5,2};
    x6 = match_block{6,2};
    x7 = match_block{7,2};
    x8 = match_block{8,2};
    x9 = match_block{9,2};
    x10 = match_block{10,2};
    [x10,x9,x8,x7,x6,x5,x4,x3,x2,x1] = ndgrid(x10,x9,x8,x7,x6,x5,x4,x3,x2,x1);
    combinatorics = [x1(:) x2(:) x3(:) x4(:) x5(:) x6(:) x7(:) x8(:) x9(:) x10(:)]; %所有排列组合情况
    for i = 1:size(combinatorics,1)
        if ismember(sortrows(match_block_all')',sortrows(combinatorics(i,:)')','rows')
            Match_YesNo = 'y';
            break
        end
    end
else
    error('group文件指定的某区块行数（键接原子作为中心原子展开书写行数）超过10，请检查或修改代码！')
end
end

