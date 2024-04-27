%scrit file name delnonbackbone
%purpose:
%This program is used to delet the groups,atoms which can not contribute to
%the orientation of backbone of polycarbonate,some unreasonable bonds are
%ruled out. The sequence of deletion operation is of importance!
function [backboneBO,delatomid]=delnonbackbone(tarBOinformcopy,element)
delatomid=[];
elenum={'C';'H';'O';'N'};%确定C,H,O,N对应的数值
for i=1:4
    for j=1:4
        if strcmp(element(j),elenum{i,1})
            elenum{i,2}=j;
        end
    end
end

[row,~]=size(tarBOinformcopy);%处理甲基——删除甲基碳键级行
k=1;
while row
    if tarBOinformcopy{k,2}==elenum{1,2}%匹配为碳
        if tarBOinformcopy{k,3}==4
            Cbondatom=[0,0,0,0];[row2,~]=size(tarBOinformcopy);
            for i=1:4 
                for j=1:row2
                    if tarBOinformcopy{k,i+3}==tarBOinformcopy{j,1}%匹配到原子id
                        Cbondatom(i)=tarBOinformcopy{j,2};%记录匹配原子type数值
                    end
                end
                if Cbondatom(i)==0
                    error('\n未匹配到甲基碳所连接的某一原子id！！！');
                end
            end
            if sum(Cbondatom)==elenum{1,2}+elenum{2,2}*3%通过所连接原子类型数值和检验是否为甲基碳
                delatomid(length(delatomid)+1,1)=tarBOinformcopy{k,1};%储存删除原子的id
            end
        end
        row=row-1;
        k=k+1;
    else
        k=k+1;
        row=row-1;
    end
end

[row,~]=size(tarBOinformcopy);%删除甲基碳所在行
while row
    if ismember(tarBOinformcopy{row,1},delatomid)
        tarBOinformcopy(row,:)=[];%删除甲基碳键级行
        row=row-1;
    else
        row=row-1;
    end
end


[row,~]=size(tarBOinformcopy);%处理氧连2个原子，其中一个与另一个氢连接，不合适的水之羟基与其他产物连接
k=1;
while row
    if tarBOinformcopy{k,2}==elenum{3,2} && tarBOinformcopy{k,3}==2%匹配为氧，只连接两个原子
        row2=size(tarBOinformcopy,1);len=size(tarBOinformcopy,1);control=0;
        for i=4:5
            for j=1:row2
                if tarBOinformcopy{j,1}==tarBOinformcopy{k,i} && tarBOinformcopy{j,2}==elenum{2,2}%连接一个氢
                    if i==4%此时tarBOinformcopy{k,5}可能为碳
                        for jj=1:row2
                            if tarBOinformcopy{jj,1}==tarBOinformcopy{k,5} && tarBOinformcopy{jj,2}==elenum{1,2} && tarBOinformcopy{jj,3}==4%苯环季碳
                                delatomid(length(delatomid)+1,1)=tarBOinformcopy{k,1};
                                tarBOinformcopy(k,:)=[];
                                control=1;
                                break;
                            end
                        end
                    else %此时tarBOinformcopy{k,4}可能为碳
                        for jj=1:row2
                            if tarBOinformcopy{jj,1}==tarBOinformcopy{k,4} && tarBOinformcopy{jj,2}==elenum{1,2} && tarBOinformcopy{jj,3}==4%苯环季碳
                                delatomid(length(delatomid)+1,1)=tarBOinformcopy{k,1};
                                tarBOinformcopy(k,:)=[];
                                control=1;
                                break;
                            end
                        end
                    end
                    if control==1
                        break;
                    end
                end
            end
            if control==1
                break;
            end
        end
        row=row-1;
        if len==size(tarBOinformcopy,1);%如果满足条件1而不满足条件2，没有删除氧，处理下一个
            k=k+1;
        end
    else
        row=row-1;
        k=k+1;
    end
end


[row,~]=size(tarBOinformcopy);%处理甲基氢
k=1;%处理H——删除行
while row
    if tarBOinformcopy{k,2}==elenum{2,2}%匹配为氢
        delatomid(length(delatomid)+1,1)=tarBOinformcopy{k,1};
        tarBOinformcopy(k,:)=[];%删除含氢键级行
        row=row-1;
    else
        k=k+1;
        row=row-1;
    end
end


[row,~]=size(tarBOinformcopy);%处理羰基氧——删除行
k=1;
while row
    if tarBOinformcopy{k,2}==elenum{3,2} && tarBOinformcopy{k,3}==1%匹配为氧，只连接一个原子
        row2=size(tarBOinformcopy,1);len=size(tarBOinformcopy,1);
        for i=1:row2
            if tarBOinformcopy{i,1}==tarBOinformcopy{k,4} && tarBOinformcopy{i,2}==elenum{1,2};%原子id匹配为碳
                delatomid(length(delatomid)+1,1)=tarBOinformcopy{k,1};
                tarBOinformcopy(k,:)=[];%删除含羰基氧键级行
                break;
            end
        end
        row=row-1;
        if len==size(tarBOinformcopy,1);%如果满足条件1而不满足条件2，不是羰基氧没有删除氧，处理下一个
            k=k+1;
        end
    else
        k=k+1;
        row=row-1;
    end
end


[row,~]=size(tarBOinformcopy);%处理氧连3个原子，不合适的水等与其他产物连接
k=1;
while row
    if tarBOinformcopy{k,2}==elenum{3,2} && tarBOinformcopy{k,3}==3%匹配为氧，连接三个原子
        delatomid(length(delatomid)+1,1)=tarBOinformcopy{k,1};
        tarBOinformcopy(k,:)=[];
        row=row-1;
    else
        row=row-1;
        k=k+1;
    end
end




[row,~]=size(tarBOinformcopy);%处理氧连2个原子，其中一个与另一个氧连接，可能是不合适的氧气等与其他产物连接，可能删除过氧自由基
k=1;
while row
    if tarBOinformcopy{k,2}==elenum{3,2} && tarBOinformcopy{k,3}==2%匹配为氧，只连接两个原子
        row2=size(tarBOinformcopy,1);control=0;len=size(tarBOinformcopy,1);
        for i=4:5
            for j=1:row2
                if tarBOinformcopy{k,i}==tarBOinformcopy{j,1} && tarBOinformcopy{j,2}==elenum{3,2}%其中一个为氧
                    delatomid(length(delatomid)+1,1)=tarBOinformcopy{k,1};
                    tarBOinformcopy(k,:)=[];
                    control=1;
                    break;
                end
            end
            if control==1%只要找到就跳出所有搜索循环
                break;
            end
        end
        row=row-1;
        if len==size(tarBOinformcopy,1);%如果满足条件1而不满足条件2，没有删除氧，处理下一个
            k=k+1;
        end
    
    else
        row=row-1;
        k=k+1;
    end     
end

backboneBO=tarBOinformcopy;

