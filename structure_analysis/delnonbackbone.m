%scrit file name delnonbackbone
%purpose:
%This program is used to delet the groups,atoms which can not contribute to
%the orientation of backbone of polycarbonate,some unreasonable bonds are
%ruled out. The sequence of deletion operation is of importance!
function [backboneBO,delatomid]=delnonbackbone(tarBOinformcopy,element)
delatomid=[];
elenum={'C';'H';'O';'N'};%ȷ��C,H,O,N��Ӧ����ֵ
for i=1:4
    for j=1:4
        if strcmp(element(j),elenum{i,1})
            elenum{i,2}=j;
        end
    end
end

[row,~]=size(tarBOinformcopy);%����׻�����ɾ���׻�̼������
k=1;
while row
    if tarBOinformcopy{k,2}==elenum{1,2}%ƥ��Ϊ̼
        if tarBOinformcopy{k,3}==4
            Cbondatom=[0,0,0,0];[row2,~]=size(tarBOinformcopy);
            for i=1:4 
                for j=1:row2
                    if tarBOinformcopy{k,i+3}==tarBOinformcopy{j,1}%ƥ�䵽ԭ��id
                        Cbondatom(i)=tarBOinformcopy{j,2};%��¼ƥ��ԭ��type��ֵ
                    end
                end
                if Cbondatom(i)==0
                    error('\nδƥ�䵽�׻�̼�����ӵ�ĳһԭ��id������');
                end
            end
            if sum(Cbondatom)==elenum{1,2}+elenum{2,2}*3%ͨ��������ԭ��������ֵ�ͼ����Ƿ�Ϊ�׻�̼
                delatomid(length(delatomid)+1,1)=tarBOinformcopy{k,1};%����ɾ��ԭ�ӵ�id
            end
        end
        row=row-1;
        k=k+1;
    else
        k=k+1;
        row=row-1;
    end
end

[row,~]=size(tarBOinformcopy);%ɾ���׻�̼������
while row
    if ismember(tarBOinformcopy{row,1},delatomid)
        tarBOinformcopy(row,:)=[];%ɾ���׻�̼������
        row=row-1;
    else
        row=row-1;
    end
end


[row,~]=size(tarBOinformcopy);%��������2��ԭ�ӣ�����һ������һ�������ӣ������ʵ�ˮ֮�ǻ���������������
k=1;
while row
    if tarBOinformcopy{k,2}==elenum{3,2} && tarBOinformcopy{k,3}==2%ƥ��Ϊ����ֻ��������ԭ��
        row2=size(tarBOinformcopy,1);len=size(tarBOinformcopy,1);control=0;
        for i=4:5
            for j=1:row2
                if tarBOinformcopy{j,1}==tarBOinformcopy{k,i} && tarBOinformcopy{j,2}==elenum{2,2}%����һ����
                    if i==4%��ʱtarBOinformcopy{k,5}����Ϊ̼
                        for jj=1:row2
                            if tarBOinformcopy{jj,1}==tarBOinformcopy{k,5} && tarBOinformcopy{jj,2}==elenum{1,2} && tarBOinformcopy{jj,3}==4%������̼
                                delatomid(length(delatomid)+1,1)=tarBOinformcopy{k,1};
                                tarBOinformcopy(k,:)=[];
                                control=1;
                                break;
                            end
                        end
                    else %��ʱtarBOinformcopy{k,4}����Ϊ̼
                        for jj=1:row2
                            if tarBOinformcopy{jj,1}==tarBOinformcopy{k,4} && tarBOinformcopy{jj,2}==elenum{1,2} && tarBOinformcopy{jj,3}==4%������̼
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
        if len==size(tarBOinformcopy,1);%�����������1������������2��û��ɾ������������һ��
            k=k+1;
        end
    else
        row=row-1;
        k=k+1;
    end
end


[row,~]=size(tarBOinformcopy);%����׻���
k=1;%����H����ɾ����
while row
    if tarBOinformcopy{k,2}==elenum{2,2}%ƥ��Ϊ��
        delatomid(length(delatomid)+1,1)=tarBOinformcopy{k,1};
        tarBOinformcopy(k,:)=[];%ɾ�����������
        row=row-1;
    else
        k=k+1;
        row=row-1;
    end
end


[row,~]=size(tarBOinformcopy);%�����ʻ�������ɾ����
k=1;
while row
    if tarBOinformcopy{k,2}==elenum{3,2} && tarBOinformcopy{k,3}==1%ƥ��Ϊ����ֻ����һ��ԭ��
        row2=size(tarBOinformcopy,1);len=size(tarBOinformcopy,1);
        for i=1:row2
            if tarBOinformcopy{i,1}==tarBOinformcopy{k,4} && tarBOinformcopy{i,2}==elenum{1,2};%ԭ��idƥ��Ϊ̼
                delatomid(length(delatomid)+1,1)=tarBOinformcopy{k,1};
                tarBOinformcopy(k,:)=[];%ɾ�����ʻ���������
                break;
            end
        end
        row=row-1;
        if len==size(tarBOinformcopy,1);%�����������1������������2�������ʻ���û��ɾ������������һ��
            k=k+1;
        end
    else
        k=k+1;
        row=row-1;
    end
end


[row,~]=size(tarBOinformcopy);%��������3��ԭ�ӣ������ʵ�ˮ����������������
k=1;
while row
    if tarBOinformcopy{k,2}==elenum{3,2} && tarBOinformcopy{k,3}==3%ƥ��Ϊ������������ԭ��
        delatomid(length(delatomid)+1,1)=tarBOinformcopy{k,1};
        tarBOinformcopy(k,:)=[];
        row=row-1;
    else
        row=row-1;
        k=k+1;
    end
end




[row,~]=size(tarBOinformcopy);%��������2��ԭ�ӣ�����һ������һ�������ӣ������ǲ����ʵ��������������������ӣ�����ɾ���������ɻ�
k=1;
while row
    if tarBOinformcopy{k,2}==elenum{3,2} && tarBOinformcopy{k,3}==2%ƥ��Ϊ����ֻ��������ԭ��
        row2=size(tarBOinformcopy,1);control=0;len=size(tarBOinformcopy,1);
        for i=4:5
            for j=1:row2
                if tarBOinformcopy{k,i}==tarBOinformcopy{j,1} && tarBOinformcopy{j,2}==elenum{3,2}%����һ��Ϊ��
                    delatomid(length(delatomid)+1,1)=tarBOinformcopy{k,1};
                    tarBOinformcopy(k,:)=[];
                    control=1;
                    break;
                end
            end
            if control==1%ֻҪ�ҵ���������������ѭ��
                break;
            end
        end
        row=row-1;
        if len==size(tarBOinformcopy,1);%�����������1������������2��û��ɾ������������һ��
            k=k+1;
        end
    
    else
        row=row-1;
        k=k+1;
    end     
end

backboneBO=tarBOinformcopy;

