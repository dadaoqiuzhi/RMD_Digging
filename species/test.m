disp('��ӭʹ�ñ�����--by ��ǿ@�Ĵ���ѧ�߷��ӿ�ѧ�빤��ѧԺ����ܽ��ڿ����飬liubinqiang@163.com')
disp('������species_analysis����󣬿��ô˳��������Ҫ����Ĳ��ֲ��ﺬ����ʱ��仯�Ľ��,�������ļ�ǰ���зǷ���ʽ��Ϣ')
species=input('������������ʽ������ʽӦ���ļ���һ�£���Ҫע��Ԫ��˳�򣩣�������ʽ���ÿո�ֿ����Ҵ����ţ�\n');
species=upper(species);
fprintf('\n\nspecies_capture����������,��ȴ�...\n\n')
species=strtrim(species);
species=strsplit(species);
outputdatast=outputdata(1,:);
datamatch=[];
for i=1:length(species)
    for j=1:length(outputdatast)
        if strcmp(species{i},upper(outputdatast{j}))
            datamatch(1,i)=i;
            datamatch(2,i)=j;
        end
    end
end

[checkrow,checkcol]=size(datamatch);%���ƥ�����
if checkcol~=length(species)
    if isempty(datamatch)
        fprintf('��ȫû��ƥ�䵽���ԭʼ�ļ�����δ�иò���������������С�ļ�飡����');
        return;
    else
        fprintf('���ֲ���û��ƥ�䵽��ԭʼ�ļ�����δ�иò���������������С�ļ�飡����');
        return;
    end
end
if sum(ismember(datamatch,0))>=1
    fprintf('�������ƥ����󣬵���С������Ϊ0��������������ʽ�Ų飬��չ����ռ�����������ļ�������');
    return;
end

    
outputdatanew={};
for k=1:3
    outputdatanew(:,k)=outputdata(:,k);
end
for j=1:length(species)
    outputdatanew(:,j+3)=outputdata(:,datamatch(2,j));
end
disp('species_capture�����������outputdatanew��')

expoans=input('�����Ƿ񵼳����ݵ�excel��y������n�������������š�y/n:\n');
expoans=lower(expoans);
if expoans=='y'
    [dataoutrow,dataoutcol]=size(outputdatanew);%��������
    dataoutputrow=strcat('A','1');
    dataoutcolchar=char(65+dataoutcol-1);
    dataoutputcol=strcat(dataoutcolchar,num2str(dataoutrow));
    filename='output_mydata.xlsx';
    xlswrite(filename,outputdatanew,dataoutputrow:dataoutputcol)
    disp('species_capture��������Ѿ�������excel:output_mydata��')
end
fprintf('\n\nspecies_capture�������н���\n\n')