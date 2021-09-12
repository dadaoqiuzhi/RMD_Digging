disp('欢迎使用本程序--by 刘强@四川大学高分子科学与工程学院李光宪教授课题组，liubinqiang@163.com')
disp('运行了species_analysis程序后，可用此程序输出想要输出的部分产物含量随时间变化的结果,保留了文件前三列非分子式信息')
species=input('请输入产物分子式，分子式应与文件中一致（主要注意元素顺序），各分子式间用空格分开，且带引号：\n');
species=upper(species);
fprintf('\n\nspecies_capture程序运行中,请等待...\n\n')
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

[checkrow,checkcol]=size(datamatch);%检查匹配情况
if checkcol~=length(species)
    if isempty(datamatch)
        fprintf('完全没有匹配到产物，原始文件可能未有该产物或者输入错误，请小心检查！！！');
        return;
    else
        fprintf('部分产物没有匹配到，原始文件可能未有该产物或者输入错误，请小心检查！！！');
        return;
    end
end
if sum(ismember(datamatch,0))>=1
    fprintf('产物存在匹配错误，导致小标索引为0，请逐个输入分子式排查，清空工作空间产生的无用文件！！！');
    return;
end

    
outputdatanew={};
for k=1:3
    outputdatanew(:,k)=outputdata(:,k);
end
for j=1:length(species)
    outputdatanew(:,j+3)=outputdata(:,datamatch(2,j));
end
disp('species_capture分析结果存在outputdatanew中')

expoans=input('请问是否导出数据到excel表？y导出，n不导出，带引号。y/n:\n');
expoans=lower(expoans);
if expoans=='y'
    [dataoutrow,dataoutcol]=size(outputdatanew);%导出数据
    dataoutputrow=strcat('A','1');
    dataoutcolchar=char(65+dataoutcol-1);
    dataoutputcol=strcat(dataoutcolchar,num2str(dataoutrow));
    filename='output_mydata.xlsx';
    xlswrite(filename,outputdatanew,dataoutputrow:dataoutputcol)
    disp('species_capture分析结果已经导出到excel:output_mydata中')
end
fprintf('\n\nspecies_capture程序运行结束\n\n')