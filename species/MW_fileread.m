%scrit file name MW_fileread
%purpose:
%This program is used to analysis species file
%version 1;2023.8.27

outputdata_temp={};
dataline=fgetl(rawdata);
datacell=textscan(dataline,'%s','delimiter','\n');
datacellchar=char(datacell{1});
datadel=strrep(datacellchar,'#','');
datarep=strtrim(datadel);
datasplit=strsplit(datarep);
datacellnum=length(datasplit);
indexapp=[];
for j=1:datacellnum
    outputdata_temp{1,j}=datasplit{j};
    indexapp(j)=j;
end
datalinenum=fgetl(rawdata);
datalinenum=strtrim(datalinenum);
datalinenum=strread(datalinenum);
for i=1:length(indexapp)
    outputdata_temp{2,i}=datalinenum(i);
end

