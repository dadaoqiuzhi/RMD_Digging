%scrit file name loglammps
%purpose:
%This program is used to analysis log.lammps file
%version 1;2018.6.20
disp('##################################################################################################################################')
disp('Welcome!--by Qiang Liu @Institute of Nuclear Physics and Chemistry, China Academy of Engineering Physics; Email: liubinqiang@163.com');
disp('Repository adress of the Source code on github: https://github.com/dadaoqiuzhi/RMD_Digging');
disp('References: 1.Fuel 287 (2021) 119484. 2.ACS Appl. Mat. Interfaces 13(34) (2021) 41287-41302. 3.ACS Appl. Mat. Interfaces 2022, 14.(4), 5959-5972.')
disp('4.ACS Materials Letters 2023, 2174-2188. More work is coming!')
disp('##################################################################################################################################')

tic;
dataname = input('Input file name should be processed: \n','s');
coloum_num = input('Please input the total column number of data expected to be exported, deleting text in advance is not necessary anymore: \n');
fprintf('loglammps is running, please wait...\n')

rawdata = [];
fid = fopen(dataname,'r');
while ~feof(fid)
    dataline = fgetl(fid);
    while isempty(dataline) && ~feof(fid) 
        dataline = fgetl(fid);
    end
    datacell = textscan(dataline,'%s','delimiter','\n');
    datacellchar = char(datacell{1});
    datarep = strtrim(datacellchar);
    datasplit = strsplit(datarep);
    data = num_obtain(datasplit,coloum_num);
    if ismatrix(data) && ~ischar(data)
        if isempty(rawdata)
            rawdata(size(rawdata,1)+1,:) = data;
        else
            if length(data) == size(rawdata,2)
                rawdata = del_repeat(rawdata,data);
                rawdata(size(rawdata,1)+1,:) = data;
            end
        end
    end
end
fclose(fid);

datacolumn=input('Column Number of the data that should be processed, multi numbers should be seperated by whitespace; "all" for all data: \n','s');
datainput=[];
if strcmpi(datacolumn,'all')==1
    datainput=rawdata;
else
    datacol=str2num(datacolumn);
    for i=1:length(datacol)
        datainput(:,i)=rawdata(:,datacol(i));
    end
end




timestep=input('Input timestep value (fs), whose unit will be converted to be ps: \n');
avesteps=input('Step/Span number used to average the data, should be a positive integer. these data can be also averagely treated by statiave code (recommend!): \n');
dataoutput=[];thermoper=datainput(2,1)-datainput(1,1);
[datarow,datacol]=size(datainput);
imax=(datainput(datarow,1))/thermoper+1;re=mod(imax,avesteps);
imaxend=(imax-re)/avesteps;
for i=1:imaxend
    if (i-1)*avesteps+1<=size(datainput,1)
        dataoutput(i,1)=datainput((i-1)*avesteps+1,1);
    end
end
for j=2:datacol
    for i=1:imaxend
        if avesteps*i<=size(datainput,1)
            dataoutput(i,j)= mean(datainput((i-1)*avesteps+1:avesteps*i,j));
        end
    end
end
dataoutput(:,1)=dataoutput(:,1)*timestep/1000;

Elapsedtime = toc;
fprintf('loglammps is end. Total task time: %.2f s\n',Elapsedtime)
disp('The averaged data is saved indataoutput')

% [dataoutrow,dataoutcol]=size(dataoutput);
% dataoutputrow=strcat('A','1');
% dataoutcolchar=char(65+dataoutcol-1);
% dataoutputcol=strcat(dataoutcolchar,num2str(dataoutrow));
% filename='output_mydata.xlsx';
% xlswrite(filename,dataoutput,dataoutputrow:dataoutputcol)

clear datacol datacolumn dataname datarow i j rawdata re thermoper timestep avesteps coloum_num data datacell datacellchar
clear dataline datarep datasplit fid imax imaxend ans Elapsedtime
