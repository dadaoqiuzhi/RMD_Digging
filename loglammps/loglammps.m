%scrit file name loglammps
%purpose:
%This program is used to analysis the log.lammps file
%version 1;2018.6.20
disp('##################################################################################################################################')
disp('Welcome!--by Qiang Liu @Institute of Nuclear Physics and Chemistry, China Academy of Engineering Physics; Email: liubinqiang@163.com');
disp('Repository adress of the Source code on github: https://github.com/dadaoqiuzhi/RMD_Digging');
disp('References: 1.Fuel 287 (2021) 119484. 2.ACS Appl. Mat. Interfaces 13(34) (2021) 41287-41302. More work is coming!')
disp('##################################################################################################################################')
disp('load importdata, dlmread functions for pure numerical data; textread and textscan functions for text data, fopen for complex cases');
disp('Please delete the text data in the beginning and end of the log.file, and only numerical data are left.')

dataname=input('Input file name should be processed: \n','s');
rawdata=load(dataname,'-ascii');
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

timestep=input('Input timestep value (fs), whose unit will be converted to be ps: \n');%average value will be calculated
avesteps=input('Step/Span number used to average the data, should be a positive integer. these data can be also averagely treated by statiave code (recommend!): \n');
dataoutput=[];thermoper=datainput(2,1)-datainput(1,1);
[datarow,datacol]=size(datainput);
imax=(datainput(datarow,1))/thermoper+1;re=mod(imax,avesteps);
imaxend=(imax-re)/avesteps;
for i=1:imaxend
    dataoutput(i,1)=datainput((i-1)*avesteps+1,1);
end
for j=2:datacol
    for i=1:imaxend
       dataoutput(i,j)= mean(datainput((i-1)*avesteps+1:avesteps*i,j));
    end
end

dataoutput(:,1)=dataoutput(:,1)*timestep;
% [dataoutrow,dataoutcol]=size(dataoutput);%export data
% dataoutputrow=strcat('A','1');
% dataoutcolchar=char(65+dataoutcol-1);
% dataoutputcol=strcat(dataoutcolchar,num2str(dataoutrow));
% filename='output_mydata.xlsx';
% xlswrite(filename,dataoutput,dataoutputrow:dataoutputcol)

clear datacol datacolumn dataname datarow i j rawdata re thermoper timestep avesteps

disp('The averaged data are stored in dataoutput matrix, can be exported to excel (output_mydata) by modifying the code block which is commented out')