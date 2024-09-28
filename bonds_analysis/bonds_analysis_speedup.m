%scrit file name bonds_analysis_speedup
%purpose:
%This program is used to analyze bonds file and is more fast than
%bonds_analysis program
%version 1;2018.6.25
disp('##################################################################################################################################')
disp('Welcome!--by Qiang Liu @Institute of Nuclear Physics and Chemistry, China Academy of Engineering Physics; Email: liubinqiang@163.com');
disp('Repository adress of the Source code on github: https://github.com/dadaoqiuzhi/RMD_Digging');
disp('References: 1.Fuel 287 (2021) 119484. 2.ACS Appl. Mat. Interfaces 13(34) (2021) 41287-41302. 3.ACS Appl. Mat. Interfaces 2022, 14.(4), 5959-5972.')
disp('4.ACS Materials Letters 2023, 2174-2188. More work is coming!')
disp('##################################################################################################################################')
fprintf('This program will read the BO information of a specified or all (not recommended) trajectories, the text in the 2-4 rows of each trajectory is omitted.\n')
dataname=input('Please input the file name to be processed: \n','s');
trajper=input('Please input the output frequency of BO information (Positive integer): \n');
tartrajectory=input('\nPlease input the timestep of the specified trajectory: \n');
atomnum=atom_num_autoread(dataname);
disp('bonds_analysis_speedup is running, please wait...')
tartrajectory={tartrajectory(1)};
if mod(tartrajectory{1},trajper)~=0
    control=0;
    fprintf('\nThis trajectory is not existed, please check it!!!\n')
    return;
else
    control=1;
end


readline=0;
gap=8+atomnum;
rawdata=fopen(dataname,'r');
dataline=fgetl(rawdata);
readline=readline+1;
datacell=textscan(dataline,'%s','delimiter','\n');
datacellchar=char(datacell{1});
datadel=strrep(datacellchar,'#','');
datarep=strtrim(datadel);
datasplit=strsplit(datarep);
if str2num(datasplit{1,2})==tartrajectory{1}
    control=0;
else
    while control
        i=1;
        unfound=1;
        while unfound
            dataline=fgetl(rawdata);
            while isempty(dataline)
                dataline=fgetl(rawdata);
            end
            readline=readline+1;
            i=i+1;
            if i==gap+1
                unfound=0;
                break;
            end
        end
        if mod(readline-1,gap)==0
            datacell=textscan(dataline,'%s','delimiter','\n');
            datacellchar=char(datacell{1});
            datadel=strrep(datacellchar,'#','');
            datarep=strtrim(datadel);
            datasplit=strsplit(datarep);
            if str2num(datasplit{1,2})==tartrajectory{1}
                control=0;
            end
            a0 = ftell(rawdata);
            dataline=fgetl(rawdata);
            a1 = ftell(rawdata);
            dataline=fgetl(rawdata);
            a2 = ftell(rawdata);
            datacell=textscan(dataline,'%s','delimiter','\n');
            datacellchar=char(datacell{1});
            datadel=strrep(datacellchar,'#','');
            datarep=strtrim(datadel);
            datasplit=strsplit(datarep);
            atomnum = str2double(datasplit{length(datasplit)});
            fseek(rawdata,-(a2-a0),'cof');
        else
            disp('Not timestep row, please check it!!!')
            return;
        end
    end
end

found=6;
while found
    dataline=fgetl(rawdata);
    readline=readline+1;
    found=found-1;
end

bondoutdata={};
bondoutdata{1,1}='Timestep';
bondoutdata{1,2}=tartrajectory{1};
for i=3:15
    bondoutdata{1,i}=[];
end

line=2;
while atomnum
    dataline=fgetl(rawdata);
    readline=readline+1;
    atomnum=atomnum-1;
    if atomnum<0
        break;
    end
    datacell=textscan(dataline,'%s','delimiter','\n');
    datacellchar=char(datacell{1});
    datarep=strtrim(datacellchar);
    datasplit=strsplit(datarep);
    bondnumdata={};
    bondnumdata(1,1:3)=datasplit(1,1:3);
    if ~strcmp(datasplit{1,3},'0')
        for i=1:str2num(datasplit{1,3})
            bondnumdata(1,i+3)=datasplit(1,i+3);
        end
        bondnumdata(1,8)=datasplit(1,i+4);
        k=i+5;
        for j=1:str2num(datasplit{1,3})
            bondnumdata(1,j+8)=datasplit(1,k);
            k=k+1;
        end
        bondnumdata(1,13:15)=datasplit(1,k:k+2);
    else
        bondnumdata(1,8)=datasplit(1,4);
        bondnumdata(1,13)=datasplit(1,5);
        bondnumdata(1,14)=datasplit(1,6);
        bondnumdata(1,15)=datasplit(1,7);
    end
    for kk=1:length(bondnumdata)
        if isempty(bondnumdata{kk})
            bondnumdata{kk}='NaN';
        else
            bondnumdata{kk}=str2num(bondnumdata{kk});
        end
    end
    for kk=1:length(bondnumdata)
        bondoutdata{line,kk}=bondnumdata{kk};
    end
    line=line+1;
end
fclose(rawdata);
fprintf('\nbonds_analysis is successfully finished, BO information is saved in bondoutdata')

outputans=input('Export results? y/n?: \n','s');
outputans=lower(outputans);
if outputans=='y'
    [dataoutrow,dataoutcol]=size(bondoutdata);
    dataoutputrow=strcat('A','1');
    dataoutcolchar=char(65+dataoutcol-1);
    dataoutputcol=strcat(dataoutcolchar,num2str(dataoutrow));
    filename='output_mydata.xlsx';
    xlswrite(filename,bondoutdata,dataoutputrow:dataoutputcol)
    fprintf('\nbonds_analysis is successfully finished. BO information is exported into the excel:output_mydata\n')
end
fprintf('\nbonds_analysis is successfully finished. BO information is saved in bondoutdata.\n')

if size(bondoutdata,2) > 15 %delet redundant blank column
    bondoutdata(:,16:end)=[];
end

clear ans atomnum bondnumdata control datacell datacellchar datadel dataline dataname datarep datasplit found gap i j k kk line 
clear outputans rawdata tartrajectory trajper unfound dataoutrow dataoutcol dataoutputrow dataoutcolchar dataoutputcol filename
clear a0 a1 a2
