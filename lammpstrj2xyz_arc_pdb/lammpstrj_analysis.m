%scrit file name lammpstrj_analysis
%purpose:
%This program is used to analysis *.lammpstrj file of REAXC
%version 1;2018.7.16
disp('##################################################################################################################################')
disp('Welcome!--by Qiang Liu @Institute of Nuclear Physics and Chemistry, China Academy of Engineering Physics; Email: liubinqiang@163.com');
disp('Repository adress of the Source code on github: https://github.com/dadaoqiuzhi/RMD_Digging');
disp('References: 1.Fuel 287 (2021) 119484. 2.ACS Appl. Mat. Interfaces 13(34) (2021) 41287-41302. More work is coming!')
disp('##################################################################################################################################')
fprintf('This program is used to read a specified trajectory in *.lammpstrj\n')
dataname=input('\nFilename name of *.lammpstrj file: \n','s');
trajper=input('\nPlease input the output frequency of BO information and trajectory file (Positive integer, see bonds or lammpstrj file):\n');
tartrajectory=input('\nPlease input the timestep of the specified trajectory: \n');
tartrajectory={tartrajectory(1)};
atomnum=input('Please input atom number: \n');
if mod(tartrajectory{1},trajper)~=0
    control=0;
    fprintf('\nnonexistent trajectory, please check it!!!\n')
    return;
else
    control=1;
end
disp('lammpstrj_analysis is running, please wait...')

readline=0;
gap=9+atomnum;
rawdata=fopen(dataname,'r');
dataline=fgetl(rawdata);
readline=readline+1;
dataline=fgetl(rawdata);
readline=readline+1;
datacell=textscan(dataline,'%s','delimiter','\n');
datacellchar=char(datacell{1});
datarep=strtrim(datacellchar);
if str2num(datarep)==tartrajectory{1}
    control=0;
else
    while control
        i=1;
        unfound=1;
        while unfound
            dataline=fgetl(rawdata);
            readline=readline+1;
            i=i+1;
            if i==gap+1
                unfound=0;
                break;
            end
        end
        if mod(readline-2,gap)==0
            datacell=textscan(dataline,'%s','delimiter','\n');
            datacellchar=char(datacell{1});
            datarep=strtrim(datacellchar);
            if str2num(datarep)==tartrajectory{1}
                control=0;
            end
        else
            disp('Not a timestep line, please check it!!!')
            return;
        end
    end
end

found=7;
while found
    dataline=fgetl(rawdata);
    readline=readline+1;
    found=found-1;
end

trjdata=[];line=1;
while atomnum
    dataline=fgetl(rawdata);
    readline=readline+1;
    atomnum=atomnum-1;
    if atomnum<=0
        break;
    end
    datacell=textscan(dataline,'%s','delimiter','\n');
    datacellchar=char(datacell{1});
    datarep=strtrim(datacellchar);
    datasplit=strsplit(datarep);
    for i=1:length(datasplit)
        trjdata(line,i)=str2num(datasplit{1,i});
    end
    line=line+1;
end
fclose(rawdata);
disp('\nlammpstrj_analysis is successfully finished.')
fprintf('\nAtomic coordination information of the specified trajectory is saved in trjdata\n');

outputans=input('\nExport data to Excel? Much time is required for large data and the Excel should be closed.y/n: \n','s');
outputans=lower(outputans);
if outputans=='y'
    [dataoutrow,dataoutcol]=size(trjdata);
    dataoutputrow=strcat('A','1');
    dataoutcolchar=char(65+dataoutcol-1);
    dataoutputcol=strcat(dataoutcolchar,num2str(dataoutrow));
    filename='output_mydata.xlsx';
    xlswrite(filename,trjdata,dataoutputrow:dataoutputcol)
    fprintf('\nlammpstrj_analysis is successfully finished, BO information etc. are exported to Excel:output_mydata.\n')
end
clear atomnum control datacell datacellchar dataline dataname datarep datasplit found gap i line rawdata tartrajectory trajper unfound ans
clear outputans dataoutrow dataoutcol dataoutputrow dataoutcolchar dataoutputcol filename 