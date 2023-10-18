%scrit file name bonds_analysis
%purpose:
%This program is used to analyze bonds file
%version 1;2018.6.25
disp('##################################################################################################################################')
disp('Welcome!--by Qiang Liu @Institute of Nuclear Physics and Chemistry, China Academy of Engineering Physics; Email: liubinqiang@163.com');
disp('Repository adress of the Source code on github: https://github.com/dadaoqiuzhi/RMD_Digging');
disp('References: 1.Fuel 287 (2021) 119484. 2.ACS Appl. Mat. Interfaces 13(34) (2021) 41287-41302. 3.ACS Appl. Mat. Interfaces 2022, 14.(4), 5959-5972.')
disp('4.ACS Materials Letters 2023, 2174-2188. More work is coming!')
disp('##################################################################################################################################')
fprintf('This program will read the BO information of a specified or all (not recommended) trajectories, the text in the 2-4 rows of each trajectory is omitted.\n')
dataname=input('Please input the file name to be processed: \n','s');
readmethod=input('Please select the way to read BO information: 1.specified trajectory 2.all trajectories. Select the No.: \n');
trajper=input('Please input the output frequency of BO information (Positive integer): \n');


if readmethod==1
    tartrajectory=input('\nPlease input the timestep of the specified trajectory: \n');
    disp('species_analysis is running, please wait...')
    tartrajectory={tartrajectory(1)};
   
    if mod(tartrajectory{1},trajper)~=0
        control=0;
        fprintf('\nThis trajectory is not existed, please check it!!!\n')
    else
        control=5;
    end
    bondoutdata={};line=1;
    readline=0;
    trajectorynum=0;
    rawdata=fopen(dataname,'r');
    while control
        if ~feof(rawdata)
            dataline=fgetl(rawdata);
            readline=readline+1;
            if ~isempty(dataline)
                datacell=textscan(dataline,'%s','delimiter','\n');
                datacellchar=char(datacell{1});
                datadel=strrep(datacellchar,'#','');
                datarep=strtrim(datadel);
                if ~isempty(datarep)
                    datasplit=strsplit(datarep);
                    datadelimiter={'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};
                    datasplitupper=upper(datasplit{1});
                    [C,matches]=strsplit(datasplitupper,datadelimiter,'CollapseDelimiters',false);
                    
                    
                    if isempty(C{1})
                        if strcmp(datasplit{1},'Timestep') 
                            if str2num(datasplit{2})==tartrajectory{1}
                                bondoutdata{line,1}=datasplit{1,1};
                                bondoutdata{line,2}=str2num(datasplit{1,2});
                                line=line+1;
                                trajectorynum=trajectorynum+1;

                                while control
                                    if ~feof(rawdata)
                                        dataline=fgetl(rawdata);
                                        readline=readline+1;
                                        if ~isempty(dataline)
                                            datacell=textscan(dataline,'%s','delimiter','\n');
                                            datacellchar=char(datacell{1});
                                            datadel=strrep(datacellchar,'#','');
                                            datarep=strtrim(datadel);
                                            if ~isempty(datarep)
                                                datasplit=strsplit(datarep);
                                                datasplitupper=upper(datasplit{1});
                                                [C,matches]=strsplit(datasplitupper,datadelimiter,'CollapseDelimiters',false);
                                                if isempty(matches)
                                                    bondnumdata={};
                                                    bondnumdata(1,1:3)=datasplit(1,1:3);
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
                                                else
                                                    control=control-1;
                                                end
                                            end
                                        end
                                    else
                                        control=0;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end




if readmethod==2
    disp('species_analysis is running, please wait...')
    bondoutdata={};line=1;
    trajectorynum=0;
    rawdata=fopen(dataname,'r');
    while ~feof(rawdata)
        dataline=fgetl(rawdata);
        if ~isempty(dataline)
            datacell=textscan(dataline,'%s','delimiter','\n');
            datacellchar=char(datacell{1});
            datadel=strrep(datacellchar,'#','');
            datarep=strtrim(datadel);
            if ~isempty(datarep)
                datasplit=strsplit(datarep);
                datadelimiter={'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};
                datasplitupper=upper(datasplit{1});
                [C,matches]=strsplit(datasplitupper,datadelimiter,'CollapseDelimiters',false);
                
                
                if isempty(C{1})
                    if strcmp(datasplit{1},'Timestep')
                        bondoutdata{line,1}=datasplit{1,1};
                        bondoutdata{line,2}=str2num(datasplit{1,2});
                        line=line+1;
                        trajectorynum=trajectorynum+1;
                    end
                end
                
                if isempty(matches)
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
            end
        end
    end
end
line=line-1;
fclose(rawdata);
fprintf('\nbonds_analysis is successfully finished, BO information is saved inbondoutdata')

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
fprintf('\nbonds_analysis is successfully finished. BO information is saved inbondoutdata.\n')

clear bondnumdata C datacell datacellchar datadel datadelimiter dataline dataname datarep datasplit datasplitupper
clear i j k kk matches rawdata outputans dataoutrow dataoutcol dataoutcolchar dataoutputrow dataoutputcol filename
clear tartrajectory readmethod datainsplit control trajper
