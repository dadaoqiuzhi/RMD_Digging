%scrit file name lammpstrj_analysis
%purpose:
%This program is used to analysis *.lammpstrj file of REAXC
%version 1;2018.7.16
disp('##################################################################################################################################')
disp('Welcome!--by Qiang Liu @Institute of Nuclear Physics and Chemistry, China Academy of Engineering Physics; Email: liubinqiang@163.com');
disp('Repository adress of the Source code on github: https://github.com/dadaoqiuzhi/RMD_Digging');
disp('References: 1.Fuel 287 (2021) 119484. 2.ACS Appl. Mat. Interfaces 13(34) (2021) 41287-41302. 3.ACS Appl. Mat. Interfaces 2022, 14.(4), 5959-5972.')
disp('4.ACS Materials Letters 2023, 2174-2188. More work is coming!')
disp('##################################################################################################################################')
fprintf('This program is used to read trajectories from *.lammpstrj file\n')

if mod(tartrajectory{1},trajper)~=0%
    control=0;
    fprintf('\nNonexistent trajectory, please check it!\n')
    return;
else
    control=1;
end
disp('lammpstrj_analysis is running, please wait...')

readline=0;
gap=9+atomnum;
rawdatatrj=fopen(datanametrj,'r');
dataline=fgetl(rawdatatrj);
readline=readline+1;
dataline=fgetl(rawdatatrj);
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
            dataline=fgetl(rawdatatrj);
            readline=readline+1;
            i=i+1;
            if i==gap+1
                unfound=0;
                break;
            end
        end
        if mod(readline-2,gap)==0%
            datacell=textscan(dataline,'%s','delimiter','\n');
            datacellchar=char(datacell{1});
            datarep=strtrim(datacellchar);
            if str2num(datarep)==tartrajectory{1}
                control=0;
            end
        else
            disp('This is not a timestep line, please check the file or code')
            return;
        end
    end
end

found=7;
while found
    dataline=fgetl(rawdatatrj);
    readline=readline+1;
    found=found-1;
end

trjdata=[];line=1;atomnumcopy=atomnum;
while atomnumcopy
    dataline=fgetl(rawdatatrj);
    readline=readline+1;
    atomnumcopy=atomnumcopy-1;
    datacell=textscan(dataline,'%s','delimiter','\n');
    datacellchar=char(datacell{1});
    datarep=strtrim(datacellchar);
    datasplit=strsplit(datarep);
    for i=1:length(datasplit)
        trjdata(line,i)=str2num(datasplit{1,i});
    end
    line=line+1;
end
fclose(rawdatatrj);
fprintf('\nlammpstrj_analysis is end')
fprintf('\nThe target atomic coordinate information of the interested trajectory is saved in trjdata\n');


% clear atomnum control datacell datacellchar dataline dataname datarep datasplit found gap i line rawdatatrj tartrajectory trajper unfound ans
% clear outputans dataoutrow dataoutcol dataoutputrow dataoutcolchar dataoutputcol filename 