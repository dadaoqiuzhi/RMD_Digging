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

fprintf('This program is used to read a specified trajectory in *.lammpstrj\n')
dataname=input('\nFilename name of *.lammpstrj file: \n','s');
trajper=input('\nPlease input the output frequency of BO information and trajectory file (Positive integer, see bonds or lammpstrj file):\n');
tartrajectory=input('\nPlease input the timestep of the specified trajectory: \n');
tartrajectory={tartrajectory(1)};
atomnum=atom_num_autoread(dataname);
if mod(tartrajectory{1},trajper)~=0
    control=0;
    fprintf('\nnonexistent trajectory, please check it!!!\n')
    return;
else
    control=1;
end
fprintf('\n强烈建议不要缩放坐标（使用dump_modify scale no)，容易导致计算误差，增大成图的瑕疵')
fprintf('\nScaled coordinate is not recommended (use "dump_modify scale no" to avoid)')
BOXsize=input('\nDoes the coordinate is scaled in the *.lammpstrj file, y/n: \n','s');BOXsize=lower(BOXsize);
if ~ismember(BOXsize,{'y','n'})
    error('Illegal BOXsize parameters, please check it!!!');
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
boxsize=[];
while found
    dataline=fgetl(rawdata);
    readline=readline+1;
    found=found-1;
	if found<=3 && found>=1 
        datacell=textscan(dataline,'%s','delimiter','\n');
        datacellchar=char(datacell{1});
        datarep=strtrim(datacellchar);
        datasplit=strsplit(datarep);
        for i=1:2
            boxsize(4-found,i)=str2num(datasplit{1,i});
        end
    end
end
PBCa=boxsize(1,2)-boxsize(1,1);
PBCb=boxsize(2,2)-boxsize(2,1);
PBCc=boxsize(3,2)-boxsize(3,1);

datacell=textscan(dataline,'%s','delimiter','\n');
datacellchar=char(datacell{1});
datarep=strtrim(datacellchar);
datasplit=strsplit(datarep);
coord_position=[];
if strcmpi('y',BOXsize) && sum(ismember({'ATOMS', 'id','type','xs','ys','zs'},datasplit))==6
    coord_tag={'ATOMS','type','xs','ys','zs'};
else
    coord_tag={'ATOMS','type','x','y','z'};
end
for i=1:5
    if sum(ismember(datasplit,coord_tag(i)))==1
        coord_position(length(coord_position)+1)=find(strcmp(datasplit,coord_tag(i)));
    end
end
if length(coord_position)~=5
    error('Something about atom id，type，x，y，z is lost, please check if the scale answer is right!')
end
if min(coord_position)==coord_position(1)
    coord_position(1)=coord_position(1)-1;
    for j=2:5
        coord_position(j)=coord_position(j)-2;
    end
elseif coord_position(1)>coord_position(2) && coord_position(1)<coord_position(3)
    for j=1:2
        coord_position(j)=coord_position(j)-1;
    end
    for j=3:5
        coord_position(j)=coord_position(j)-2;
    end
elseif max(coord_position)==coord_position(1)
    for j=1:5
        coord_position(j)=coord_position(j)-1;
    end
end

trjdata=[];line=1;
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
	trjdata(line,:)=coord_position_get(coord_position,datasplit);
    if strcmpi('y',BOXsize)
        trjdata(line,3)=boxsize(1,1)+trjdata(line,3)*PBCa;
        trjdata(line,4)=boxsize(2,1)+trjdata(line,4)*PBCb;
        trjdata(line,5)=boxsize(3,1)+trjdata(line,5)*PBCc;
    end
    line=line+1;
end
fclose(rawdata);
disp('\nlammpstrj_analysis is successfully finished.')
fprintf('\nAtomic coordinate information of the specified trajectory is saved in trjdata\n');

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
clear outputans dataoutrow dataoutcol dataoutputrow dataoutcolchar dataoutputcol filename BOXsize boxsize PBCa PBCb PBCc coord_tag coord_position
clear j