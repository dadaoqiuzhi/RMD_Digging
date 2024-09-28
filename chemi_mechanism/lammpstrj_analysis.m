%scrit file name lammpstrj_analysis
%purpose:
%This program is used to analysis *.lammpstrj file of REAXC
%version 1;2018.7.16
%fprintf('本程序用于读取*.lammpstrj文件中某一轨迹\n')
disp('##################################################################################################################################')
disp('Welcome!--by Qiang Liu @Institute of Nuclear Physics and Chemistry, China Academy of Engineering Physics; Email: liubinqiang@163.com');
disp('References: 1.Fuel 287 (2021) 119484. 2.ACS Appl. Mat. Interfaces 13(34) (2021) 41287-41302. 3.ACS Appl. Mat. Interfaces 2022, 14.(4), 5959-5972.')
disp('4.ACS Materials Letters 2023, 2174-2188. More work is coming!')
disp('##################################################################################################################################')

if ~ismember(tartrajectory{1,1},cell2mat(outputdatanew(2:end,1)))
    control=0;
    fprintf('\nNonexisent trajectory, please check it!!!\n')
    return;
else
    control=1;
end
fprintf('\nlammpstrj_analysis is running, please wait...')

readline=0;
gap=9+atomnum;

dataline=fgetl(rawdatatrj);
readline=readline+1;
dataline=fgetl(rawdatatrj);
readline=readline+1;
datacell=textscan(dataline,'%s','delimiter','\n');
datacellchar=char(datacell{1});
datarep=strtrim(datacellchar);
if str2num(datarep)==tartrajectory{1}
    control=0;
elseif str2num(datarep)>tartrajectory{1}
    if choi==2 || choi==4
        fprintf('\nMinimum timestep in lammpstrj %d is >= target timestep%d, mainly due to the beginning timestep in lammpstrj is not 0,\nplease check it\n',str2num(datarep),tartrajectory{1});
        error('Timestep exceed the first frame No. in lammpstrj.* file, which is nonexixtent. Please check it!!!');
    end
    fprintf('nMinimum timestep in dump.*%d > target timestep %d, mainly due to the beginning timestep in dump.* (species) is not 0. Skip and continue to search the next one\n',str2num(datarep),tartrajectory{1});
    warndlg('Warning: dMinimum timestep in dump.* > target timestep!Skip and continue');
    return;
elseif str2num(datarep)<tartrajectory{1}
    while control
        i=1;
        unfound=1;
        while unfound
            dataline=fgetl(rawdatatrj);
            while isempty(dataline)%可能存在ITEM: TIMESTEP行前空行
                dataline=fgetl(rawdatatrj);
            end
            readline=readline+1;
            i=i+1;
            if i==gap+1
                unfound=0;
                break;
            end
        end
        if mod(readline-2,gap)==0
            if feof(rawdatatrj)
                fprintf('\nThe last line of lammpstrj.* file is reached. Please check if the atom number or timestep is correct.\n The timestep in species.* and lammpstrj.* maybe inconsistent, eg. 29,429...\n')
                error('As above, please check it!')
            end
            datacell=textscan(dataline,'%s','delimiter','\n');
            datacellchar=char(datacell{1});
            datarep=strtrim(datacellchar);
            if str2num(datarep)==tartrajectory{1}
                control=0;
            elseif str2num(datarep)>tartrajectory{1}
                error('The frame No. appointed according to species.* file exceeds that in lammpstrj.* file, please check it!')
            end
        else
            disp('This is not a line with timestep information, please check it!!!')
            return;
        end
    end
end
%
dataline=textscan(rawdatatrj,'%q',1,'headerlines',2,'delimiter','\n');
boxsize=[];
dataline=fgetl(rawdatatrj);
datacell=textscan(dataline,'%s','delimiter','\n');
datacellchar=char(datacell{1});
datarep=strtrim(datacellchar);%
datasplit=strsplit(datarep);
for i=1:length(datasplit)
    boxsize(1,i)=str2num(datasplit{1,i});
end
dataline=fgetl(rawdatatrj);
datacell=textscan(dataline,'%s','delimiter','\n');
datacellchar=char(datacell{1});
datarep=strtrim(datacellchar);
datasplit=strsplit(datarep);
for i=1:length(datasplit)
    boxsize(2,i)=str2num(datasplit{1,i});
end
dataline=fgetl(rawdatatrj);
datacell=textscan(dataline,'%s','delimiter','\n');
datacellchar=char(datacell{1});
datarep=strtrim(datacellchar);
datasplit=strsplit(datarep);
for i=1:length(datasplit)
    boxsize(3,i)=str2num(datasplit{1,i});
end
if formatout==2 || formatout==3
    %
    if strcmp(PBCchoi,'ON')
        PBC='PBC=ON';
        PBCa=boxsize(1,2)-boxsize(1,1);
        PBCb=boxsize(2,2)-boxsize(2,1);
        PBCc=boxsize(3,2)-boxsize(3,1);
    elseif strcmp(PBCchoi,'OFF')
        if formatout==2
            PBC='PBC=OFF';
        elseif formatout==3
            PBCa=0.00;PBCb=PBCa;PBCc=PBCa;
            PBCalpha=0.00;PBCbeta=PBCalpha;PBCgamma=PBCalpha;
        end
    else
        disp('Illegal periodic boundary condition in PBCchoi, please check it!!!');
        return;
    end
end

dataline=textscan(rawdatatrj,'%q',1,'headerlines',0,'delimiter','\n');
datacell=dataline;
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

trjdata=[];line=1;atomnumcopy=atomnum;
while atomnumcopy
    dataline=fgetl(rawdatatrj);
    readline=readline+1;
    atomnumcopy=atomnumcopy-1;
    if atomnumcopy<=0
        break;
    end
    datacell=textscan(dataline,'%s','delimiter','\n');
    datacellchar=char(datacell{1});
    datarep=strtrim(datacellchar);
    datasplit=strsplit(datarep);
    trjdata(line,:)=coord_position_get(coord_position,datasplit);
    line=line+1;
end
fprintf('\nlammpstrj_analysis is successfully finished')
fprintf('\nAtom coordinate of the target trajecrory is saved in trjdata');


clear atomnumcopy control datacell datacellchar dataline datasplit found gap i line unfound ans
clear outputans dataoutrow dataoutcol dataoutputrow dataoutcolchar dataoutputcol 