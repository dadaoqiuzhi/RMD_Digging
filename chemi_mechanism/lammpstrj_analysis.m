%scrit file name lammpstrj_analysis
%purpose:
%This program is used to analysis *.lammpstrj file of REAXC
%version 1;2018.7.16
disp('##################################################################################################################################')
disp('Welcome!--by Qiang Liu @Institute of Nuclear Physics and Chemistry, China Academy of Engineering Physics; Email: liubinqiang@163.com');
disp('Repository adress of the Source code on github: https://github.com/dadaoqiuzhi/RMD_Digging');
disp('References: 1.Fuel 287 (2021) 119484. 2.ACS Appl. Mat. Interfaces 13(34) (2021) 41287-41302. More work is coming!')
disp('##################################################################################################################################')
if mod(tartrajectory{1,1},trajper)~=0
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
        error('Timestep of reactants is exceed the first trajectory in lammpstrj, which is nonexixtent. Please check it!!!');
    end
    fprintf('\nMinimum timestep in dump.*%d > target timestep %d, mainly due to the beginning timestep in dump.* (species) is not 0. Skip and continue to search the next one\n',str2num(datarep),tartrajectory{1});
    warndlg('Warning: dMinimum timestep in dump.* > target timestep!Skip and continue');
    return;
elseif str2num(datarep)<tartrajectory{1}
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
        if mod(readline-2,gap)==0
            datacell=textscan(dataline,'%s','delimiter','\n');
            datacellchar=char(datacell{1});
            datarep=strtrim(datacellchar);
            if str2num(datarep)==tartrajectory{1}
                control=0;
            end
        else
            disp('This is not a line with timestep, please check it!!!')
            return;
        end
    end
end

dataline=textscan(rawdatatrj,'%q',1,'headerlines',2,'delimiter','\n');
boxsize=[];
dataline=fgetl(rawdatatrj);
datacell=textscan(dataline,'%s','delimiter','\n');
datacellchar=char(datacell{1});
datarep=strtrim(datacellchar);
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
    for i=1:length(datasplit)
        trjdata(line,i)=str2num(datasplit{1,i});
    end
    line=line+1;
end
fprintf('\nlammpstrj_analysis is successfully finished.')
fprintf('\nBO of the target trajecrory is saved in trjdata.');


clear atomnumcopy control datacell datacellchar dataline datasplit found gap i line unfound ans
clear outputans dataoutrow dataoutcol dataoutputrow dataoutcolchar dataoutputcol 