%scrit file name atom_num_autoreadtrj
%purpose:
%This program is used to read atom number in the simulation box
%automatically
%version 1;2023.09.05

function atomnum=atom_num_autoreadtrj(datanametrj)
rawdatatrj=fopen(datanametrj,'r');
for i=1:3
    dataline=fgetl(rawdatatrj);
end
datacell=textscan(dataline,'%s','delimiter','\n');
datacellchar=char(datacell{1});
datarep=strtrim(datacellchar);
if strcmpi(datarep,'ITEM: NUMBER OF ATOMS')
    fprintf('\nThe atom number line is found, ready to read it\n')
end
dataline=fgetl(rawdatatrj);
datacell=textscan(dataline,'%s','delimiter','\n');
datacellchar=char(datacell{1});
datarep=strtrim(datacellchar);
datarep=str2num(datarep);
if datarep==round(datarep) && datarep>0
    fprintf('The atom number is found, which is %d\n',datarep)
    atomnum=datarep;
end
fclose(rawdatatrj);
end
