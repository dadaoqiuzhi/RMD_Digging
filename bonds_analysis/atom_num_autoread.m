%scrit file name atom_num_autoread
%purpose:
%This program is used to read atom number in the simulation box
%automatically
%version 1;2024.04.04

function atomnum=atom_num_autoread(dataname)
fprintf('\nRead atom number from bond.* file, please wait...\n')
rawdata=fopen(dataname,'r');
for i=1:3
    dataline=fgetl(rawdata);
end
datacell=textscan(dataline,'%s','delimiter','\n');
datacellchar=char(datacell{1});
datarep=strtrim(datacellchar);
if contains(datarep,'# Number of particles')
    fprintf('\nThe atom number line is found, ready to read it\n')
end
datadelimiter={'# Number of particles'};
[C,~]=strsplit(datarep,datadelimiter,'CollapseDelimiters',false);
atomnum=str2num(C{2});
fprintf('The atom number is found, which is%d\n',atomnum)
fclose(rawdata);%
end
