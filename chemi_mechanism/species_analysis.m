%scrit file name species_analysis
%purpose:
%This program is used to analysis species file
%version 1;2018.6.21
disp('Welcome!--by Qiang Liu @Institute of Nuclear Physics and Chemistry, China Academy of Engineering Physics; Email: liubinqiang@163.com');
disp('Repository adress of the Source code on github: https://github.com/dadaoqiuzhi/RMD_Digging');

fprintf('\nspecies_analysis is running, Please wait...')
outputdata={};line=1;
rawdata=fopen(datanamespe,'r');
dataline=fgetl(rawdata);
datacell=textscan(dataline,'%s','delimiter','\n');
datacellchar=char(datacell{1});
datadel=strrep(datacellchar,'#','');
datarep=strtrim(datadel);
datasplit=strsplit(datarep);
datacellnum=length(datasplit);
indexapp=[];
for j=1:datacellnum
    outputdata{1,j}=datasplit{j};
    indexapp(j)=j;
end
datalinenum=fgetl(rawdata);
datalinenum=strtrim(datalinenum);
datalinenum=strread(datalinenum);
for i=1:length(indexapp)
    outputdata{2,i}=datalinenum(i);
end


line=3;
while ~feof(rawdata)
    dataline=fgetl(rawdata);
    if ~isempty(dataline)
        datacell=textscan(dataline,'%s','delimiter','\n');
        datacellchar=char(datacell{1});
        datadel=strrep(datacellchar,'#','');
        datarep=strtrim(datadel);
        datasplit=strsplit(datarep);
        datacellnum=length(datasplit);
        datafirstrow=outputdata(1,:);
        [indexcol,indexovlp,indexapp]=membercheck(datasplit,datafirstrow);
        [~,datacol]=size(outputdata);
        for k=1:length(indexapp)
            outputdata{1,k+datacol}=datasplit{indexapp(k)};
        end
        
        datalinenum=fgetl(rawdata);
        datalinenum=strtrim(datalinenum);
        datalinenum=strread(datalinenum);
        for i=1:length(indexcol)
            outputdata{line,indexcol(i)}=datalinenum(indexovlp(i));
        end
        for i=1:length(indexapp)
            outputdata{line,i+datacol}=datalinenum(indexapp(i));
        end
    end
    line=line+1;
end
fclose(rawdata);

[datarow,datacol]=size(outputdata);
for i=2:datarow
    for j=1:datacol
        if isempty(outputdata{i,j})
        outputdata{i,j}=0;
        end
    end
end

fprintf('\nSpecies is successfully sorted and imported in outputdata!')
fprintf('species_analysis is successfully finished.')

clear datacell datacellchar datacellnum datacol datadel datafirstrow dataline datalinenum datanow filename ans
clear dataoutcol dataoutcolchar dataoutputcol dataoutputrow datarow datarep datarow datasec datasplit statans
clear i j k line rawdata charcor datadelimiter dataoutrow matches outputans datanamespe indexapp indexcol indexovlp

