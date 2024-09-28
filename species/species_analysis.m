%scrit file name species_analysis
%purpose:
%This program is used to analysis species file
%version 1;2018.6.21
disp('##################################################################################################################################')
disp('Welcome!--by Qiang Liu @Institute of Nuclear Physics and Chemistry, China Academy of Engineering Physics; Email: liubinqiang@163.com');
disp('Repository adress of the Source code on github: https://github.com/dadaoqiuzhi/RMD_Digging');
disp('References: 1.Fuel 287 (2021) 119484. 2.ACS Appl. Mat. Interfaces 13(34) (2021) 41287-41302. More work is coming!')
disp('##################################################################################################################################')
disp('This procedure sorts the species evolution with time according to species file')

dataname=input('Please input the file name to be processed: \n','s');
tic
disp('species_analysis is running, please wait...')
outputdata={};line=1;
rawdata=fopen(dataname,'r');
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
    if ~isempty(dataline) && ischar(dataline) && length(dataline) > 1
        datacell=textscan(dataline,'%s','delimiter','\n');
        datacellchar=char(datacell{1});
        datadel=strrep(datacellchar,'#','');
        datarep=strtrim(datadel);
        datasplit=strsplit(datarep);
        datacellnum=length(datasplit);
        datafirstrow=outputdata(1,:);
        [indexcol,indexovlp,indexapp]=membercheck(datasplit,datafirstrow);
        
        [datarow,datacol]=size(outputdata);
        for k=1:length(indexapp)
            outputdata{1,k+datacol}=datasplit{indexapp(k)};
        end
        
        datalinenum=fgetl(rawdata);
        datalinenum=strtrim(datalinenum);
        datalinenum=strread(datalinenum);
        %datalinenum=strsplit(datalinenum);
        for i=1:length(indexcol)
            outputdata{line,indexcol(i)}=datalinenum(indexovlp(i));
        end
        for i=1:length(indexapp)
            outputdata{line,i+datacol}=datalinenum(indexapp(i));
        end
    else
        break
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

msgbox('Species is successfully sorted and imported!');
fprintf('\nspecies_capture can screen out the interested species data, species_classfy can obtain the interested species with specifically structure characterization\n')
fprintf('\nData is saved in outputdata\n')
fprintf('\nspecies_analysis is successfully finished\n\n')
Elapsedtime = toc;
fprintf('\nTotal task time: %.2f s\n',Elapsedtime)

clear datacell datacellchar datacellnum datacol datadel datafirstrow dataline datalinenum datanow filename ans
clear dataoutcol dataoutcolchar dataoutputcol dataoutputrow datarow datarep datarow datasec datasplit statans
clear i j k line rawdata charcor datadelimiter dataoutrow matches outputans dataname indexapp indexcol indexovlp
clear Elapsedtime

