%scrit file name species_capture
%purpose:
%This program is used to analysis species file
%version 1;2018.6.22
disp('Welcome!--by Qiang Liu @Institute of Nuclear Physics and Chemistry, China Academy of Engineering Physics; Email: liubinqiang@163.com');
disp('Repository adress of the Source code on github: https://github.com/dadaoqiuzhi/RMD_Digging');

disp('When species_analysis is executed, this procedure can obtain the interested species')
species=input('Please input the molecular formula, should be lin line with the species file, especially the element sequence, \nmultiple molecular formula can be seperated by white space: \n','s');
species=upper(species);
fprintf('\n\nspecies_capture is running, please wait...\n\n')
species=strtrim(species);
species=strsplit(species);
outputdatast=outputdata(1,:);
datamatch=[];
for i=1:length(species)
    for j=1:length(outputdatast)
        if strcmpi(species{i},(outputdatast{j}))
            datamatch(1,i)=i;
            datamatch(2,i)=j;
        end
    end
end

[checkrow,checkcol]=size(datamatch);
if checkcol~=length(species)
    if isempty(datamatch)
        fprintf('not match the interested species, please check it!');
        return;
    else
        fprintf('some species are not matched, please check it!');
        return;
    end
end
if sum(ismember(datamatch,0))>=1
    fprintf('match error, please check it and try to clear up the work space!');
    return;
end

outputdatanew={};
for k=1:3
    outputdatanew(:,k)=outputdata(:,k);
end
for j=1:length(species)
    outputdatanew(:,j+3)=outputdata(:,datamatch(2,j));
end
disp('Results of species_capture are saved in outputdatanew')

expoans=input('Export to excel? y/n: \n','s');
expoans=lower(expoans);
if expoans=='y'
    [dataoutrow,dataoutcol]=size(outputdatanew);
    dataoutputrow=strcat('A','1');
    dataoutcolchar=char(65+dataoutcol-1);
    dataoutputcol=strcat(dataoutcolchar,num2str(dataoutrow));
    filename='output_mydata.xlsx';
    xlswrite(filename,outputdatanew,dataoutputrow:dataoutputcol)
    disp('Results are species_capture are exported into excel: output_mydata')
end
fprintf('\n\nspecies_capture is successfully finished\n\n')
disp('Results of species_capture are saved in outputdatanew')
clear datamatch dataoutcol dataoutrow dataoutcolchar dataoutputcol  dataoutputrow i j k outputdatast species
clear filename expoans checkrow checkcol