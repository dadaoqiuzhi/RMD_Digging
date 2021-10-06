%scrit file name species_capture
%purpose:
%This program is used to analysis species file
%version 1;2018.6.22
disp('##################################################################################################################################')
disp('Welcome!--by Qiang Liu @Institute of Nuclear Physics and Chemistry, China Academy of Engineering Physics; Email: liubinqiang@163.com');
disp('Repository adress of the Source code on github: https://github.com/dadaoqiuzhi/RMD_Digging');
disp('References: 1.Fuel 287 (2021) 119484. 2.ACS Appl. Mat. Interfaces 13(34) (2021) 41287-41302. More work is coming!')
disp('##################################################################################################################################')
species=upper(species);
fprintf('\nspecies_capture is running, please wait...')
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
    else
        fprintf('some species are not matched, please check it!');
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
disp('Results of species_capture are saved in outputdatanew.')

fprintf('\nspecies_capture is successfully finished\n\n')
disp('Results of species_capture are saved in outputdatanew')

clear datamatch dataoutcol dataoutrow dataoutcolchar dataoutputcol  dataoutputrow i j k outputdatast
clear filename expoans checkrow checkcol