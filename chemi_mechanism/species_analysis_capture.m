%scrit file name species_capture
%purpose:
%This program is used to analysis species file
%version 1;2018.6.22
disp('##################################################################################################################################')
disp('Welcome!--by Qiang Liu @Institute of Nuclear Physics and Chemistry, China Academy of Engineering Physics; Email: liubinqiang@163.com');
disp('Repository adress of the Source code on github: https://github.com/dadaoqiuzhi/RMD_Digging');
disp('References: 1.Fuel 287 (2021) 119484. 2.ACS Appl. Mat. Interfaces 13(34) (2021) 41287-41302. 3.ACS Appl. Mat. Interfaces 2022, 14.(4), 5959-5972.')
disp('4.ACS Materials Letters 2023, 2174-2188. More work is coming!')
disp('##################################################################################################################################')
fprintf('\nspecies_analysis_capture is running, please wait...\n')

if exist('outputdatanew','var')
    if strcmpi(sperunans2,'n')
        outputdatanew_frame={};
        outputdatanew_frame=outputdatanew(:,1:3);
    end
end
species=strtrim(species);
species=strsplit(species);
outputdata_copy={};
rawdata=fopen(datanamespe,'r');
dataline=fgetl(rawdata);
datacell=textscan(dataline,'%s','delimiter','\n');
datacellchar=char(datacell{1});
datadel=strrep(datacellchar,'#','');
datarep=strtrim(datadel);
datasplit=strsplit(datarep);
datacellnum=length(datasplit);
species_copy(1,1:3)={'Timestep','No_Moles','No_Specs'};
for i=1:size(species,2)
    species_copy{1,3+i}=species{i};
end
outputdatanew={};
outputdatanew=species_copy(1,:);
outputdata_copy(1,:)=species_copy(1,:);
outputdata_copy(2,:)=num2cell(zeros([1 length(species_copy)]));
for j=1:datacellnum
    for k=1:size(outputdata_copy,2)
        if strcmpi(datasplit{j},outputdata_copy{1,k})
            outputdata_copy{2,k}=j;
        end
    end
end
datalinenum=fgetl(rawdata);
datalinenum=strtrim(datalinenum);
datalinenum=strread(datalinenum);
for i=1:size(outputdata_copy,2)
    if outputdata_copy{2,i}~=0
        outputdata_copy{2,i}=datalinenum(outputdata_copy{2,i});
    end
end
outputdatanew(2,:)=outputdata_copy(2,:);


while ~feof(rawdata)
    dataline=fgetl(rawdata);
    if ~isempty(dataline)
        datacell=textscan(dataline,'%s','delimiter','\n');
        datacellchar=char(datacell{1});
        datadel=strrep(datacellchar,'#','');
        datarep=strtrim(datadel);
        datasplit=strsplit(datarep);
        datacellnum=length(datasplit);
        
        outputdata_copy(2,:)=num2cell(zeros([1 length(species_copy)]));
        for j=1:datacellnum
            for k=1:size(outputdata_copy,2)
                if strcmpi(datasplit{j},outputdata_copy{1,k})
                    outputdata_copy{2,k}=j;
                end
            end
        end
        datalinenum=fgetl(rawdata);
        datalinenum=strtrim(datalinenum);
        datalinenum=strread(datalinenum);
        for i=1:size(outputdata_copy,2)
            if outputdata_copy{2,i}~=0
                outputdata_copy{2,i}=datalinenum(outputdata_copy{2,i});
            end
        end
        outputdatanew(size(outputdatanew,1)+1,:)=outputdata_copy(2,:);
    end
end
fclose(rawdata);

for i=1:size(outputdatanew,2) 
    if sum(cell2mat(outputdatanew(2:end,i)))==0
        fprintf('Species (%s) is not matched£¬please check it\n',outputdatanew{1,i});
    end
end

if exist('outputdatanew','var')
    if exist('outputdatanew_frame','var')
        outputdatanew(:,1:3)=outputdatanew_frame;%
    end
end
