%scrit file name Bond_Group_Analyzer
%purpose:
%This program is used to statistically analyze the evolved number of bond and group by bonds file
%version 1;2024.04.05

disp('##################################################################################################################################')
disp('Welcome!--by Qiang Liu @Institute of Nuclear Physics and Chemistry, China Academy of Engineering Physics; Email: liubinqiang@163.com');
disp('Repository adress of the Source code on github: https://github.com/dadaoqiuzhi/RMD_Digging');
disp('References: 1.Fuel 287 (2021) 119484. 2.ACS Appl. Mat. Interfaces 13(34) (2021) 41287-41302. 3.ACS Appl. Mat. Interfaces 2022, 14.(4), 5959-5972.')
disp('4.ACS Materials Letters 2023, 2174-2188. More work is coming!')

fprintf('\nThis program calculates and analyzes the changes in the number of specific chemical bonds and functional groups based on the bond files')
fprintf('Please update the corresponding search_bond_inform.txt or search_group_inform.txt files according to your needs and examples')
fprintf('\nPlease update the corresponding search_bond_inform.txt file based on the atom type numerical number according to the requirements and examples, \nwith one line for each chemical bond type')
fprintf('\nPlease update the corresponding search_group_inform.txt file based on the requirements and examples')
dataname=input('\n请输入要处理的bonds.*文件名：\n','s');
atomnum=atom_num_autoreadbond(dataname);
fprintf('\n1.Analyze chemical bonds  2.Analyze functional groups  3.Both 1 and 2')
Ana_choi = input('\nPlease select the analysis option1-3：\n');

Species_restrain = input('\nAre there restrictions on analyzing the structure of molecules, such as molecular weight, composition, etc? y/n(not applicable now):\n','s');

fprintf('\nPlease select the method number for exporting frames: \n1. Export all frames')
fprintf('\n2. Frame number determined by arithmetic difference method (including the first and last frames)')
fprintf('\n3. Manually input frames with time steps\n')
bond_frame = input('');
if bond_frame == 1
    fprintf('\nAnalyze and export all frames')
elseif bond_frame == 2
    bond_frame_num = input('\nPlease enter common difference value, positive integer：');
elseif bond_frame == 3
    bond_frame_num = input('\nPlease enter the time step matrix, such as[400 10000 400000]：\n');
    bond_frame_num = sort(bond_frame_num);
end

tic 
fprintf('\nBond_Group_Analyzer program is running, please wait...')

if Ana_choi == 1
    bond_target = load('search_bond_inform.txt');
    if isempty(bond_target)
        error('The chemical bond information is not specified in the search_bond_inform.txt file, please check it!')
    end
    bond_output = [];
elseif Ana_choi == 2
    n_bond = input('\nIf the bonded atom is a single bond, should it be used as the central atom in topological structure analysis to further \nclarify the search and determine the structure (such as -H，-F/Cl/Br, etc）？y/n，usually n：\n','s');
    if strcmpi(n_bond,'y')
        n_bond = 0;
    elseif strcmpi(n_bond,'n')
        n_bond = 1;
    else
        error('Incorrect input, please check it!')
    end
    group_target = dlmread('search_group_inform.txt');
    if isempty(group_target)
        error('The group information is not specified in the search_group_inform.txt file, please check it!')
    end
    group_output = [];
elseif Ana_choi == 3
    n_bond = input('\nIf the bonded atom is a single bond, should it be used as the central atom in topological structure analysis to further \nclarify the search and determine the structure (such as -H, -F/Cl/Br, etc.)？y/n：\n','s');
    if strcmpi(n_bond,'y')
        n_bond = 1;
    elseif strcmpi(n_bond,'n')
        n_bond = 0;
    else
        error('Incorrect input, please check it!')
    end
    fprintf('\nPlease update the corresponding search_bond_inform.txt and search_group_inform.txt files based on the requirements and examples')
    bond_target = load('search_bond_inform.txt');
    group_target = dlmread('search_group_inform.txt');
    if isempty(bond_target)
        error('The chemical bond information is not specified in the search_bond_inform.txt file, please check it!')
    end
    if isempty(group_target)
        error('The group information is not specified in the search_group_inform.txt file, please check!')
    end
    bond_output = [];
    group_output = [];
else 
    error('Incorrect input, please check it!')
end


frame_note = 1;
rawdata=fopen(dataname,'r');
Timestep_OK = 0;
while ~feof(rawdata)
    bondoutdata = [];
    dataline=fgetl(rawdata);
    if ~isempty(dataline)
        datacell=textscan(dataline,'%s','delimiter','\n');
        datacellchar=char(datacell{1});
        datarep=strtrim(datacellchar);
        if ~isempty(datarep) && contains(datarep,'Timestep')
            datasplit=strsplit(datarep);
            Timestep = str2num(datasplit{3});
            if frame_note == 1
                Timestep_ini = Timestep;
                frame_note = frame_note + 1;
            end
            a0 = ftell(rawdata);
            dataline=fgetl(rawdata);
            a1 = ftell(rawdata);
            dataline=fgetl(rawdata);
            a2 = ftell(rawdata);
            datacell=textscan(dataline,'%s','delimiter','\n');
            datacellchar=char(datacell{1});
            datadel=strrep(datacellchar,'#','');
            datarep=strtrim(datadel);
            datasplit=strsplit(datarep);
            atomnum = str2double(datasplit{length(datasplit)});
            fseek(rawdata,-(a2-a0),'cof');
        end
        if bond_frame == 1
            Timestep_OK = 1;
        elseif bond_frame == 2
            if rem(Timestep-Timestep_ini,bond_frame_num) == 0
                Timestep_OK = 1;
            end
        elseif bond_frame == 3
            if ismember(Timestep,bond_frame_num)
                Timestep_OK = 1;
            else
                Timestep_OK = 0;
            end
        end
        
        if Timestep_OK == 1
            for i=1:6
                dataline=fgetl(rawdata);
            end
            for i=1:atomnum 
                dataline=fgetl(rawdata);
                datacell=textscan(dataline,'%s','delimiter','\n');
                datacellchar=char(datacell{1});
                datarep=strtrim(datacellchar);
                datasplit=strsplit(datarep);
                if isscalar(datasplit) 
                    if strcmpi(datasplit,'#')
                        break
                    end
                end
                bondnumdata = datasplit(1:str2num(datasplit{3})+3);
                bondnumdata = str2num(char(bondnumdata))';
                bondoutdata(size(bondoutdata,1)+1,1:size(bondnumdata,2)) = bondnumdata;
            end
            
            
            
            if Ana_choi == 1 
                if frame_note == 2 
                    id_type = bondoutdata(:,1:2);
                    frame_note = frame_note + 1;
                end
                bond_target=IdType_Bond(bond_target,bondoutdata,id_type);
                fprintf('\nThe chemical bond information analysis in frame %d has been completed',Timestep)
                bond_output(size(bond_output,1)+1,1) = Timestep;
                bond_output(size(bond_output,1),2:1+size(bond_target,1)) = bond_target(:,5)';
                
                
                if ~feof(rawdata)
                    dataline = fgetl(rawdata);
                else
                    fprintf('\nArriving at the end of the bonds.* file')
                    break
                end
                
            elseif Ana_choi == 2 
                if frame_note == 2 
                    id_type = bondoutdata(:,1:2);
                    frame_note = frame_note + 1;
                end
                
                group_num = find(ismember(group_target(:,1),9999));
                if max(group_num) ~= size(group_target,1)
                    error('\nThere is an input format issue with the search_group_inform.txt file. Please refer to the search_group_inform_example.txt file')
                end
                for i = 1:size(group_num,1) 
                    group_num(i,2) = group_num(i,1);
                    if i == 1
                        group_num(i,1) = 1;
                    else
                        group_num(i,1) = group_num(i-1,2)+1;
                    end
                end
                group_output(size(group_output,1)+1,1) = Timestep;
                for i = 1:size(group_num,1) 
                    group_target_copy = group_target(group_num(i,1):group_num(i,2),:);
                    group_output(size(group_output,1),1+i) = IdType_topology_Group(group_target_copy,bondoutdata,id_type,n_bond);
                    
                end
                
                
                
                
                
                
                
                
                
                
                if ~feof(rawdata)
                    dataline = fgetl(rawdata);
                else
                    fprintf('\nArriving at the end of the bonds.* file')
                    break
                end
                 
            elseif Ana_choi == 3 
               if frame_note == 2 
                    id_type = bondoutdata(:,1:2);
                    frame_note = frame_note + 1;
                end
                bond_target=IdType_Bond(bond_target,bondoutdata,id_type);
                fprintf('\nThe chemical bond information analysis in frame %d has been completed',Timestep)
                bond_output(size(bond_output,1)+1,1) = Timestep;
                bond_output(size(bond_output,1),2:1+size(bond_target,1)) = bond_target(:,5)';
                
                
                
                
                
                
                
                
                if ~feof(rawdata)
                    dataline = fgetl(rawdata);
                else
                    fprintf('\nArriving at the end of the bonds.* file')
                    break
                end
                
                
            end
            Timestep_OK = 0; 
         
        else
            check1 = 0;
            for i = 1:atomnum+7 
                if ~feof(rawdata)
                    dataline = fgetl(rawdata);
                else
                    check1 = check1 + 1;
                    fprintf('\nArriving at the end of the bonds.* file')
                    break
                end
            end
            if check1 == 1
                break
            end
        end
        if bond_frame == 3 && Timestep == max(bond_frame_num) 
            fprintf('\nThe analysis task for the last specified frame has been completed, about to exit...')
            break
        end
        
    else
        fprintf('\n')
        warning('Blank line data appears, which may cause read errors')
    end
end
fclose(rawdata);




if Ana_choi == 1
    msgbox('Chemical bond analysis task is completed!');
    fprintf('\nBond_Group_Analyzer program is sinished')
    fprintf('\nThe data is stored in bond_output, and the column order of the chemical bond data is consistent with the row order of search_bond_inform.txt')
elseif Ana_choi == 2
    msgbox('Group analysis task is completed!');
elseif Ana_choi == 3
    msgbox('Chemical bond and functional group analysis tasks are completed!');
    fprintf('\nBond_Group_Analyzer program is sinished')
    fprintf('\nThe data is stored in bond_output, and the column order of the chemical bond data is consistent with the row order of search_bond_inform.txt')
    fprintf('\nThe data is stored in group_output, and the column order of the chemical bond data is consistent with the row order of search_group_inform.txt')
end


fprintf('\n')
Elapsedtime = toc;
fprintf('\nDuration of this run：%.2f s\n',Elapsedtime)

clear Ana_choi atomnum bond_frame bond_frame_num bondnumdata bondoutdata datacell datacellchar dataline dataname datarep
clear datasplit frame_note group_num i id_type rawdata Timestep Timestep_ini Timestep_OK a0 a1 a2 bond_target Elapsedtime
clear Species_restrain datadel ans