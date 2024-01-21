%scrit file name species_bonds_lammpstrj_framecheck
%purpose:
%This program is used to check and modify the frame number in the species
%file (outputnew) compared with the bonds and lammpstrj files
%version 1;2023.9.4

if ~exist('fram_num_check','var') && strcmpi(rerun_ans,'y')
    fprintf('Check the consistency of frame No. between species.* (outputnew), bonds.* and lammpstrj.* files.\nSometimes energy minimization causes the inconsistent record problem\n')
    fprintf('species_bonds_lammpstrj_framecheck is running, please wait...\n')
    fprintf('Firstly, check frame No. between bonds.* and species.* files, please wait...\n')
    fram_num_check=[];
    check_control=check_control_origin;
    control=1;
    rawdata=fopen(datanamebond,'r');
    readline=0;
    gap=8+atomnum;
    dataline=fgetl(rawdata);
    readline=readline+1;
    datacell=textscan(dataline,'%s','delimiter','\n');
    datacellchar=char(datacell{1});
    datadel=strrep(datacellchar,'#','');
    datarep=strtrim(datadel);
    datasplit=strsplit(datarep);
    fram_num_check(size(fram_num_check,1)+1,1)=str2num(datasplit{1,2});
    check_control=check_control-2;
    while check_control
        if check_control==0
            control=0;
        else
            while control
                i=1;
                unfound=1;
                while unfound
                    dataline=fgetl(rawdata);
                    readline=readline+1;
                    i=i+1;
                    if i==gap+1
                        unfound=0;
                        break;
                    end
                end
                if mod(readline-1,gap)==0
                    if feof(rawdata)
                        fprintf('\nReach the last line of bond.* file, please check the atom number or timestep. Sometimes timestep is not recorded as it is, eg. 29,429... other than 0,400...\n')
                        error('Check the aforementioned information!')
                    end
                    datacell=textscan(dataline,'%s','delimiter','\n');
                    datacellchar=char(datacell{1});
                    datadel=strrep(datacellchar,'#','');
                    datarep=strtrim(datadel);
                    datasplit=strsplit(datarep);
                    fram_num_check(size(fram_num_check,1)+1,1)=str2num(datasplit{1,2});
                    if check_control==0
                        control=0;
                    else
                        check_control=check_control-1;
                    end
                else
                    disp('Not a timestep line，please check the file or code!')
                    return;
                end
            end
        end
    end
    fclose(rawdata);
    fprintf('The appointed frame No. in bonds.* file (%d in total) has been read\n',check_control_origin)
    
    clear atomnumcopy ans bondnumdata control datacell datacellchar datadel dataline datarep datasplit found gap i j k kk line
    clear outputans unfound dataoutrow dataoutcol dataoutputrow dataoutcolchar dataoutputcol filename check_control
    
    
    
    
    fprintf('Then, check frame No. between lammpstrj.* and species.* files, please wait...\n')
    fram_num_check2=[];
    check_control=check_control_origin;
    control=1;
    rawdatatrj=fopen(datanametrj,'r');
    readline=0;
    gap=9+atomnum;
    
    dataline=fgetl(rawdatatrj);
    readline=readline+1;
    dataline=fgetl(rawdatatrj);
    readline=readline+1;
    datacell=textscan(dataline,'%s','delimiter','\n');
    datacellchar=char(datacell{1});
    datarep=strtrim(datacellchar);
    fram_num_check2(size(fram_num_check2,1)+1,1)=str2num(datarep);
    check_control=check_control-1;
    while check_control
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
            if feof(rawdatatrj)
                fprintf('Reach the last line of lammpstrj.* file, please check the atom number or timestep. Sometimes timestep is not recorded as it is, eg. 29,429... other than 0,400...\n')
                error('Check the aforementioned information!')
            end
            datacell=textscan(dataline,'%s','delimiter','\n');
            datacellchar=char(datacell{1});
            datarep=strtrim(datacellchar);
            fram_num_check2(size(fram_num_check2,1)+1,1)=str2num(datarep);
            if check_control==0
                control=0;
            else
                check_control=check_control-1;
            end
        else
            disp('Not a timestep line，please check the file or code!')
            return;
        end
    end
    fclose(rawdatatrj);
    fprintf('The appointed frame No. in lammpstrj.* file (%d in total) has been read\n',check_control_origin)
    
    fram_num_check(:,2)=fram_num_check2;

    fram_num_check(:,3)=cell2mat(outputdatanew(2:1+check_control_origin,1));
    fprintf('%d frame No. in the bonds.*, lammpstrj.* and species.* files are as follows: \n',check_control_origin)
    disp(fram_num_check)
    checkframe1=0;
    checkframe2=0;
    for checkframe=1:check_control_origin
        if ~ismember(fram_num_check(checkframe,3),fram_num_check(:,1))
            checkframe1=checkframe1+1;
        end
        if ~ismember(fram_num_check(checkframe,3),fram_num_check(:,2))
            checkframe2=checkframe2+1;
        end
    end
    
    
    if checkframe1>=3 || checkframe2>=3 || checkframe1>=check_control_origin/2 || checkframe2>=check_control_origin/2
        fprintf('As above, frame No. in the bonds.*, lammpstrj.* and species.* files are inconsistent, \nplease treat the frame No. of species.* file (outputnew) by the following method:\n')
        msgbox('Please treat the frame No. of species.* file (outputnew) by the following method')
        framecheck=input('Please treat the frame No. of species.* file (outputnew) by the following method: \n1.Substract an integer 2.Add an integer 3.No action\n');
        num_modify=input('Please input the unsigned number, eg.8, 29, or 377, etc.: \n');
        if framecheck==1
            for checknum=2:size(outputdatanew,1)
                outputdatanew{checknum,1}=outputdatanew{checknum,1}-num_modify;
            end
            tartrajectory=tartrajectory-num_modify;
            promptans3=promptans3-num_modify;
            fprintf('The initial frame No.is changed to %d\n',tartrajectory)
        elseif framecheck==2
            for checknum=2:size(outputdatanew,1)
                outputdatanew{checknum,1}=outputdatanew{checknum,1}+num_modify;
            end
            tartrajectory=tartrajectory+num_modify;
            promptans3=promptans3+num_modify;
            fprintf('The initial target frame No.is changed to %d\n%d',tartrajectory)
        elseif framecheck==3
            fprintf('No action is carried out\n')
        else
            error('Illegal method selection, please check it!')
        end
        fprintf('The frame No. problem in species.* file (outputnew) has been disposed\n')
    end
    
    clear atomnumcopy control datacell datacellchar dataline datasplit found gap i line unfound ans checkframe1 checkframe2 framecheck
    clear outputans dataoutrow dataoutcol dataoutputrow dataoutcolchar dataoutputcol check_control checkframe num_modify fram_num_check2
else
    fprintf('The frame No. of species.* file (outputnew) has been disposed and reused\n')
end

if rerun_ans2==1 
    species_frame_check2='n';
    if strcmpi(species_frame_checkans,'y')
        for iijj=3:size(outputdatanew,1)
            if mod(outputdatanew{iijj,1}-outputdatanew{iijj-1,1},trajper(1))~=0
                fprintf('From the %d th frame, viz from %d timestep the arithmetical progression is damaged\n',iijj,outputdatanew{iijj,1})
                fprintf('Frame No. of the %d th frame to %d th frame are:%d   %d   %d\n',iijj-1,iijj+1,outputdatanew{iijj-1,1},outputdatanew{iijj,1},outputdatanew{iijj+1,1})
                msgbox('The arithmetical progression is damaged in species.* file. How to cope with the problem?')
                species_frame_check3=input('If to treat the destructive arithmetical progression problem, y/n：\n','s');
                species_frame_check2='y';
                ijk=iijj;
                break
            else
                species_frame_check3='n';
            end
        end
        if strcmpi(species_frame_check2,'y')
            for jjkk=iijj+1:size(outputdatanew,1)
                if mod(outputdatanew{jjkk,1}-outputdatanew{jjkk-1,1},trajper(1))~=0
                    fprintf('For the destructive arithmetical progression, from %d th frame on, viz %d timestep, the arithmetical progression is damage again\n',jjkk,outputdatanew{jjkk,1})
                    msgbox('The arithmetical progression is damage again!')
                    break
                end
            end
        end
    else
        species_frame_check3='n';
        fprintf('Not to check the arithmetical progression problem for the modified frame No. of species.* file, good luck!\n')
    end
    if strcmpi(species_frame_check3,'y')
        for jjkk=ijk:size(outputdatanew,1)
            if tartrajectory==outputdatanew{jjkk,1}
                tartrajectory=outputdatanew{jjkk-1,1}+trajper(1);
                promptans3=outputdatanew{jjkk-1,1}+trajper(1);
                fprintf('The initial target frame No. is modified to %d once more\n',tartrajectory)
            end
            outputdatanew{jjkk,1}=outputdatanew{jjkk-1,1}+trajper(1);
        end
        fprintf('The arithmetical progression problem is repaired for species.* file (outputnew), but not always guaranteed for 100% sure\n')
    end
    
    for jjkk=3:size(outputdatanew,1)
        if tartrajectory>outputdatanew{jjkk-1,1} && tartrajectory<outputdatanew{jjkk,1}
            tartrajectory=outputdatanew{jjkk-1,1};
            promptans3=outputdatanew{jjkk-1,1};
        end
    end
end

clear iijj species_frame_check2 jjkk species_frame_check3 ijk 
