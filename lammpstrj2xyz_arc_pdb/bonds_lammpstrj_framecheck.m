%scrit file name bonds_lammpstrj_framecheck
%purpose:
%This program is used to check and modify the frame number in the species
%file (outputnew) compared with the bonds and lammpstrj files
%version 1;2023.9.4

if strcmpi(rerun_ans,'y')
    fprintf('Check the consistency of frame No. between bonds.* and lammpstrj.* files.\nSometimes energy minimization causes the inconsistent record problem\n')
    fprintf('species_bonds_lammpstrj_framecheck is running, please wait...\n')
    fprintf('Firstly, check frame No. in bonds.* file, please wait...\n')
    fram_num_check=[];
    check_control=check_control_origin;
    control=1;
    rawdata=fopen(bonddataname,'r');
    atomnum=atom_num_autoreadbond(bonddataname);
    dataline=fgetl(rawdata);
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
                    while isempty(dataline)
                        dataline=fgetl(rawdata);
                    end
                    i=i+1;
                    gap=8+atomnum;
                    if i==gap+1
                        unfound=0;
                        break;
                    end
                end
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
                a0 = ftell(rawdata);
                dataline=fgetl(rawdata);%前进两行获取原子数
                a1 = ftell(rawdata);
                dataline=fgetl(rawdata);
                a2 = ftell(rawdata);
                datacell=textscan(dataline,'%s','delimiter','\n');
                datacellchar=char(datacell{1});
                datadel=strrep(datacellchar,'#','');
                datarep=strtrim(datadel);
                datasplit=strsplit(datarep);
                atomnum = str2double(datasplit{length(datasplit)});
                fseek(rawdata,-(a2-a0),'cof');%回退两行
            end
        end
    end
    fclose(rawdata);
    fprintf('The appointed frame No. in bonds.* file (%d in total) has been read\n',check_control_origin)
    
    clear atomnumcopy ans bondnumdata control datacell datacellchar datadel dataline datarep datasplit found gap i j k kk line
    clear outputans unfound dataoutrow dataoutcol dataoutputrow dataoutcolchar dataoutputcol filename check_control

    fprintf('Then, check frame No. in lammpstrj.* file, please wait...\n')
    fram_num_check2=[];
    check_control=check_control_origin;
    control=1;
    rawdatatrj=fopen(trjdataname,'r');
    atomnum=atom_num_autoreadtrj(trjdataname);
    dataline=fgetl(rawdatatrj);
    dataline=fgetl(rawdatatrj);
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
            while isempty(rawdatatrj)
                dataline=fgetl(rawdata);
            end
            i=i+1;
            gap=9+atomnum;
            if i==gap+1
                unfound=0;
                break;
            end
        end
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
        a0 = ftell(rawdatatrj);
        dataline=fgetl(rawdatatrj);%前进两行获取原子数
        a1 = ftell(rawdatatrj);
        dataline=fgetl(rawdatatrj);
        a2 = ftell(rawdatatrj);
        datacell=textscan(dataline,'%s','delimiter','\n');
        datacellchar=char(datacell{1});
        datadel=strrep(datacellchar,'#','');
        datarep=strtrim(datadel);
        datasplit=strsplit(datarep);
        atomnum = str2double(datasplit{length(datasplit)});
        fseek(rawdatatrj,-(a2-a0),'cof');%回退两行
    end
    fclose(rawdatatrj);%
    fprintf('The appointed frame No. in lammpstrj.* file (%d in total) has been read\n',check_control_origin)
    
    fram_num_check(:,2)=fram_num_check2;
    fprintf('bonds_lammpstrj_framecheck is end\n')
    
    
    fprintf('%d frame No. in the bonds.* and lammpstrj.* are as follows: \n',check_control_origin)
    disp(fram_num_check)
    checkframe1=0;
    for checkframe=1:check_control_origin
        if ~ismember(fram_num_check(checkframe,2),fram_num_check(:,1))
            checkframe1=checkframe1+1;
        end
    end
    
    %
    if checkframe1>=3 || checkframe1>=check_control_origin/2
        fprintf('As above, frame No. in the bonds.* and lammpstrj.* files are inconsistent, \nplease treat the frame No. of lammpstrj.* file by the following method:\n')
        msgbox('Treat the frame No. of lammpstrj.* file by the following method')
        framecheck=input('Please treat the frame No. of lammpstrj.* file (outputnew) by the following method: \n1.Substract an integer 2.Add an integer 3.No action\n');
        num_modify=input('Please input the unsigned number, eg.8, 29, or 377, etc.: \n');
    else
        num_modify=0;
        fprintf('No action is carried out for frame No. in lammpstrj。* file, continue...\n')
    end
    
    clear atomnumcopy control datacell datacellchar dataline datasplit found gap i line unfound ans checkframe1 checkframe2 framecheck
    clear outputans dataoutrow dataoutcol dataoutputrow dataoutcolchar dataoutputcol check_control checkframe fram_num_check2
else
	num_modify = 0;
    fprintf('The frame No. problem in lammpstrj.* file has been disposed and reused\n')
end
 
clear  a0 a1 a2