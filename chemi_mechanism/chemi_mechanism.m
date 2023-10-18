%scrit file name chemi_mechanism2
%purpose:
%This program is used to analysis chemical reaction mechanism according to
%molecular formula with the datafile of species and bonds
%include modified species_analysis, species_capture,
%bonds_analysis_speedup, bondorder_deepmining, lammpstrj_analysis and
%xyz_car_pdb_filemaker 
%version 1;2018.10.24
disp('##################################################################################################################################')
disp('Welcome!--by Qiang Liu @Institute of Nuclear Physics and Chemistry, China Academy of Engineering Physics; Email: liubinqiang@163.com');
disp('Repository adress of the Source code on github: https://github.com/dadaoqiuzhi/RMD_Digging');
disp('References: 1.Fuel 287 (2021) 119484. 2.ACS Appl. Mat. Interfaces 13(34) (2021) 41287-41302. 3.ACS Appl. Mat. Interfaces 2022, 14.(4), 5959-5972.')
disp('4.ACS Materials Letters 2023, 2174-2188. More work is coming!')
disp('##################################################################################################################################')
fprintf('This program is used to\n1.Export products or reactants files\n2.Analyze and export reactants-products files\n3.Analyze special groups, bonds and species (in development)\n4.Analyze where the interested species go\n')
choi=input('Please select the option No.: \n');
if choi==1 || choi==2 || choi==4
    formatout=input('\nPlease select the file format No.: 1.xyz 2.car  3.pdb\n');
    datanamespe=input('File name of species file: \n','s');
    species=input('Please input the molecular formula, element sequence should be consistent with the in.* file, multi input should be seperated by whith space: \n','s');
    species=upper(species);
end

datanamebond=input('File name of bonds file: \n','s');
if choi==1 || choi==2 || choi==4
    datanametrj=input('\nFile name of lammpstrj file: \n','s');
end
trajper=input('\nPlease input the output frequency of BO information and trajectory file (Positive integer, see bonds or lammpstrj file):  \n');

fprintf('Automatically read atom number from the *.lammpstrj file, please wait...')
atomnum=atom_num_autoread(datanametrj); 

check_control_origin=input('\nPlease input the number for frame No. check of bonds.* and *.lammpstrj files compared with species file��\navoiding mismatch induced error, >=5 is suggested: \n');
if check_control_origin <4
    warning('Maybe impossible to detect the dismatch problem for bonds.* and *.lammpstrj files compared with species file, >=5 is suggested!')
end
if exist('fram_num_check','var')
    rerun_ans=input('The check file for bonds.* and *.lammpstrj files is existent,, rerun check program (extra time is needed)? y/n: \n','s');
    if strcmpi(rerun_ans,'n')
        fprintf('According to your choice, try to input the corrected trajectory frame No., including output frequency and equal difference correction!\n')
        msgbox('According to your choice, try to input the corrected trajectory frame No.!')
    elseif ~ismember(rerun_ans,{'y','Y','n','N'})
        error('Illegal choice for rerun check of bonds.* and *.lammpstrj files compared with species file, please check it!')
    end
else
    rerun_ans=input('If to rerun check program (extra time is needed) for trajectory frame No. in bonds.* and *.lammpstrj files \ncompared with species file? y/n: \n','s');
end
rerun_ans2=1;

species_frame_checkans=input('Whether to check the arithmetic sequence of the trajectory frame No. in species file (output frequency). \nSometimes systematic error causes the inconsistency of common difference. y/n: \n','s');

elementsequence=input('\nPlease input atom type like C,H.O,N, seperated by white space, corresponding to 1,2,3,4...n (see *.data or in.*, \nespecially for element mapping:\n','s');
elementsequence=strtrim(elementsequence);elementsequence=strsplit(elementsequence);elementsequence=upper(elementsequence);
elemax=0;
for i=1:length(elementsequence)
    if length(elementsequence{i})>elemax
        elemax=length(elementsequence{i});
    end
end
eleswapans=input('\nDoes there exist element mapping?y/n: \n','s');
if strcmpi(eleswapans,'y')
    fprintf('\nPlease input element mapping cell, eg.{''Si'',''C'';''S'',''O''},\nsingle quotes are used in practical,element in data file--expected mapping elements: \n');
    eleswap=input('');
    eleswap=upper(eleswap);
else
    eleswap={'Nan'};
end
if elemax==2
    if formatout==1 || formatout==2
        fprintf('\nDifferent number system is adopted according to the atom number (ASCII)\n');
        if 262143>=atomnum && atomnum>32767
            fprintf('64 base number system is recommended for atom id');base=64;
        elseif 32767>=atomnum && atomnum>4095
            fprintf('32 base number system is recommended for atom id');base=32;
        elseif 4095>=atomnum && atomnum>999
            fprintf('16 base number system is recommended for atom id');base=16;
        elseif 999>=atomnum
            fprintf('10 base number system (decimalism) is recommended for atom id');
        elseif atomnum<0 || atomnum>32767
            error('Atom number is less than 0 or larger than  262143. If larger, please check it and modify code accordingly!')
        end
    elseif formatout==3
        fprintf('\nDifferent number system is adopted according to the atom number (ASCII)\n');
        if 65535>=atomnum && atomnum>16383
            fprintf('256 base number system is recommended for atom id');base=256;
        elseif 16383>=atomnum && atomnum>4095
            fprintf('128 base number system is recommended for atom id');base=128;
        elseif 4095>=atomnum && atomnum>1023
            fprintf('64 base number system is recommended for atom id');base=64;
        elseif 1023>=atomnum && atomnum>255
            fprintf('32 base number system is recommended for atom id');base=32;
        elseif 255>=atomnum && atomnum>99
            fprintf('16 base number system is recommended for atom id');base=16;
        elseif 99>=atomnum
            fprintf('10 base number system (decimalism) is recommended for atom id');
        elseif atomnum<0 || atomnum>65535
            error('Atom number is less than 0 or larger than  65535. If larger, please check it and modify code accordingly!')
        end
    end
elseif elemax==1
    if formatout==1 || formatout==2
        fprintf('\nDifferent number system is adopted according to the atom number (ASCII)\n');
        if 16777215>=atomnum && atomnum>1048575
            fprintf('64 base number system is recommended for atom id');base=64;
        elseif 1048575>=atomnum && atomnum>65535
            fprintf('32 base number system is recommended for atom id');base=32;
        elseif 65535>=atomnum && atomnum>9999
            fprintf('16 base number system is recommended for atom id');base=16;
        elseif 9999>=atomnum
            fprintf('10 base number system (decimalism) is recommended for atom id');
        elseif atomnum<0 || atomnum>16777215
            error('Atom number is less than 0 or larger than 16777215. If larger, please check it and modify code accordingly!')
        end
    elseif formatout==3
        fprintf('\n����ԭ�Ӹ������ò�ͬ���ƣ�ASCII���룩��Ԫ�����ϳ�ԭ�����ƣ��ܳ��Ȳ�����4���ַ���Ԫ������������ַ�����ԭ��idӦ������3���ַ�\n');
        if 262143>=atomnum && atomnum>32767
            fprintf('64 base number system is recommended for atom id');base=64;
        elseif 32767>=atomnum && atomnum>4095
            fprintf('32 base number system is recommended for atom id');base=32;
        elseif 4095>=atomnum && atomnum>999
            fprintf('16 base number system is recommended for atom id');base=16;
        elseif 999>=atomnum
            fprintf('10 base number system (decimalism) is recommended for atom id');
        elseif atomnum<0 || atomnum>262143
            error('Atom number is less than 0 or larger than 262143. If larger, please check it and modify code accordingly!')
        end
    end
else
    error('Length of atom type is NOT 1 or 2, please check it!');
end
base=input('\n��������õĽ�����������Ϊ�������Ҳ�С���Ƽ��Ľ�������\n');
unwrapans=input('\n�Ƿ�ִ������unwrap�����Խ�����Ժ���֮���ghost����Կ��ӻ���Ӱ��,y/n?: \n','s');

if choi==1
    maxchoi=input('\n������֡��ֵ��������ֵ������ʾ����ͬ�������в��ֳ�������,����ֵΪ10��\n');
end
if choi==1 || choi==2 || choi==4
    if formatout==2
        PBCchoi=input('\n�����������Ա߽���������ѡ�ON/OFF��\n','s');PBCchoi=upper(PBCchoi);
        if strcmp(PBCchoi,'ON')%д��file header��Ϣ
            PBC='PBC=ON';
            disp('�����Ա߽�����ʾ����PBC   33.4531   33.4531   33.4531   90.0000   90.0000   90.0000 (P1)��');
            PBCalpha=input('�����������Ա߽�����alpha������MS��lattice parameters�õ�����λС����\n');
            disp('�����Ա߽�����ʾ����PBC   33.4531   33.4531   33.4531   90.0000   90.0000   90.0000 (P1)��');
            PBCbeta=input('�����������Ա߽�����beta,����MS��lattice parameters�õ�����λС����\n');
            disp('�����Ա߽�����ʾ����PBC   33.4531   33.4531   33.4531   90.0000   90.0000   90.0000 (P1)��');
            PBCgamma=input('�����������Ա߽�����gamma,����MS��lattice parameters�õ�����λС����\n');
            disp('�����Ա߽�����ʾ����PBC   33.4531   33.4531   33.4531   90.0000   90.0000   90.0000 (P1)��');
            spacegroupname=input('�������������������Ա߽�����group name,ʾ������P1������\n','s');spacegroupname=upper(spacegroupname);
        elseif strcmp(PBCchoi,'OFF')
            PBC='PBC=OFF';
        else
            disp('PBCchoi�����Ա߽���������Ƿ�������');
            return;
        end
    elseif formatout==3
        PBCchoi=input('\n�����������Ա߽���������ѡ�ON/OFF��\n','s');PBCchoi=upper(PBCchoi);
        if strcmp(PBCchoi,'ON')
            PBCalpha=input('\n�����������Ա߽�����alpha������MS��lattice parameters�õ�����λС����\n');
            PBCbeta=input('\n�����������Ա߽�����beta,����MS��lattice parameters�õ�����λС����\n');
            PBCgamma=input('\n�����������Ա߽�����gamma,����MS��lattice parameters�õ�����λС����\n');
            spacegroupname=input('\n�������������������Ա߽�����group name,ʾ������*.arc�ǡ���P1��������*.pdb�ǡ�P 1����\n','s');spacegroupname=upper(spacegroupname);
        elseif strcmp(PBCchoi,'OFF')
            warndlg('��������Լ����������ȡ�����漰�����Ժ��ӳߴ磬ֱ����ȡ��')
        end
    end
    fprintf('\nǿ�ҽ��鲻Ҫ�������꣨ʹ��dump_modify scale no)�����׵��¼����������ͼ��覴�')
    BOXsize=input('\n����*.lammpstrj��ԭ��������Ƿ������ŵģ�dump_modify scale yes(no)����(�����ţ���y/n:\n','s');BOXsize=lower(BOXsize);
    if ~ismember(BOXsize,{'y','n'})
        error('�Ƿ���BOXsize����');
    end
    
    if exist('outputdatanew','var')
        fprintf('\n���棺species_analysis_capture����δ���У���Ϊ�����ռ����outputdatanew����,��ȷ���ǵ�ǰ�����ļ��ĵ�������');
        msgbox('���棺species_analysis_capture����δ����,��Ϊ�����ռ����outputdatanew����,�Ƿ���������species_analysis_capture��������');
        sperunans=input('\nspecies_analysis_capture����δ����,��Ϊ�����ռ����outputdatanew����,�Ƿ���������species_analysis_capture��������y/n��\n','s');
        if strcmpi(sperunans,'y')
            if strcmpi(rerun_ans,'n')%�����¼��֡����������������
                sperunans2=input('�Ƿ������޸��Ѵ��ڵ�species֡���(outputnew��)��ͬһ�ļ���֡��Ų�һ�����޸�������о���������ʱ����Ҫ�޸ģ�y/n��\n','s');%ͬһ�ļ��������У�����species֡���ʱ����ÿ����������ʱ
            else
                sperunans2='y';
            end
            species_analysis_capture %����Ŀ��species�ļ�
            msgbox('species_analysis_capture���н���');
        elseif strcmpi(sperunans,'n')
            fprintf('\n����ʹ�ù����ռ����еĴ洢�˲�����Ϣ��outputdatanew����');
            species=strtrim(species);
            species=strsplit(species);
        end
    end
    if ~exist('outputdatanew','var')
        species_analysis_capture %����Ŀ��species�ļ�
    end
    %����Ŀ�����ʽ�����仯��Ϣ������ȷ����������֡
    [row,~]=size(outputdatanew);
    residueseqname=species{1};%�л���ȷ��ΪĿ�����ʽ
    
    tic %��ʼ��ʱ
    
    prompt=0;rowcopy=row-1;irow=2;%������ʾ����,�ҵ��ܹ���Ҫ��������Ŀ
    while rowcopy
        if (irow==2 && outputdatanew{irow,4}~=0) ||(irow>2 && outputdatanew{irow,4}~=0 && outputdatanew{irow,4}~=outputdatanew{irow-1,4})
            prompt=prompt+1;
        end
        rowcopy=rowcopy-1;irow=irow+1;
    end
    
    if choi==1
        if prompt>maxchoi
            fprintf('\n��outputdatanew���ܹ���%d֡���㵼��Ҫ��,�ļ���Ŀ����֡��ֵ�����鵼������ͼ��',prompt);
            promptans=input('\n�Ƿ񵼳�����ͼ��,y/n��\n','s');
            if strcmpi(promptans,'y')
                promptans2=input('��ѡ��ȷ������֡�������ִ��룺\n1.����������֡\n2.����֡���ĵȲ������β��֡������֡������ͬ��\n3.���Ƶ���������֡�ĵȲ������β��֡����\n4.�ֶ�����ʱ�䲽��֡\n');
            end
        else
            promptans='n';
        end
    elseif choi==2 || choi==4
        promptans='y';promptans2=4;
    end
    
    
    if strcmpi(promptans,'y')%����ѡ��������������֡�������ٵ�֡
        frame={};%�洢������ѡ������֡��ʱ�䲽
        if promptans2==1
            maxdata=0;
            for irow=2:row
                if outputdatanew{irow,4}~=0
                    if outputdatanew{irow,4}>maxdata
                        maxdata=outputdatanew{irow,4};
                        [~,lenframe]=size(frame);frame{1,lenframe+1}=outputdatanew{irow,1};frame{2,lenframe+1}=outputdatanew{irow,4};
                    end
                end
            end
            [~,lenframe]=size(frame);
            for i=1:lenframe
                if frame{2,i}==0
                    frame(:,i)=[];
                end
            end
        elseif promptans2==2
            promptans3=input('\n��������Ҫ������֡��,����2��\n');
            counter=0;%�������Ŀprompt,���������֡����
            modcounter=mod(prompt-2,promptans3);
            moddivid=(prompt-2-modcounter)/(promptans3-2);
            if prompt<promptans3
                error('�ɵ�����֡������Ҫ���֡��');
            end
            for irow=2:row
                if (irow==2 && outputdatanew{irow,4}~=0) ||(irow>2 && outputdatanew{irow,4}~=0 && outputdatanew{irow,4}~=outputdatanew{irow-1,4})%ֻ�����в��������仯��֡
                    counter=counter+1;
                    if (counter==1 || counter==prompt) || (counter-1-modcounter>0 && mod(counter-1-modcounter,moddivid)==0)
                        [~,lenframe]=size(frame);frame{1,lenframe+1}=outputdatanew{irow,1};frame{2,lenframe+1}=outputdatanew{irow,4};
                    end
                end
            end
        elseif promptans2==3
            maxdata=0;counter=0;
            for irow=2:row
                if outputdatanew{irow,4}~=0
                    if outputdatanew{irow,4}>maxdata
                        maxdata=outputdatanew{irow,4};
                        counter=counter+1;%ͳ������Ҫ��ĵ���������֡��
                    end
                end
            end
            countercopy=counter;
            fprintf('\n����Ҫ��Ĺ���%d֡',countercopy);
            promptans3=input('\n��������Ҫ������֡����\n');
            modcounter=mod(counter-2,promptans3);
            moddivid=(counter-2-modcounter)/(promptans3-2);
            maxdata=0;counter=0;
            for irow=2:row
                if outputdatanew{irow,4}~=0
                    if outputdatanew{irow,4}>maxdata
                        maxdata=outputdatanew{irow,4};
                        counter=counter+1;
                        if (counter==1 || counter==countercopy) || (counter-1-modcounter>0 && mod(counter-1-modcounter,moddivid)==0)
                            [~,lenframe]=size(frame);frame{1,lenframe+1}=outputdatanew{irow,1};frame{2,lenframe+1}=outputdatanew{irow,4};
                        end
                    end
                end
            end
        elseif promptans2==4
            if choi==1
                promptans3=input('\n��������Ҫ������֡����ʱ�䲽���ÿպŸ�����\n','s');
            
            elseif choi==2 || choi==4
                fprintf('��ͨ��species_analysis��species_capture��ȷ����Ҫ�о���Ŀ�����ǰ��Ӧ���ض�ʱ�䲽');
                msgbox('species������ϣ�������������ڵ�֡����ʱ�䲽');
                promptans3=input('\n��������Ҫ�����Ĳ������ڵ�֡����ʱ�䲽���ÿպŸ�����\n','s');
            end
            
            promptans3=strtrim(promptans3);promptans3=strsplit(promptans3);promptans3=str2double(promptans3);
            frame=num2cell(promptans3);lenframe=length(frame);
            tartrajectory=promptans3;
            
            if exist('fram_num_check','var') && strcmpi(rerun_ans,'y') && rerun_ans2==1
                clear fram_num_check
                species_bonds_lammpstrj_framecheck%���species�ļ�֡���(outputnew��)��bonds��files�ļ��Ƿ�һ�£���һ��ʱ�޸�
                rerun_ans2=rerun_ans2+1;%֮��ѭ���������ټ��
            elseif exist('fram_num_check','var') && strcmpi(rerun_ans,'n')
                fprintf('species֡��outputnew�У�����Ѵ�����������\n')
            elseif ~exist('fram_num_check','var')
                species_bonds_lammpstrj_framecheck
                rerun_ans2=rerun_ans2+1;%֮��ѭ���������ټ�� 
            end
            if rerun_ans2==1
                species_bonds_lammpstrj_framecheck
                rerun_ans2=rerun_ans2+1;
            end
            frame=num2cell(promptans3);%���¸�ֵ�����޸����֡���
            
            for irow=1:lenframe
                if ~ismember(frame{irow},cell2mat(outputdatanew(2:end,1)))
                    error('�ֶ�����֡������һ֡�����ڣ����飡����');
                end
            end
            for irow=2:row
                for i=1:lenframe
                    if frame{1,i}==outputdatanew{irow,1}
                        frame{2,i}=outputdatanew{irow,4};
                    end
                end
            end
        else
            fprintf('�Ƿ��ĵ���֡�������ִ������룬���飡����');
        end
        
        
    elseif strcmpi(promptans,'n')
        frame={};rowcopy=row-1;irow=2;
        while rowcopy
            if (irow==2 && outputdatanew{irow,4}~=0) ||(irow>2 && outputdatanew{irow,4}~=0 && outputdatanew{irow,4}~=outputdatanew{irow-1,4})
                [~,lenframe]=size(frame);frame{1,lenframe+1}=outputdatanew{irow,1};frame{2,lenframe+1}=outputdatanew{irow,4};
            end
            rowcopy=rowcopy-1;irow=irow+1;
        end
        [~,lenframe]=size(frame);
        for i=1:lenframe
            if frame{2,i}==0
                frame(:,i)=[];
            end
        end
    end
end


if choi==1 && strcmpi(promptans,'y')
    rowcopy=0;[~,lenframe]=size(frame);
    rawdatatrj=fopen(datanametrj,'r');%�������ж�ȡ����lammpstrj�ļ�
    for irow=1:lenframe%ѭ���������ֳ�����ͼ
        rowcopy=rowcopy+1;
        fprintf('\n\nchemi_mechanism���������У���ȴ�...');
        fprintf('\n�ܹ���%d֡����ǰ���ɵ�%d֡car�ļ�',lenframe,rowcopy);
        tartrajectory={};
        tartrajectory=frame{1,irow};%Ŀ�����������켣
            
        if exist('fram_num_check','var') && strcmpi(rerun_ans,'y') && rerun_ans2==1
            clear fram_num_check
            species_bonds_lammpstrj_framecheck %���species�ļ�֡���(outputnew��)��bonds��files�ļ��Ƿ�һ�£���һ��ʱ�޸�
            rerun_ans2=rerun_ans2+1;%֮��ѭ���������ټ��
        elseif exist('fram_num_check','var') && strcmpi(rerun_ans,'n')
            fprintf('species֡��outputnew�У�����Ѵ�����������\n')
        elseif ~exist('fram_num_check','var')
            species_bonds_lammpstrj_framecheck
            rerun_ans2=rerun_ans2+1;%֮��ѭ���������ټ��
        end
        if rerun_ans2==1
            species_bonds_lammpstrj_framecheck
            rerun_ans2=rerun_ans2+1;
        end
        
        bonds_analysis_speedup%��ȡĿ��켣������Ϣ
        bondorder_deepmining%������ʽ��������ʣ��������Ŀ�����ʽͼ��
        
        [~,col]=size(tarelenummatch);bondset=[];
        for j=1:col/2%��װtarelenummatch�з���ʽ����Ŀ�����ʽ�ȶ�
            molecomp='';
            for k=1:length(elementsequence)
                if tarelenummatch{k,2*j}==1
                    molecomp=strcat(molecomp,tarelenummatch{k,2*j-1});
                elseif tarelenummatch{k,2*j}>1
                    molecomp=strcat(molecomp,tarelenummatch{k,2*j-1},num2str(tarelenummatch{k,2*j}));
                end
            end
            if strcmp(molecomp,species{1})%�ҵ�Ŀ��������ڵĵ�i�������tarBOinform��'#'�ָ���֮��
                bondsetdata=length(bondset);bondset(bondsetdata+1)=j;
            end
        end
        bondsetdata=length(bondset);
        if bondsetdata~=frame{2,irow}%����species�ļ���¼�ķ�������bondsͳ�Ƶ���Щ���ݾ�Ȼ�Բ��ϣ�����
            if bondsetdata==0
                fprintf('\n\n���棡����ͨ��bonds�ļ�δƥ�䵽Ŀ��������species�ļ���¼�˸�����(����Ϊ%d)��ԭ�����\n��in�ļ���fix reax/c/bondsΪ�ض�Nevery�����������Ϣ����fix reax/c/species��������ΪNeveryȡ������Nrepeatƽ��\n��Nfreq�����ƽ�����������һԭ����ͬ��ԭ��ӳ��Ϊ��ͣ��ԭ��ģ���ԭ��type���',frame{2,irow});
                fprintf('\n�������������in�ļ���fix reax/c/species��NrepeatΪ1������Nevery(species)=Nfreq(species)=Nevery(bonds),��ȫƥ����������Ϣ\n�滻ת��Ϊ���߱��������͡�\n���¼������ԭ��type\n');
                msgbox('δƥ�䵽Ŀ�����,����Ԫ��˳���Ƿ���ȷ����ѡ����ʽ��');
				errorexe=input('�Ƿ���Դ���y/n:\n','s');
				msgbox('�����⵽������ڣ��봦��');
                if strcmpi(errorexe,'y')
                    continue;%���Ա�֡
                elseif strcmpi(errorexe,'n')
                    error('��μ��Ϸ�������������飡����');
                end
            else
                fprintf('\n��ע��:bonds�ļ���ƥ�䵽��Ŀ���������%d)��species�ļ�������(%d)��һ��!!!\n',bondsetdata,frame{2,irow});
            end
        else
            fprintf('������Ŀ�����������%d)��������species�ļ�������(%d)һ��',bondsetdata,frame{2,irow});
        end
        [tarraw,~]=size(tarBOinform);ii=1;jj=1;
        while tarraw%��tarBOinform�з�Ŀ�����ļ�����Ϣɾ��
            if ismember(ii,bondset)
                if strcmp(tarBOinform{jj,1},'#')
                    ii=ii+1;
                end
                jj=jj+1;
                tarraw=tarraw-1;
            else
                if strcmp(tarBOinform{jj,1},'#')
                    ii=ii+1;
                end
                tarBOinform(jj,:)=[];
                tarraw=tarraw-1;
            end
        end
        fprintf('\n�ɹ��޳�tarBOinform�з�Ŀ�����ļ�����Ϣ\n');
        lammpstrj_analysis%��ȡĿ��켣������Ϣ
        %����ù켣�����пɱ�ʾΪĿ�����ʽ�Ĳ���car�ļ�������tarBOinform�еļ�����Ϣ�����¼��bondset�У�������Ϣ��trjdata��
        %���������unwrap
        if strcmpi(unwrapans,'y')
            load('BondRadii.mat');
            trjdata=PBC_Unwrap(tarBOinform,trjdata,BOXsize,boxsize,element,BondRadii);
        end
        if str2num(datarep)<=tartrajectory{1}%dump.*�ļ�֡ȷ������Ҫ������֡��dump.*������0ʱ�䲽��ʼ��¼���´���,С����С֡�������ڣ����޷�д�ļ�
            if formatout==1
                datanamecar=strcat(species{1},'-',num2str(frame{1,irow}),'-',num2str(frame{2,irow}),'-base',num2str(base),'.xyz');%��Ŀ������-ʵ����-���ʸ���Ϊcar�ļ���
            elseif formatout==2
                datanamecar=strcat(species{1},'-',num2str(frame{1,irow}),'-',num2str(frame{2,irow}),'-base',num2str(base),'.car');%��Ŀ������-ʵ����-���ʸ���Ϊcar�ļ���
            elseif formatout==3
                datanamecar=strcat(species{1},'-',num2str(frame{1,irow}),'-',num2str(frame{2,irow}),'-base',num2str(base),'.pdb');%��Ŀ������-ʵ����-���ʸ���Ϊcar�ļ���
            end
            xyz_car_pdb_filemaker%����xyz_car�ļ�
            seekBOinform%д�����������Ϣ
        else
            fclose(rawdatatrj);%�ر�lammpstrj�ļ�
            rawdatatrj=fopen(datanametrj,'r');%�������ж�ȡ����lammpstrj�ļ�,���´�
        end
    end
    
    
    
    
    
elseif choi==1 && strcmpi(promptans,'n')%�����������ڲ��ȵ�֡���ܶ࣡����
    rowcopy=0;[~,lenframe]=size(frame);%��Լ��������ƴ��
    rawdatatrj=fopen(datanametrj,'r');%�������ж�ȡ����lammpstrj�ļ�
    for irow=2:lenframe%���β��ɸ�irow,ѭ������
        rowcopy=rowcopy+1;
        fprintf('\n\nchemi_mechanism���������У���ȴ�...');
        fprintf('\n�ܹ���%d֡����ǰ���ɵ�%d֡car�ļ�\n',prompt,rowcopy);
        tartrajectory=frame{1,irow};%Ŀ�����������켣
        
        if exist('fram_num_check','var') && strcmpi(rerun_ans,'y') && rerun_ans2==1
            clear fram_num_check
            species_bonds_lammpstrj_framecheck%���species�ļ�֡���(outputnew��)��bonds��files�ļ��Ƿ�һ�£���һ��ʱ�޸�
            rerun_ans2=rerun_ans2+1;%֮��ѭ���������ټ��
        elseif exist('fram_num_check','var') && strcmpi(rerun_ans,'n')
            fprintf('species֡��outputnew�У�����Ѵ�����������\n')
        elseif ~exist('fram_num_check','var')
            species_bonds_lammpstrj_framecheck
            rerun_ans2=rerun_ans2+1;%֮��ѭ���������ټ��
        end
        if rerun_ans2==1
            species_bonds_lammpstrj_framecheck
            rerun_ans2=rerun_ans2+1;
        end
        bonds_analysis_speedup%��ȡĿ��켣������Ϣ
        bondorder_deepmining%������ʽ��������ʣ��������Ŀ�����ʽͼ��
        
        [~,col]=size(tarelenummatch);bondset=[];
        for j=1:col/2%��װtarelenummatch�з���ʽ����Ŀ�����ʽ�ȶ�
            molecomp='';
            for k=1:length(elementsequence)
                if tarelenummatch{k,2*j}==1
                    molecomp=strcat(molecomp,tarelenummatch{k,2*j-1});
                elseif tarelenummatch{k,2*j}>1
                    molecomp=strcat(molecomp,tarelenummatch{k,2*j-1},num2str(tarelenummatch{k,2*j}));
                end
            end
            if strcmp(molecomp,species{1})%�ҵ�Ŀ��������ڵĵ�i�������tarBOinform��'#'�ָ���֮��
                bondsetdata=length(bondset);bondset(bondsetdata+1)=j;
            end
        end
        bondsetdata=length(bondset);
        if bondsetdata~=outputdatanew{irow,4}%����species�ļ���¼�ķ�������bondsͳ�Ƶ���Щ���ݾ�Ȼ�Բ��ϣ�����
            if bondsetdata==0
                fprintf('\n\n���棡����ͨ��bonds�ļ�δƥ�䵽Ŀ��������species�ļ���¼�˸�����(����Ϊ%d)��ԭ�����\n��in�ļ���fix reax/c/bondsΪ�ض�Nevery�����������Ϣ����fix reax/c/species��������ΪNeveryȡ������Nrepeatƽ��\n��Nfreq�����ƽ�����������һԭ����ͬ��ԭ��ӳ��Ϊ��ͣ��ԭ��ģ���ԭ��type���',outputdatanew{irow,4});
                fprintf('\n�������������in�ļ���fix reax/c/species��NrepeatΪ1������Nevery(species)=Nfreq(species)=Nevery(bonds),��ȫƥ����������Ϣ\n�滻ת��Ϊ���߱���������\n');
                errorexe=input('�Ƿ���Դ���y/n:\n','s');
				msgbox('�����⵽������ڣ��봦������');
                if strcmpi(errorexe,'y')
                    continue;%���Ա�֡
                elseif strcmpi(errorexe,'n')
                    error('��μ��Ϸ�������������飡����');
                end
            else
                fprintf('\n��ע��:bonds�ļ���ƥ�䵽��Ŀ���������%d)��species�ļ�������(%d)��һ��!!!\n',bondsetdata,outputdatanew{irow,4});
            end
        else
            fprintf('������Ŀ�����������%d)��������species�ļ�������(%d)һ��',bondsetdata,outputdatanew{irow,4});
        end
        [tarraw,~]=size(tarBOinform);ii=1;jj=1;
        while tarraw%��tarBOinform�з�Ŀ�����ļ�����Ϣɾ��
            if ismember(ii,bondset)
                if strcmp(tarBOinform{jj,1},'#')
                    ii=ii+1;
                end
                jj=jj+1;
                tarraw=tarraw-1;
            else
                if strcmp(tarBOinform{jj,1},'#')
                    ii=ii+1;
                end
                tarBOinform(jj,:)=[];
                tarraw=tarraw-1;
            end
        end
        fprintf('\n�ɹ��޳�tarBOinform�з�Ŀ�����ļ�����Ϣ\n');
        lammpstrj_analysis%��ȡĿ��켣������Ϣ
        %����ù켣�����пɱ�ʾΪĿ�����ʽ�Ĳ���car�ļ�������tarBOinform�еļ�����Ϣ�����¼��bondset�У�������Ϣ��trjdata��
        %���������unwrap
        if strcmpi(unwrapans,'y')
            load('BondRadii.mat');
            trjdata=PBC_Unwrap(tarBOinform,trjdata,BOXsize,boxsize,element,BondRadii);
        end
        if str2num(datarep)<=tartrajectory{1}%dump.*�ļ�֡ȷ������Ҫ������֡��dump.*������0ʱ�䲽��ʼ��¼���´���,С����С֡�������ڣ����޷�д�ļ�
            if formatout==1
                datanamecar=strcat(species{1},'-',num2str(frame{1,irow}),'-',num2str(frame{2,irow}),'-base',num2str(base),'.xyz');%��Ŀ������-ʵ����-���ʸ���Ϊcar�ļ���
            elseif formatout==2
                datanamecar=strcat(species{1},'-',num2str(frame{1,irow}),'-',num2str(frame{2,irow}),'-base',num2str(base),'.car');%��Ŀ������-ʵ����-���ʸ���Ϊcar�ļ���
            elseif formatout==3
                datanamecar=strcat(species{1},'-',num2str(frame{1,irow}),'-',num2str(frame{2,irow}),'-base',num2str(base),'.pdb');%��Ŀ������-ʵ����-���ʸ���Ϊcar�ļ���
            end
            xyz_car_pdb_filemaker%����xyz_car�ļ�
            seekBOinform%д�����������Ϣ
        else
            fclose(rawdatatrj);%�ر�lammpstrj�ļ�
            rawdatatrj=fopen(datanametrj,'r');%�������ж�ȡ����lammpstrj�ļ�,���´�
        end
    end
    msgbox('����car�ļ��������');
    fclose(rawdatatrj);%�ر�lammpstrj�ļ�
    


    
    
    
elseif choi==2 || choi==4%ȷ���ò��������仯��������֡ʱ�䲽�������ݵ�һ��������Ķ�����Ӧ�ĳ�����Ż��߻�ѧ��
    [~,outcol]=size(outputdatanew);%����̫��ʱ�����ļ���Ŀ
    for i=4:outcol
        if strcmpi(species{1},(outputdatanew{1,i}))
            outcol=i;
            break
        end
    end
    if choi==2
        fprintf('\nĿ������%s�����ܹ���%d������ʼ%d����',species{1},outputdatanew{floor((promptans3(1)-outputdatanew{2,1})/trajper+2),outcol},outputdatanew{2,outcol});
        numstop=input('\n�Ƿ����Ƶ�����Ӧ��-����car�ļ�����Ŀ��������ƣ��ﵽ����Ŀ��ǿ����ֹ��Ѱ����y/n��\n','s');
        if strcmpi(numstop,'y')
            fprintf('\n�����벻����%d���ļ�������Ŀ',outputdatanew{floor((promptans3(1)-outputdatanew{2,1})/trajper+2),outcol})
            numstop=input(':\n');
        else
            numstop=outcol;
        end
    elseif choi==4
        fprintf('\n���ڷ�Ӧ��%d����һ��ȫ����ʧת������Сʣ��%d��������ʣ%d����������ݲ��ﵼ������ȷ��car�ļ�����Ŀ���ﵽ����Ŀ��ǿ����ֹ��Ѱ����',outputdatanew{floor((promptans3(1)-outputdatanew{2,1})/trajper+2),outcol},min(cell2mat(outputdatanew(2:end,outcol))),outputdatanew{end,outcol});
        numstop=input('\n�������ļ�������Ŀ:\n');
    end
    fprintf('\nΪ���������Ч�ʣ�����ֻ������Ŀ�仯��֡�������������ܲ��仯Ҳ����ǡ��������������������������Ǹ��ʺ�С');%���Ч�ʣ�ֻ��������Ŀ�仯��֡
    seekacc=input('\n�Ƿ�����ֻ������Ŀ�仯��֡������2���ӳ����ǵݼ���Ӧ���֡������4�ӳ����ǵ��������֡��y/n:\n','s');
    if ~strcmpi(seekacc,'y') && ~strcmpi(seekacc,'n')
        error('�Ƿ�������������Ŀ�仯��֡ѡ�����룬���飡����');
    end
    
    tartrajectory=frame{1,1};
    if exist('fram_num_check','var') && strcmpi(rerun_ans,'y') && rerun_ans2==1
        clear fram_num_check
        species_bonds_lammpstrj_framecheck%���species�ļ�֡���(outputnew��)��bonds��files�ļ��Ƿ�һ�£���һ��ʱ�޸�
        rerun_ans2=rerun_ans2+1;%֮��ѭ���������ټ��
    elseif exist('fram_num_check','var') && strcmpi(rerun_ans,'n')
        fprintf('species֡��outputnew�У�����Ѵ�����������\n')
        species_bonds_lammpstrj_framecheck%���species�ļ�֡���(outputnew��)��bonds��files�ļ��Ƿ�һ�£���һ��ʱ�޸�
        rerun_ans2=rerun_ans2+1;%֮��ѭ���������ټ��
    elseif ~exist('fram_num_check','var')
        species_bonds_lammpstrj_framecheck
        rerun_ans2=rerun_ans2+1;%֮��ѭ���������ټ��
    end
    if rerun_ans2==1
        species_bonds_lammpstrj_framecheck
        rerun_ans2=rerun_ans2+1;
    end
    
    fprintf('\n��ȡĿ��켣������Ϣ\n')
    bonds_analysis_speedup%��ȡĿ��켣������Ϣ
    bondorder_deepmining%������ʽ��������ʣ��������Ŀ�����ʽ���ں����ȶ����ͼ
    tarBOinformfullcopy=tarBOinform;%ȫ���ƣ������޳��ǲ����Ӱ�췴Ӧ�ﷴ�������Ѱ��Ӧԭ�����ڲ���
    [~,col]=size(tarelenummatch);bondset=[];
    for j=1:col/2%��װtarelenummatch�з���ʽ����Ŀ�����ʽ�ȶ�
        molecomp='';
        for k=1:length(elementsequence)
            if tarelenummatch{k,2*j}==1
                molecomp=strcat(molecomp,tarelenummatch{k,2*j-1});
            elseif tarelenummatch{k,2*j}>1
                molecomp=strcat(molecomp,tarelenummatch{k,2*j-1},num2str(tarelenummatch{k,2*j}));
            end
        end
        if strcmp(molecomp,species{1})%�ҵ�Ŀ��������ڵĵ�i�������tarBOinform��'#'�ָ���֮��
            bondsetdata=length(bondset);bondset(bondsetdata+1)=j;
        end
    end
    bondsetdata=length(bondset);
    if bondsetdata~=frame{2,1}%����species�ļ���¼�ķ�������bondsͳ�Ƶ���Щ���ݾ�Ȼ�Բ��ϣ�����
        if bondsetdata==0
            fprintf('\n\nͨ��bonds�ļ�δƥ�䵽Ŀ��������species�ļ���¼�˸�����(����Ϊ%d)��ԭ�����\n��in�ļ���fix reax/c/bondsΪ�ض�Nevery�����������Ϣ����fix reax/c/species��������ΪNeveryȡ������Nrepeatƽ��\n��Nfreq�����ƽ�����������һԭ����ͬ��ԭ��ӳ��Ϊ��ͬ��ԭ��ģ���ԭ��type˳���������',frame{2,1});
            fprintf('\n�������������in�ļ���fix reax/c/species��NrepeatΪ1������Nevery(species)=Nfreq(species)=Nevery(bonds),��ȫƥ����������Ϣ\n�滻ת��Ϊ���߱���������');
            msgbox('δƥ�䵽Ŀ�����,����Ԫ��˳���Ƿ���ȷ����ѡ����ʽ��');
			errorexe=input('\n�Ƿ���Դ���y/n:\n','s');
            if strcmpi(errorexe,'y')
%                 continue;%���Ա�֡
            elseif strcmpi(errorexe,'n')
                error('��μ��Ϸ�������������飡����');
            end
        else
            fprintf('\n��ע��:bonds�ļ���ƥ�䵽��Ŀ���������%d)��species�ļ�������(%d)��һ��!!!\n',bondsetdata,frame{2,1});
        end
    else
        fprintf('������Ŀ�����������%d)��������species�ļ�������(%d)һ��',bondsetdata,frame{2,1});
    end
    
    [tarraw,~]=size(tarBOinform);ii=1;jj=1;
    while tarraw%��tarBOinform�з�Ŀ�����ļ�����Ϣɾ��
        if ismember(ii,bondset)
            if strcmp(tarBOinform{jj,1},'#')
                ii=ii+1;
            end
            jj=jj+1;
            tarraw=tarraw-1;
        else
            if strcmp(tarBOinform{jj,1},'#')
                ii=ii+1;
            end
            tarBOinform(jj,:)=[];
            tarraw=tarraw-1;
        end
    end
    fprintf('\n�ɹ��޳�tarBOinform�з�Ŀ�����ļ�����Ϣ\n');
    rawdatatrj=fopen(datanametrj,'r');%�������ж�ȡ����lammpstrj�ļ�
    lammpstrj_analysis%��ȡ����켣������Ϣ
    %���������unwrap
    if strcmpi(unwrapans,'y')
        load('BondRadii.mat');
        trjdata=PBC_Unwrap(tarBOinform,trjdata,BOXsize,boxsize,element,BondRadii);
    end
    fclose(rawdatatrj);%�ر�lammpstrj�ļ�
    trjdatacopy=trjdata;%����Ŀ�����������Ϣ���������ǰ�˷�Ӧ�ﱻ����
    tarBOinformcopy=tarBOinform;%�����޳��˷�Ŀ�������Ϣ���Ŀ����������Ϣ���������ǰ�˷�Ӧ�ﱻ����
    tartrajectorycopy=tartrajectory{1};
    
    m=0;productnum=[];%��¼���е���Ӧ��Ĳ���Ͳ�������ԭ���ڲ���֡�и�����ļ�����Ϣ����Ӧ������trjdatacopy��
    [row5,~]=size(tarBOinformfullcopy);
    for i=1:row5%��¼����֡�и�������ʼ�к�����
        if strcmp(tarBOinformfullcopy{i,1},'#')
            m=m+1;
            if m==1
                productnum(m,1)=1;
                productnum(m,2)=i;
                productnum(m,3)=productnum(m,2)-productnum(m,1)+1;
            else
                productnum(m,1)=productnum(m-1,2)+1;
                productnum(m,2)=i;
                productnum(m,3)=productnum(m,2)-productnum(m,1)+1;
            end
        end
    end
    
    [tarraw,~]=size(tarBOinformcopy);loopnum=1;pairnum=0;%��ǰ���ݷ�Ӧ��֡
    loop=1;trajectorynote=[];trajectorynote(1,1)=tartrajectorycopy(1);trajectorynote(1,2)=tartrajectorycopy(1);%��ź���֡��ʱ�䲽��
    while tarraw%Ѱ��Ŀ��֡����ķ�Ӧ��֡����Ӧ���ڲ�ͬʱ�䲽�ķֿ�����xyz_car�ļ�
        if choi==2
            if strcmpi(seekacc,'n')
                tartrajectoryact=tartrajectorycopy(1)-loopnum*trajper;
            elseif strcmpi(seekacc,'y')
                control=1;
                while control
                    tartrajectoryact=tartrajectorycopy(1)-loop*trajper;
                    loop=loop+1;
                    if tartrajectoryact/trajper>=1 && outputdatanew{(tartrajectoryact-outputdatanew{2,1})/trajper+2,4}<outputdatanew{(trajectorynote(1,2)-outputdatanew{2,1})/trajper+2,4}
                        fprintf('\n�ҵ��µĲ���֡%d�����������Ŀ%d����һ֡%d��Ѱ����Ŀ%d��\n',outputdatanew{(tartrajectoryact-outputdatanew{2,1})/trajper+2,1},outputdatanew{(tartrajectoryact-outputdatanew{2,1})/trajper+2,4},outputdatanew{(trajectorynote(1,2)-outputdatanew{2,1})/trajper+2,1},outputdatanew{(trajectorynote(1,2)-outputdatanew{2,1})/trajper+2,4});
                        trajectorynote(1,1)=trajectorynote(1,2);trajectorynote(1,2)=tartrajectoryact;
                        control=0;
                        break;
                    elseif tartrajectoryact/trajper<1
                        error('��Ӧ��֡%d������һ֡����֡�����ڣ����飡����',tartrajectoryact);
                    end
                end
            end
            
        elseif choi==4
            if strcmpi(seekacc,'n')
                tartrajectoryact=tartrajectorycopy(1)+loopnum*trajper;
            elseif strcmpi(seekacc,'y')
                [row,~]=size(outputdatanew);
                control=1;
                while control
                    tartrajectoryact=tartrajectorycopy(1)+loop*trajper;
                    loop=loop+1;
                    if tartrajectoryact<=outputdatanew{row,1} && outputdatanew{(tartrajectoryact-outputdatanew{2,1})/trajper+2,4}<outputdatanew{(trajectorynote(1,2)-outputdatanew{2,1})/trajper+2,4}
                        fprintf('\n�ҵ��µĲ���֡%d�����������Ŀ%d����һ֡%d��Ѱ����Ŀ%d��\n',outputdatanew{(tartrajectoryact-outputdatanew{2,1})/trajper+2,1},outputdatanew{(tartrajectoryact-outputdatanew{2,1})/trajper+2,4},outputdatanew{(trajectorynote(1,2)-outputdatanew{2,1})/trajper+2,1},outputdatanew{(trajectorynote(1,2)-outputdatanew{2,1})/trajper+2,4});
                        trajectorynote(1,1)=trajectorynote(1,2);trajectorynote(1,2)=tartrajectoryact;
                        control=0;
                        break;
                    elseif tartrajectoryact>outputdatanew{row,1}
                        error('��Ӧ��֡%d�������֡����֡�����ڣ����飡����',tartrajectoryact);
                    end
                end
            end
            [row,~]=size(outputdatanew);
            if tartrajectoryact>outputdatanew{row,1}
                fprintf('\n�������ݵ����ʱ�䲽����Ѱ��ֹ')
                break;
            end
        end
        
        if tartrajectoryact/trajper>=1
            tartrajectory=tartrajectoryact;
            if exist('fram_num_check','var') && strcmpi(rerun_ans,'y') && rerun_ans2==1
                clear fram_num_check
                species_bonds_lammpstrj_framecheck%���species�ļ�֡���(outputnew��)��bonds��files�ļ��Ƿ�һ�£���һ��ʱ�޸�
                rerun_ans2=rerun_ans2+1;%֮��ѭ���������ټ��
            elseif exist('fram_num_check','var') && strcmpi(rerun_ans,'n')
                fprintf('species֡��outputnew�У�����Ѵ�����������\n')
            elseif ~exist('fram_num_check','var')
                species_bonds_lammpstrj_framecheck
                rerun_ans2=rerun_ans2+1;%֮��ѭ���������ټ��
            end
            if rerun_ans2==1
                species_bonds_lammpstrj_framecheck
                rerun_ans2=rerun_ans2+1;
            end
            
            bonds_analysis_speedup%��ȡĿ��켣������Ϣ
            bondorder_deepmining%������ʽ��������ʼ�����Ϣ����÷�Ӧ֡tarBOinform
            %��tarBOinformcopy�в����ڷ�Ӧ��֡��tarBOinform���ҵ���������
            
            frameproact=[];%���tarBOinformcopy�в���˳�����Ӧ��Ӧ������֡ʱ�䲽�ͷ�Ӧ��tarBOinform����������
            frameproact=react_blocklocate(tarBOinformcopy,tarBOinform,tartrajectory);%�õ���ʱ�䲽��Ӧ��tarBOinform����������
            
            if sum(frameproact(:,5))~=0%ȷ���з�Ӧ��
                if choi==2
                    fprintf('\n\n\n����������%d�鷴Ӧ��(%d)-����(%d)ƥ����,����Ϊ�����飬������%d��(������),�ܹ�%d��(����ͬһ��Ӧ��֡���ֶ����Ӧʱ),���Ƶ���%d��',loopnum,tartrajectory{1,1},tartrajectorycopy(1),pairnum+1,outputdatanew{floor((promptans3(1)-outputdatanew{2,1})/trajper+2),outcol},numstop);
                elseif choi==4
                    fprintf('\n\n\n����������%d�����(%d)-������(%d)ƥ����,����Ϊ�����飬������%d��(������),�����%d��(��ȫ��ת��Ϊ��Ĳ���ʱ),���Ƶ���%d��',loopnum,tartrajectory{1,1},tartrajectorycopy(1),pairnum+1,outputdatanew{floor((promptans3(1)-outputdatanew{2,1})/trajper+2),outcol},numstop);
                end
                loopnum=loopnum+1;
                %д��car�ļ�������������������Ϣ����д��ȡ������ʱ�䣻��д��Ӧ����д����
                tarBOinformcopy3={};%��¼��Ӧ�������Ϣ����Ӧ���껹δ��ȡ
                
                frameproact=react_blocklocate(tarBOinformcopy,tarBOinform,tartrajectory);%�õ���ʱ�䲽��Ӧ��tarBOinform����������/
                [lenframe,~]=size(frameproact);pronum=0;
                for i=1:lenframe
                    if frameproact(i,5)~=0
                        pronum=pronum+1;
                        tarBOinformcopy(frameproact(i,2):frameproact(i,3),:)=[];%����ҵ���Ӧ��Ĳ��������Ϣ
                        for j=i+1:lenframe%����ɾ�������к�frameproact�о�¹��tarBOinformcopy��ɾ�������µĸ�Ŀ�����������к�
                            frameproact(j,2)=frameproact(j,2)-frameproact(j,4);
                            frameproact(j,3)=frameproact(j,3)-frameproact(j,4);
                        end
                        [~,colframe]=size(frameproact);[numzero,~]=ismember(frameproact(i,:),0);%ȷ����Ӧ����Ŀ
                        numreactant=(colframe-5-sum(numzero))/2;
                        if numreactant>0%�ҵ���Ӧ�������Ϣ���洢��tarBOinformcopy3��
                            for j=1:numreactant
                                [row3,~]=size(tarBOinformcopy3);
                                tarBOinformcopy3(row3+1:row3+frameproact(i,2*j+5)-frameproact(i,2*j+4)+1,:)=tarBOinform(frameproact(i,2*j+4):frameproact(i,2*j+5),:);
                            end
                        else
                            if choi==2
                                error('�з�Ӧ�����У�����û���ҵ�������ţ����飡����')
                            elseif choi==4
                                error('�в������У�����û���ҵ�������ţ����飡����')
                            end
                        end
                    end
                end
                
                [row4,~]=size(tarBOinformcopy3);[row5,~]=size(tarBOinformfullcopy);%�ҵ���Ӧ������ԭ���ڲ�������֡���ڵ����ʣ�ע����������������������Ӧ����������ԭ������Ȼ��һ�����ذ���Ŀ�����
                blockNO=[];
                for k=1:row4
                    m=1;%����ÿһ���µķ�Ӧ�����飬��ʼ�µļ����ȶ�
                    for ii=1:row5
                        if strcmp(tarBOinformfullcopy{ii,1},'#')
                            m=m+1;
                        end
                        if ~strcmp(tarBOinformcopy3{k,1},'#') && ~strcmp(tarBOinformfullcopy{ii,1},'#') && tarBOinformcopy3{k,1}==tarBOinformfullcopy{ii,1}%ע������ֵ����ŵ�ASCII����ȣ�����'#'==35��Ϊ��ģ�����
                            if ~ismember(m,blockNO)%�����д����Ŀ������ظ�
                                lenbloc=length(blockNO);blockNO(lenbloc+1)=m;%�ҵ���Ӧ������ԭ���ڲ���֡�Ŀ�����
                            end
                        end
                    end
                end
                blockNO=sort(blockNO);
                
                
                tarBOinformcopy2={};
                for i=1:length(blockNO)%���¼�¼����֡������Ϣ������Ŀ�����ͷ�Ӧ������ԭ�����ڲ���
                    [row6,~]=size(tarBOinformcopy2);
                    tarBOinformcopy2(row6+1:row6+productnum(blockNO(i),3),:)=tarBOinformfullcopy(productnum(blockNO(i),1):productnum(blockNO(i),2),:);
                end
                
                %tartrajectory={tartrajectory(1)};%�õ�Ŀ��켣Ԫ����num��תcell��!!!
                rawdatatrj=fopen(datanametrj,'r');%�������ж�ȡ����lammpstrj�ļ�
                lammpstrj_analysis%��Ӧ��֡������Ϣ��tarBOinformcopy3�У�������Ϣ��trjdata�У����������Ϣ��tarBOinformcopy2�У�trjdata��trjdatacopy��
                 %���������unwrap
                 if strcmpi(unwrapans,'y')
                     load('BondRadii.mat');
                     trjdata=PBC_Unwrap(tarBOinform,trjdata,BOXsize,boxsize,element,BondRadii);
                 end
                fclose(rawdatatrj);%�ر�lammpstrj�ļ�
                if choi==2
                    if formatout==1
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','reactant','-base',num2str(base),'.xyz');%�ļ���Ϊ������-��Ӧ��ʱ�䲽-����ʱ�䲽-�������-��Ӧ��-�ļ���ʽ
                    elseif formatout==2
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','reactant','-base',num2str(base),'.car');%�ļ���Ϊ������-��Ӧ��ʱ�䲽-����ʱ�䲽-�������-��Ӧ��-�ļ���ʽ
                    elseif formatout==3
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','reactant','-base',num2str(base),'.pdb');%�ļ���Ϊ������-��Ӧ��ʱ�䲽-����ʱ�䲽-�������-��Ӧ��-�ļ���ʽ
                    end
                elseif choi==4
                    if formatout==1
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','product','-base',num2str(base),'.xyz');%�ļ���Ϊ������-��Ӧ��ʱ�䲽-����ʱ�䲽-�������-��Ӧ��-�ļ���ʽ
                    elseif formatout==2
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','product','-base',num2str(base),'.car');%�ļ���Ϊ������-��Ӧ��ʱ�䲽-����ʱ�䲽-�������-��Ӧ��-�ļ���ʽ
                    elseif formatout==3
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','product','-base',num2str(base),'.pdb');%�ļ���Ϊ������-��Ӧ��ʱ�䲽-����ʱ�䲽-�������-��Ӧ��-�ļ���ʽ
                    end
                end
                tarBOinform=tarBOinformcopy3;%�ֱ�д��Ӧ��Ͳ���дcar�ļ�����Ϊԭ����Ŀ�仯��arc�ļ������������ź����֡����дcar�ļ�
                xyz_car_pdb_filemaker%��ʼд�뷴Ӧ��
                
                if choi==2
                    if formatout==1
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','product','-base',num2str(base),'.xyz');%�ļ���Ϊ������-��Ӧ��ʱ�䲽-����ʱ�䲽-�������-����-�ļ���ʽ
                    elseif formatout==2
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','product','-base',num2str(base),'.car');%�ļ���Ϊ������-��Ӧ��ʱ�䲽-����ʱ�䲽-�������-����-�ļ���ʽ
                    elseif formatout==3
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','product','-base',num2str(base),'.pdb');%�ļ���Ϊ������-��Ӧ��ʱ�䲽-����ʱ�䲽-�������-����-�ļ���ʽ
                    end
                elseif choi==4
                    if formatout==1
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','reactant','-base',num2str(base),'.xyz');%�ļ���Ϊ������-��Ӧ��ʱ�䲽-����ʱ�䲽-�������-����-�ļ���ʽ
                    elseif formatout==2
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','reactant','-base',num2str(base),'.car');%�ļ���Ϊ������-��Ӧ��ʱ�䲽-����ʱ�䲽-�������-����-�ļ���ʽ
                    elseif formatout==3
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','reactant','-base',num2str(base),'.pdb');%�ļ���Ϊ������-��Ӧ��ʱ�䲽-����ʱ�䲽-�������-����-�ļ���ʽ
                    end
                end
                tarBOinform=tarBOinformcopy2;trjdata=trjdatacopy;
                xyz_car_pdb_filemaker%��ʼд�����
                
                seekBOinform%д��tarBOinformcopy3�з�Ӧ�choi=2������choi=4��������Ϣ,���ڲ���Ԫ�����꣬��MS��Ͻ�ʾ����
                
                pairnum=pairnum+1;
            else
                if choi==2
                    fprintf('\n\n\n��%d�鷴Ӧ��(%d)-����(%d)ƥ����δ�����κη�Ӧ�������%d��ƥ�䣬������һ֡���ݷ�Ӧ��\n',loopnum,tartrajectory{1,1},tartrajectorycopy(1),pairnum);
                elseif choi==4
                    fprintf('\n\n\n��%d�����(%d)-��Ӧ��(%d)ƥ����δ�����κη�Ӧ�������%d��ƥ�䣬������һ֡��Ѱ����\n',loopnum,tartrajectory{1,1},tartrajectorycopy(1),pairnum);
                end
                loopnum=loopnum+1;
            end  
        else
            error('��Ӧ��֡������һ֡����֡�����ڣ����飡����');
        end
        [tarraw,~]=size(tarBOinformcopy);%��������Ŀ�����֡���в���ķ�Ӧ��֡
        if pairnum==numstop
            tarraw=0;
		elseif pairnum>numstop
			fprintf('���棡������������ƥ��Դ�����ֹ���ƣ������ļ�������Ŀ����ֵ');
        end
    end
    if choi==2
        inspectBOinform%������choi=2��������Ϣ
        fprintf('\nĿ��֡%d���������ķ�Ӧ��-����������ϣ���%d�鷴Ӧ��-����car�ļ���',tartrajectorycopy(1),pairnum);
    elseif choi==4
        inspectBOinform%�����Ӧ�choi=4��������Ϣ
        fprintf('\nĿ��֡%d��Ӧ�������ķ�Ӧ��-����������ϣ���%d�鷴Ӧ��-����car�ļ���',tartrajectorycopy(1),pairnum);
    end
    
    msgbox('��Ӧ��-����������ϣ��ɹ�������car�ļ���');
   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
elseif choi==3
    fprintf('\n�����������Ļ��š����ɻ����߻�ѧ�����ϣ���鿴�ļ���������ָ����\n');
    fprintf('\n1.��ѧ����');
    
    
    
    
    
    
    
    
else
    disp('�Ƿ����룬δԤ����ӳ����ţ����飡����')
end
fprintf('\n\nchemi_mechanism���н���\n');
fprintf('Ŀ����š����ɻ����߻�ѧ����ʱ�䲽�仯����Ϣ������tarstepspecies��\n');
msgbox('chemi_mechanism���н���');

toc %������ʱ
fprintf('\n�������к�ʱ��%.2f s\n',toc)

clear blockNO bondset bondsetdata BOXsize choi colframe datanamebond datanametrj date elementsequence fidBO fileheader frame frameproact
clear ii irow jj lenbloc lenframe loopnum m molecomp numreactant numstop numzero outcol pairnum PBC PBCa PBCalpha PBCb PBCbeta PBCc
clear PBCchoi PBCgamma productnum prompt promptans promptans2 promptans3 pronum rawdata rawdatatrj reactname residueseqname
clear row3 row4 row5 row6 rowcopy rowfram spacegroupname species tarelenummatch tarraw tartrajectory tartrajectoryact tartrajectorycopy
clear tartrjdata title trajper loop datanamespe i j seekacc sperunans trajectorynote producttname col maxchoi maxdata outputdata_copy
clear atomid_conv base boxsize counter datanamecar eleswap eleswapans modcounter moddivid formatout datarep errorexe fram_num_check
clear rerun_ans sperunans2 outputdatanew_frame rerun_ans2 species_frame_checkans unwrapans species_frame_check2 readline MOLE
clear iijj elemax datalinenum datacellnum control checknum check_control_origin 
