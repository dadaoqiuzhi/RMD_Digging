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
fprintf('This program is used to\n1.Export structure files of products or reactants \n2.Analyze and export products->reactants files\n3.Analyze special groups, bonds and species (in development)\n4.Analyze where the interested species go (reactants->products)\n')
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
fprintf('\nPlease input the output frequency matrix of species, BO information and trajectory files (Positive integer, see species, bonds or lammpstrj file)');
trajper=input('\neg., [400 400 400]:  \n');

fprintf('Automatically read atom number from the *.lammpstrj file, please wait...')
atomnum=atom_num_autoread(datanametrj); 

check_control_origin=input('\nPlease input the number for frame No. check of bonds.* and *.lammpstrj files compared with species file，\navoiding mismatch induced error, >=5 is suggested: \n');
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
        if 238327>=atomnum && atomnum>32767
            fprintf('62 base number system is recommended for atom id');base=62;
        elseif 32767>=atomnum && atomnum>4095
            fprintf('32 base number system is recommended for atom id');base=32;
        elseif 4095>=atomnum && atomnum>999
            fprintf('16 base number system is recommended for atom id');base=16;
        elseif 999>=atomnum
            fprintf('10 base number system (decimalism) is recommended for atom id');
        elseif atomnum<=916132831 || atomnum>238327
            fprintf('\nAtom number is no more than 916132831 but larger than  238327. If larger, no element name is list in the number system, indicating 5 ASCII chars are used to encode\n')
			base=62;
		elseif atomnum > 916132831
            error('Atom number exceeds 916132831, which can not be encoded. Please check it or contact the developer!')
        end
    elseif formatout==3
        fprintf('\nDifferent number system is adopted according to the atom number (ASCII)\n');
        if 3843>=atomnum && atomnum>1023
            fprintf('62 base number system is recommended for atom id');base=62;
        elseif 1023>=atomnum && atomnum>255
            fprintf('32 base number system is recommended for atom id');base=32;
        elseif 255>=atomnum && atomnum>99
            fprintf('16 base number system is recommended for atom id');base=16;
        elseif 99>=atomnum
            fprintf('10 base number system (decimalism) is recommended for atom id');
        elseif atomnum<=14776335 || atomnum>3843
            fprintf('\nAtom number is no more than 14776335 but larger than 3843. If larger, no element name is list in the number system, indicating 4 ASCII chars are used to encode\n')
			base=62;
		elseif atomnum > 14776335
            error('Atom number exceeds 14776335, which can not be encoded. Please check it or contact the developer!')
        end
    end
elseif elemax==1
    if formatout==1 || formatout==2
        fprintf('\nDifferent number system is adopted according to the atom number (ASCII)\n');
        if 14776335>=atomnum && atomnum>1048575
            fprintf('62 base number system is recommended for atom id');base=62;
        elseif 1048575>=atomnum && atomnum>65535
            fprintf('32 base number system is recommended for atom id');base=32;
        elseif 65535>=atomnum && atomnum>9999
            fprintf('16 base number system is recommended for atom id');base=16;
        elseif 9999>=atomnum
            fprintf('10 base number system (decimalism) is recommended for atom id');
        elseif atomnum<=916132831 || atomnum>14776335
            fprintf('\nAtom number is no more than 916132831 but larger than 14776335. If larger, no element name is list in the number system, indicating 5 ASCII chars are used to encode\n')
			base=62;
		elseif atomnum > 916132831
            error('Atom number exceeds 916132831, which can not be encoded. Please check it or contact the developer!')
        end
    elseif formatout==3
        fprintf('\nDifferent number system is adopted according to the atom number (ASCII)\n');
        if 238327>=atomnum && atomnum>32767
            fprintf('62 base number system is recommended for atom id');base=62;
        elseif 32767>=atomnum && atomnum>4095
            fprintf('32 base number system is recommended for atom id');base=32;
        elseif 4095>=atomnum && atomnum>999
            fprintf('16 base number system is recommended for atom id');base=16;
        elseif 999>=atomnum
            fprintf('10 base number system (decimalism) is recommended for atom id');
        elseif atomnum<= 14776335 || atomnum>238327
            fprintf('Atom number is no more than 14776335 but larger than 238327. If larger, no element name is list in the number system, indicating 4 ASCII chars are used to encode')
			base=62;
		elseif atomnum > 14776335
            error('Atom number exceeds 14776335, which can not be encoded. Please check it or contact the developer!')
        end
    end
else
    error('Length of atom type is NOT 1 or 2, please check it!');
end
base=input('\nPlease select number system for atom id, positive integer and not less than the recommended value: \n');
unwrapans=input('\nIf to perform unwrap for atoms with coordinate affected by ghost position, causing discontinuity of structure,y/n? (Unusable Now!): \n','s');

if choi==1
    maxchoi=input('\nPlease input threshold trajectory frame value, when exceeding it different methods are suggested to limit the exportation of files\n');
end
if choi==1 || choi==2 || choi==4
    if formatout==2
        PBCchoi=input('\nPlease input periodic boundary condition, ON/OFF: \n','s');PBCchoi=upper(PBCchoi);
        if strcmp(PBCchoi,'ON')
            PBC='PBC=ON';
            disp('Periodic boundary condition, eg."PBC   33.4531   33.4531   33.4531   90.0000   90.0000   90.0000 (P1)"');
            PBCalpha=input('Periodic boundary condition, alpha, four decimal digits: ');
            disp('Periodic boundary condition, eg."PBC   33.4531   33.4531   33.4531   90.0000   90.0000   90.0000 (P1)"');
            PBCbeta=input('Periodic boundary condition, beta, four decimal digits:\n');
            disp('Periodic boundary condition, eg."PBC   33.4531   33.4531   33.4531   90.0000   90.0000   90.0000 (P1)"');
            PBCgamma=input('Periodic boundary condition, gamma, four decimal digits:\n');
            disp('Periodic boundary condition, eg."PBC   33.4531   33.4531   33.4531   90.0000   90.0000   90.0000 (P1)"');
            spacegroupname=input('Point group name, eg."(P1)": \n','s');spacegroupname=upper(spacegroupname);
        elseif strcmp(PBCchoi,'OFF')
            PBC='PBC=OFF';
        else
            disp('Illegal periodic boundary condition in PBCchoi, please check it!!!');
            return;
        end
    elseif formatout==3
        PBCchoi=input('\nPlease input periodic boundary condition, ON/OFF: \n','s');PBCchoi=upper(PBCchoi);
        if strcmp(PBCchoi,'ON')
            PBCalpha=input('\nPeriodic boundary condition, alpha, four decimal digits: \n');
            PBCbeta=input('\nPeriodic boundary condition, beta, four decimal digits:\n');
            PBCgamma=input('\nPeriodic boundary condition, gamma, four decimal digits:\n\n');
            spacegroupname=input('\nPoint group name, eg."(P1)": \n','s');spacegroupname=upper(spacegroupname);
        elseif strcmp(PBCchoi,'OFF')
            warndlg('\nNonperiodic constraint: PBC condition is not used and coordinate is directly abstracted!')
        end
    end
    fprintf('\nScaled coordinate is not recommended (use "dump_modify scale no" to avoid)')
    BOXsize=input('\nDoes the coordinate is scaled in the *.lammpstrj file, y/n: \n','s');BOXsize=lower(BOXsize);
    if ~ismember(BOXsize,{'y','n'})
        error('Illegal BOXsize parameters, please check it!!!');
    end
    
    if exist('outputdatanew','var')
        fprintf('\nWarning: species_analysis is not running due to the existence of outputdata in work space. \nPlease make sure this is expected species file');
        msgbox('Warning: species_analysis is not running due to the existence of outputdata in work space.\nRerun species_analysis to import species file?');
        sperunans=input('\nspecies_analysis is not running due to the existence of outputdata in work space.\nRerun species_analysis to import species file?y/n: \n','s');
        if strcmpi(sperunans,'y')
            if strcmpi(rerun_ans,'n')
                sperunans2=input('If to modify the frame No. from species.* file (outputnew), the same species.* file it is not recommended! y/n：\n','s');
            else
                sperunans2='y';
            end
            species_analysis_capture 
            msgbox('species_analysis is successfully finished.');
        elseif strcmpi(sperunans,'n')
            fprintf('\nSpecies data in outputdata is reused.');
            species=strtrim(species);
            species=strsplit(species);
        end
    end
    if ~exist('outputdatanew','var')
        species_analysis_capture 
    end
    
    [row,~]=size(outputdatanew);
    residueseqname=species{1};
    
    tic 
    
    prompt=0;rowcopy=row-1;irow=2;
    while rowcopy
        if (irow==2 && outputdatanew{irow,4}~=0) ||(irow>2 && outputdatanew{irow,4}~=0 && outputdatanew{irow,4}~=outputdatanew{irow-1,4})
            prompt=prompt+1;
        end
        rowcopy=rowcopy-1;irow=irow+1;
    end
    
    if choi==1
        if prompt>maxchoi
            fprintf('\n%d trajectory frame(s) can be exported in outputdatanew, which is larger the threshold trajectory frame value.\nPartial trajectory frame(s) are suggested to be exported',prompt);
            promptans=input('\nExport Partial trajectory frame(s),y/n: \n','s');
            if strcmpi(promptans,'y')
                promptans2=input('Please select the export method No.:\n1.Monotonically increasing frame(s)\n2.Frame(s) in arithmetic sequence(closed interval)\n3.Monotonically increasing frame(s) in arithmetic sequence(closed interval)\n4.Manually specify\n');
            end
        else
            promptans='n';
        end
    elseif choi==2 || choi==4
        promptans='y';promptans2=4;
    end
    
    
    if strcmpi(promptans,'y')
        frame={};
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
            promptans3=input('\nPlease input trajectory number to be exported, >2: \n');
            counter=0;
            modcounter=mod(prompt-2,promptans3);
            moddivid=(prompt-2-modcounter)/(promptans3-2);
            if prompt<promptans3
                error('Trajectory frame(s) can be exported is less than the expected input!!!');
            end
            for irow=2:row
                if (irow==2 && outputdatanew{irow,4}~=0) ||(irow>2 && outputdatanew{irow,4}~=0 && outputdatanew{irow,4}~=outputdatanew{irow-1,4})
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
                        counter=counter+1;
                    end
                end
            end
            countercopy=counter;
            fprintf('\n%d trajectory frame(s) in total meet requirements',countercopy);
            promptans3=input('\nPlease input trajectory number to be exported: \n');
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
                promptans3=input('\nPlease input the timestep expected to be exported, seperated by white space: \n','s');
            
            elseif choi==2 || choi==4
                fprintf('User can confirm the timestep by species_analysis, species_capture program.');
                msgbox('Species analasis is finished, please input the timestep of the interwsted species.');
                promptans3=input('\nPlease input the timestep expected to be exported, seperated by white space: \n','s');
            end
            
            promptans3=strtrim(promptans3);promptans3=strsplit(promptans3);promptans3=str2double(promptans3);
            frame=num2cell(promptans3);lenframe=length(frame);
            tartrajectory=promptans3;
            
            if exist('fram_num_check','var') && strcmpi(rerun_ans,'y') && rerun_ans2==1
                clear fram_num_check
                species_bonds_lammpstrj_framecheck
                rerun_ans2=rerun_ans2+1;
            elseif exist('fram_num_check','var') && strcmpi(rerun_ans,'n')
                fprintf('The frame No. from species.* file (outputnew) has been treated and reused\n')
            elseif ~exist('fram_num_check','var')
                species_bonds_lammpstrj_framecheck
                rerun_ans2=rerun_ans2+1;
            end
            if rerun_ans2==1
                species_bonds_lammpstrj_framecheck
                rerun_ans2=rerun_ans2+1;
            end
            frame=num2cell(promptans3);
            
            for irow=1:lenframe
                if ~ismember(frame{irow},cell2mat(outputdatanew(2:end,1)))
                    error('At least one of the trajectory frame(s) specified manually is not existent, please check it!!!');
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
            fprintf('Illegal method No. input, please check it!');
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
    rawdatatrj=fopen(datanametrj,'r');
    for irow=1:lenframe
        rowcopy=rowcopy+1;
        fprintf('\n\nchemi_mechanism is running, please wait...');
        fprintf('\n%d trajectory frame(s) in total. Now group %d is being processing',lenframe,rowcopy);
        tartrajectory={};
        tartrajectory=frame{1,irow};
            
        if exist('fram_num_check','var') && strcmpi(rerun_ans,'y') && rerun_ans2==1
            clear fram_num_check
            species_bonds_lammpstrj_framecheck 
            rerun_ans2=rerun_ans2+1;
        elseif exist('fram_num_check','var') && strcmpi(rerun_ans,'n')
            fprintf('The frame No. from species.* file (outputnew) has been treated and reused\n')
        elseif ~exist('fram_num_check','var')
            species_bonds_lammpstrj_framecheck
            rerun_ans2=rerun_ans2+1;
        end
        if rerun_ans2==1
            species_bonds_lammpstrj_framecheck
            rerun_ans2=rerun_ans2+1;
        end
        
        bonds_analysis_speedup
        bondorder_deepmining
        
        [~,col]=size(tarelenummatch);bondset=[];
        for j=1:col/2
            molecomp='';
            for k=1:length(elementsequence)
                if tarelenummatch{k,2*j}==1
                    molecomp=strcat(molecomp,tarelenummatch{k,2*j-1});
                elseif tarelenummatch{k,2*j}>1
                    molecomp=strcat(molecomp,tarelenummatch{k,2*j-1},num2str(tarelenummatch{k,2*j}));
                end
            end
            if strcmp(molecomp,species{1})
                bondsetdata=length(bondset);bondset(bondsetdata+1)=j;
            end
        end
        bondsetdata=length(bondset);
        if bondsetdata~=frame{2,irow}
            if bondsetdata==0
                fprintf('\n\nWarning!!!No species is found through bonds file, but species file record it(number:%d).Possible reason:\nNevery parameter of "fix reax/c/bonds" in the in.* file is not match with the Nrepeat"fix reax/c/species"\n. Another possible reason element mapping or type error',frame{2,irow});
                fprintf('\nSolution:Nrepeat of "fix reax/c/species" in the in.* file should be 1, and Nevery(species)=Nfreq(species)=Nevery(bonds)\nReplace the mapping element to the expected one before data processing\n');
                msgbox('Target species is not found!How to cope with this problem?');
				errorexe=input('Ignore this problem and continue?y/n: \n','s');
				msgbox('Problem are detected by the program check some, How to cope with this problem?');
                if strcmpi(errorexe,'y')
                    continue;
                elseif strcmpi(errorexe,'n')
                    error('Please see the aforesaid analysis and solution,please check it!!!');
                end
            else
                fprintf('\nCaution:target species number (%d) is not consistent with that of species recordation(%d)!!!\n',bondsetdata,frame{2,irow});
            end
        else
            fprintf('Target species number (%d) is consistent with that of species recordation(%d)',bondsetdata,frame{2,irow});
        end
        [tarraw,~]=size(tarBOinform);ii=1;jj=1;
        while tarraw
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
        fprintf('\nSuccessfully delete the unexpected BO information in tarBOinform\n');
        lammpstrj_analysis
        
        %unwrap
        if strcmpi(unwrapans,'y')
            load('BondRadii.mat');
            trjdata=PBC_Unwrap(tarBOinform,trjdata,BOXsize,boxsize,elementsequence,BondRadii);
        end
        if str2num(datarep)<=tartrajectory{1}
            if formatout==1
                datanamecar=strcat(species{1},'-',num2str(frame{1,irow}),'-',num2str(frame{2,irow}),'-base',num2str(base),'.xyz');%以目标物质-实践步-物质个数为car文件名
            elseif formatout==2
                datanamecar=strcat(species{1},'-',num2str(frame{1,irow}),'-',num2str(frame{2,irow}),'-base',num2str(base),'.car');%以目标物质-实践步-物质个数为car文件名
            elseif formatout==3
                datanamecar=strcat(species{1},'-',num2str(frame{1,irow}),'-',num2str(frame{2,irow}),'-base',num2str(base),'.pdb');%以目标物质-实践步-物质个数为car文件名
            end
            xyz_car_pdb_filemaker
            seekBOinform
        else
            fclose(rawdatatrj);
            rawdatatrj=fopen(datanametrj,'r');
        end
    end
    
    
    
    
    
elseif choi==1 && strcmpi(promptans,'n')
    rowcopy=0;[~,lenframe]=size(frame);
    rawdatatrj=fopen(datanametrj,'r');
    for irow=2:lenframe
        rowcopy=rowcopy+1;
        fprintf('\n\nchemi_mechanism is running, please wait...');
        fprintf('\ntrajectory frame(s) in total. Now group %d is being processing\n',prompt,rowcopy);
        tartrajectory=frame{1,irow};
        
        if exist('fram_num_check','var') && strcmpi(rerun_ans,'y') && rerun_ans2==1
            clear fram_num_check
            species_bonds_lammpstrj_framecheck
            rerun_ans2=rerun_ans2+1;
        elseif exist('fram_num_check','var') && strcmpi(rerun_ans,'n')
            fprintf('The frame No. from species.* file (outputnew) has been treated and reused\n')
        elseif ~exist('fram_num_check','var')
            species_bonds_lammpstrj_framecheck
            rerun_ans2=rerun_ans2+1;
        end
        if rerun_ans2==1
            species_bonds_lammpstrj_framecheck
            rerun_ans2=rerun_ans2+1;
        end
        bonds_analysis_speedup
        bondorder_deepmining
        
        [~,col]=size(tarelenummatch);bondset=[];
        for j=1:col/2
            molecomp='';
            for k=1:length(elementsequence)
                if tarelenummatch{k,2*j}==1
                    molecomp=strcat(molecomp,tarelenummatch{k,2*j-1});
                elseif tarelenummatch{k,2*j}>1
                    molecomp=strcat(molecomp,tarelenummatch{k,2*j-1},num2str(tarelenummatch{k,2*j}));
                end
            end
            if strcmp(molecomp,species{1})
                bondsetdata=length(bondset);bondset(bondsetdata+1)=j;
            end
        end
        bondsetdata=length(bondset);
        if bondsetdata~=outputdatanew{irow,4}
            if bondsetdata==0
                fprintf('\n\nWarning!!!No species is found through bonds file, but species file record it(number:%d).Possible reason:\nNevery parameter of "fix reax/c/bonds" in the in.* file is not match with the Nrepeat"fix reax/c/species"\n. Another possible reason element mapping or type error',outputdatanew{irow,4});
                fprintf('\nSolution:Nrepeat of "fix reax/c/species" in the in.* file should be 1, and Nevery(species)=Nfreq(species)=Nevery(bonds)\nReplace the mapping element to the expected one before data processing\n');
                errorexe=input('Ignore this problem and continue?y/n: \n','s');
				msgbox('Target species is not found! How to cope with this problem?');
                if strcmpi(errorexe,'y')
                    continue;
                elseif strcmpi(errorexe,'n')
                    error('Please see the aforesaid analysis and solution,please check it!!!');
                end
            else
                fprintf('\nCaution:target species number (%d) is not consistent with that of species recordation(%d)!!!\n',bondsetdata,outputdatanew{irow,4});
            end
        else
            fprintf('Target species number(%d) is consistent with that of species recordation(%d)',bondsetdata,outputdatanew{irow,4});
        end
        [tarraw,~]=size(tarBOinform);ii=1;jj=1;
        while tarraw
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
        fprintf('\nSuccessfully delete the unexpected BO information in tarBOinform\n\n');
        lammpstrj_analysis
        
        %unwrap
        if strcmpi(unwrapans,'y')
            load('BondRadii.mat');
            trjdata=PBC_Unwrap(tarBOinform,trjdata,BOXsize,boxsize,elementsequence,BondRadii);
        end
        if str2num(datarep)<=tartrajectory{1}
            if formatout==1
                datanamecar=strcat(species{1},'-',num2str(frame{1,irow}),'-',num2str(frame{2,irow}),'-base',num2str(base),'.xyz');%以目标物质-实践步-物质个数为car文件名
            elseif formatout==2
                datanamecar=strcat(species{1},'-',num2str(frame{1,irow}),'-',num2str(frame{2,irow}),'-base',num2str(base),'.car');%以目标物质-实践步-物质个数为car文件名
            elseif formatout==3
                datanamecar=strcat(species{1},'-',num2str(frame{1,irow}),'-',num2str(frame{2,irow}),'-base',num2str(base),'.pdb');%以目标物质-实践步-物质个数为car文件名
            end
            xyz_car_pdb_filemaker
            seekBOinform
        else
            fclose(rawdatatrj);
            rawdatatrj=fopen(datanametrj,'r');
        end
    end
    msgbox('Products data is successfully exported.');
    fclose(rawdatatrj);
    


    
    
    
elseif choi==2 || choi==4
    [~,outcol]=size(outputdatanew);
    for i=4:outcol
        if strcmpi(species{1},(outputdatanew{1,i}))
            outcol=i;
            break
        end
    end
    if choi==2
        fprintf('\nTarget species %s amount to %d in total (original number %d)',species{1},outputdatanew{floor((promptans3(1)-outputdatanew{2,1})/trajper(1)+2),outcol},outputdatanew{2,outcol});
        numstop=input('\nLimit the exportation of reactants-products? If yes, program will be forced to terminate. y/n: \n','s');
        if strcmpi(numstop,'y')
            fprintf('\nPlease input limited file number, <%d',outputdatanew{floor((promptans3(1)-outputdatanew{2,1})/trajper(1)+2),outcol})
            numstop=input(':\n');
        else
            numstop=outcol;
        end
    elseif choi==4
        fprintf('\nSince not all the species %d are depleted, please give the limited file number according to the species data, \nprogram will be forced to terminate.',outputdatanew{floor((promptans3(1)-outputdatanew{2,1})/trajper(1)+2),outcol},min(cell2mat(outputdatanew(2:end,outcol))),outputdatanew{end,outcol});
        numstop=input('\nPlease input the limited exported file number:\n');
    end
    fprintf('\nFor high efficiency, only search the trajectory frame with changed species number');%提高效率，只搜索有数目变化的帧
    seekacc=input('\nLimit the trajectory frame to be searched? For method 2 the searched trajectory frame with degressive timestep.\nFor method 4 is on the contrary. y/n: \n','s');
    if ~strcmpi(seekacc,'y') && ~strcmpi(seekacc,'n')
        error('Illegal limitation method option. Please check it!!!');
    end
    
    tartrajectory=frame{1,1};
    if exist('fram_num_check','var') && strcmpi(rerun_ans,'y') && rerun_ans2==1
        clear fram_num_check
        species_bonds_lammpstrj_framecheck
        rerun_ans2=rerun_ans2+1;
    elseif exist('fram_num_check','var') && strcmpi(rerun_ans,'n')
        fprintf('The frame No. from species.* file (outputnew) has been treated and reused\n')
        species_bonds_lammpstrj_framecheck
        rerun_ans2=rerun_ans2+1;
    elseif ~exist('fram_num_check','var')
        species_bonds_lammpstrj_framecheck
        rerun_ans2=rerun_ans2+1;
    end
    if rerun_ans2==1
        species_bonds_lammpstrj_framecheck
        rerun_ans2=rerun_ans2+1;
    end
    
    fprintf('\nAbstract BO information of target trajectory\n')
    bonds_analysis_speedup
    bondorder_deepmining
    tarBOinformfullcopy=tarBOinform;
    [~,col]=size(tarelenummatch);bondset=[];
    for j=1:col/2
        molecomp='';
        for k=1:length(elementsequence)
            if tarelenummatch{k,2*j}==1
                molecomp=strcat(molecomp,tarelenummatch{k,2*j-1});
            elseif tarelenummatch{k,2*j}>1
                molecomp=strcat(molecomp,tarelenummatch{k,2*j-1},num2str(tarelenummatch{k,2*j}));
            end
        end
        if strcmp(molecomp,species{1})
            bondsetdata=length(bondset);bondset(bondsetdata+1)=j;
        end
    end
    bondsetdata=length(bondset);
    if bondsetdata~=frame{2,1}
        if bondsetdata==0
            fprintf('\n\nWarning!!!No species is found through bonds file, but species file record it(number:%d).Possible reason:\nNevery parameter of "fix reax/c/bonds" in the in.* file is not match with the Nrepeat"fix reax/c/species"\n.Another possible reason element mapping or type error',frame{2,1});
            fprintf('\nSolution:Nrepeat of "fix reax/c/species" in the in.* file should be 1, and Nevery(species)=Nfreq(species)=Nevery(bonds)\nReplace the mapping element to the expected one before data processing\n');
            msgbox('Target species is not found!How to cope with this problem?');
			errorexe=input('Ignore this problem and continue?y/n: \n','s');
            if strcmpi(errorexe,'y')
%                 continue
            elseif strcmpi(errorexe,'n')
                error('Please see the aforesaid analysis and solution,please check it!!!');
            end
        else
            fprintf('\n\nCaution:target species number (%d) is not consistent with that of species recordation(%d)!!!\n',bondsetdata,frame{2,1});
        end
    else
        fprintf('Target species number(%d) is consistent with that of species recordation(%d)',bondsetdata,frame{2,1});
    end
    
    [tarraw,~]=size(tarBOinform);ii=1;jj=1;
    while tarraw
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
    fprintf('\nSuccessfully delete the unexpected BO information in tarBOinform\n');
    rawdatatrj=fopen(datanametrj,'r');
    lammpstrj_analysis
    %unwrap
    if strcmpi(unwrapans,'y')
        load('BondRadii.mat');
        trjdata=PBC_Unwrap(tarBOinform,trjdata,BOXsize,boxsize,elementsequence,BondRadii);
    end
    fclose(rawdatatrj);
    trjdatacopy=trjdata;
    tarBOinformcopy=tarBOinform;
    tartrajectorycopy=tartrajectory{1};
    
    m=0;productnum=[];
    [row5,~]=size(tarBOinformfullcopy);
    for i=1:row5
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
    
    [tarraw,~]=size(tarBOinformcopy);loopnum=1;pairnum=0;
    loop=1;trajectorynote=[];trajectorynote(1,1)=tartrajectorycopy(1);trajectorynote(1,2)=tartrajectorycopy(1);%存放合适帧的时间步数
    while tarraw
        if choi==2
            if strcmpi(seekacc,'n')
                tartrajectoryact=tartrajectorycopy(1)-loopnum*trajper(1);
                tartrajectoryact = FrameNoCertification (tartrajectoryact,trajper,choi,outputdatanew);
            elseif strcmpi(seekacc,'y')
                control=1;
                while control
                    tartrajectoryact=tartrajectorycopy(1)-loop*trajper(1);
                    tartrajectoryact = FrameNoCertification (tartrajectoryact,trajper,choi,outputdatanew);
                    loop=loop+1;
                    if tartrajectoryact/trajper(1)>=1 && outputdatanew{(tartrajectoryact-outputdatanew{2,1})/trajper(1)+2,4}<outputdatanew{(trajectorynote(1,2)-outputdatanew{2,1})/trajper(1)+2,4}
                        fprintf('\nNew frame is found %d, species number %d is less than that of last frame %d, which has %d',outputdatanew{(tartrajectoryact-outputdatanew{2,1})/trajper(1)+2,1},outputdatanew{(tartrajectoryact-outputdatanew{2,1})/trajper(1)+2,4},outputdatanew{(trajectorynote(1,2)-outputdatanew{2,1})/trajper(1)+2,1},outputdatanew{(trajectorynote(1,2)-outputdatanew{2,1})/trajper(1)+2,4});
                        trajectorynote(1,1)=trajectorynote(1,2);trajectorynote(1,2)=tartrajectoryact;
                        control=0;
                        break;
                    elseif tartrajectoryact/trajper(1)<1
                        error('Frame %d exceed the first timestep, meaning nonexistence, please check it!!!',tartrajectoryact);
                    end
                end
            end
            
        elseif choi==4
            if strcmpi(seekacc,'n')
                tartrajectoryact=tartrajectorycopy(1)+loopnum*trajper(1);
                tartrajectoryact = FrameNoCertification (tartrajectoryact,trajper,choi,outputdatanew);
            elseif strcmpi(seekacc,'y')
                [row,~]=size(outputdatanew);
                control=1;
                while control
                    tartrajectoryact=tartrajectorycopy(1)+loop*trajper(1);
                    tartrajectoryact = FrameNoCertification (tartrajectoryact,trajper,choi,outputdatanew);
                    loop=loop+1;
                    if tartrajectoryact<=outputdatanew{row,1} && outputdatanew{(tartrajectoryact-outputdatanew{2,1})/trajper(1)+2,4}<outputdatanew{(trajectorynote(1,2)-outputdatanew{2,1})/trajper(1)+2,4}
                        fprintf('\nNew frame is found %d, species number %d is less than that of last frame %d, which has %d',outputdatanew{(tartrajectoryact-outputdatanew{2,1})/trajper(1)+2,1},outputdatanew{(tartrajectoryact-outputdatanew{2,1})/trajper(1)+2,4},outputdatanew{(trajectorynote(1,2)-outputdatanew{2,1})/trajper(1)+2,1},outputdatanew{(trajectorynote(1,2)-outputdatanew{2,1})/trajper(1)+2,4});
                        trajectorynote(1,1)=trajectorynote(1,2);trajectorynote(1,2)=tartrajectoryact;
                        control=0;
                        break;
                    elseif tartrajectoryact>outputdatanew{row,1}
                        error('Frame %d exceed the last timestep, meaning nonexistence, please check it!!!',tartrajectoryact);
                    end
                end
            end
            [row,~]=size(outputdatanew);
            if tartrajectoryact>outputdatanew{row,1}
                fprintf('\nExceeding the last timestep, termination of the procedure.')
                break;
            end
        end
        
        if tartrajectoryact/trajper(1)>=1
            tartrajectory=tartrajectoryact;
            if exist('fram_num_check','var') && strcmpi(rerun_ans,'y') && rerun_ans2==1
                clear fram_num_check
                species_bonds_lammpstrj_framecheck
                rerun_ans2=rerun_ans2+1;
            elseif exist('fram_num_check','var') && strcmpi(rerun_ans,'n')
                fprintf('The frame No. from species.* file (outputnew) has been treated and reused\n')
            elseif ~exist('fram_num_check','var')
                species_bonds_lammpstrj_framecheck
                rerun_ans2=rerun_ans2+1;
            end
            if rerun_ans2==1
                species_bonds_lammpstrj_framecheck
                rerun_ans2=rerun_ans2+1;
            end
            
            bonds_analysis_speedup
            bondorder_deepmining
            
            frameproact=[];
            frameproact=react_blocklocate(tarBOinformcopy,tarBOinform,tartrajectory);
            
            if sum(frameproact(:,5))~=0
                if choi==2
                    fprintf('\n\n\nNow is searching group %d reactants(%d)-products(%d), this is the hit group, already hits %d groups(including this group),\n%d groups in total to be processed, limited to export %d groups',loopnum,tartrajectory{1,1},tartrajectorycopy(1),pairnum+1,outputdatanew{floor((promptans3(1)-outputdatanew{2,1})/trajper(1)+2),outcol},numstop);
                elseif choi==4
                    fprintf('\n\n\nNow is searching group %d products(%d)-reactants(%d),this is the hit group,already hits %d groups(including this group), \nat most %d froups, limited to export %d groups',loopnum,tartrajectory{1,1},tartrajectorycopy(1),pairnum+1,outputdatanew{floor((promptans3(1)-outputdatanew{2,1})/trajper(1)+2),outcol},numstop);
                end
                loopnum=loopnum+1;
                
                tarBOinformcopy3={};
                
                frameproact=react_blocklocate(tarBOinformcopy,tarBOinform,tartrajectory);
                [lenframe,~]=size(frameproact);pronum=0;
                for i=1:lenframe
                    if frameproact(i,5)~=0
                        pronum=pronum+1;
                        tarBOinformcopy(frameproact(i,2):frameproact(i,3),:)=[];
                        for j=i+1:lenframe
                            frameproact(j,2)=frameproact(j,2)-frameproact(j,4);
                            frameproact(j,3)=frameproact(j,3)-frameproact(j,4);
                        end
                        [~,colframe]=size(frameproact);[numzero,~]=ismember(frameproact(i,:),0);
                        numreactant=(colframe-5-sum(numzero))/2;
                        if numreactant>0
                            for j=1:numreactant
                                [row3,~]=size(tarBOinformcopy3);
                                tarBOinformcopy3(row3+1:row3+frameproact(i,2*j+5)-frameproact(i,2*j+4)+1,:)=tarBOinform(frameproact(i,2*j+4):frameproact(i,2*j+5),:);
                            end
                        else
                            if choi==2
                                error('Hit some reactants, but without block No.Please check it')
                            elseif choi==4
                                error('Hit some reactants, but without block No.Please check it')
                            end
                        end
                    end
                end
                
                [row4,~]=size(tarBOinformcopy3);[row5,~]=size(tarBOinformfullcopy);
                blockNO=[];
                for k=1:row4
                    m=1;
                    for ii=1:row5
                        if strcmp(tarBOinformfullcopy{ii,1},'#')
                            m=m+1;
                        end
                        if ~strcmp(tarBOinformcopy3{k,1},'#') && ~strcmp(tarBOinformfullcopy{ii,1},'#') && tarBOinformcopy3{k,1}==tarBOinformfullcopy{ii,1}
                            if ~ismember(m,blockNO)
                                lenbloc=length(blockNO);blockNO(lenbloc+1)=m;
                            end
                        end
                    end
                end
                blockNO=sort(blockNO);
                
                
                tarBOinformcopy2={};
                for i=1:length(blockNO)%
                    [row6,~]=size(tarBOinformcopy2);
                    tarBOinformcopy2(row6+1:row6+productnum(blockNO(i),3),:)=tarBOinformfullcopy(productnum(blockNO(i),1):productnum(blockNO(i),2),:);
                end
                
                
                rawdatatrj=fopen(datanametrj,'r');
                lammpstrj_analysis
                 %unwrap
                 if strcmpi(unwrapans,'y')
                     load('BondRadii.mat');
                     trjdata=PBC_Unwrap(tarBOinform,trjdata,BOXsize,boxsize,elementsequence,BondRadii);
                 end
                fclose(rawdatatrj);
                if choi==2
                    if formatout==1
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','reactant','-base',num2str(base),'.xyz');
                    elseif formatout==2
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','reactant','-base',num2str(base),'.car');
                    elseif formatout==3
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','reactant','-base',num2str(base),'.pdb');
                    end
                elseif choi==4
                    if formatout==1
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','product','-base',num2str(base),'.xyz');
                    elseif formatout==2
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','product','-base',num2str(base),'.car');
                    elseif formatout==3
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','product','-base',num2str(base),'.pdb');
                    end
                end
                tarBOinform=tarBOinformcopy3;
                xyz_car_pdb_filemaker
                
                if choi==2
                    if formatout==1
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','product','-base',num2str(base),'.xyz');
                    elseif formatout==2
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','product','-base',num2str(base),'.car');
                    elseif formatout==3
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','product','-base',num2str(base),'.pdb');
                    end
                elseif choi==4
                    if formatout==1
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','reactant','-base',num2str(base),'.xyz');
                    elseif formatout==2
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','reactant','-base',num2str(base),'.car');
                    elseif formatout==3
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','reactant','-base',num2str(base),'.pdb');
                    end
                end
                tarBOinform=tarBOinformcopy2;trjdata=trjdatacopy;
                xyz_car_pdb_filemaker
                
                seekBOinform
                
                pairnum=pairnum+1;
            else
                if choi==2
                    fprintf('\n\n\nGroup %d reactants(%d)-products(%d) hit nothing,already hit %d groups, continue to search the next frame of reactants\n',loopnum,tartrajectory{1,1},tartrajectorycopy(1),pairnum);
                elseif choi==4
                    fprintf('\n\n\nGroup %d products(%d)-reactants(%d) hit nothing, already hit %d groups, continue to search the next frame of products\n',loopnum,tartrajectory{1,1},tartrajectorycopy(1),pairnum);
                end
                loopnum=loopnum+1;
            end  
        else
            error('Frame exceed the first timestep, meaning nonexistence, please check it!!!');
        end
        [tarraw,~]=size(tarBOinformcopy);
        if pairnum==numstop
            tarraw=0;
		elseif pairnum>numstop
			fprintf('Warning!!!Hit group is larger than the limited trajectory frame(s), please check it!!!');
        end
    end
    if choi==2
        inspectBOinform
        fprintf('\nTarget trajectory frame %d: products related reactants-products is finished, generating %d groups of reactants-products files in total',tartrajectorycopy(1),pairnum);
    elseif choi==4
        inspectBOinform
        fprintf('\nTarget trajectory frame %d reactants related reactants-products is finished, generating %d groups of reactants-products files in total',tartrajectorycopy(1),pairnum);
    end
    
    msgbox('Reactants-products files are successfully analyzed and generated.');
   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
elseif choi==3
    fprintf('\nWork is on the go and it will be on line soon!\n');
    fprintf('\n1.Cheer!');
    
    
    
    
    
    
    
    
else
    disp('Illegal input, please check the subprogram No.Please check it!!!')
end
fprintf('chemi_mechanism is successfully finished');
msgbox('chemi_mechanism is successfully finished');

Elapsedtime = toc;
fprintf('\nTotal task time: %.2f s\n',Elapsedtime)

clear blockNO bondset bondsetdata BOXsize choi colframe datanamebond datanametrj date elementsequence fidBO fileheader frame frameproact
clear ii irow jj lenbloc lenframe loopnum m molecomp numreactant numstop numzero outcol pairnum PBC PBCa PBCalpha PBCb PBCbeta PBCc
clear PBCchoi PBCgamma productnum prompt promptans promptans2 promptans3 pronum rawdata rawdatatrj reactname residueseqname
clear row3 row4 row5 row6 rowcopy rowfram spacegroupname species tarelenummatch tarraw tartrajectory tartrajectoryact tartrajectorycopy
clear tartrjdata title trajper loop datanamespe i j seekacc sperunans trajectorynote producttname col maxchoi maxdata outputdata_copy
clear atomid_conv base boxsize counter datanamecar eleswap eleswapans modcounter moddivid formatout datarep errorexe fram_num_check
clear rerun_ans sperunans2 outputdatanew_frame rerun_ans2 species_frame_checkans unwrapans species_frame_check2 readline MOLE
clear iijj elemax datalinenum datacellnum control checknum check_control_origin coord_position coord_tag outputdatanew row species_copy
clear tarBOinformcopy2 tarBOinformcopy3 tarBOinformfullcopy trjdata trjdatacopy tarBOinformcopy ans
