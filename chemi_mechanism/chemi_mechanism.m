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
        fprintf('\n根据原子个数采用不同进制（ASCII编码）与元素名合成原子名称，总长度不大于4个字符，元素名最大含两个字符，则原子id应不大于3个字符\n');
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
base=input('\n请输入采用的进制数，必须为正整数且不小于推荐的进制数：\n');
unwrapans=input('\n是否执行坐标unwrap处理跨越周期性盒子之外的ghost坐标对可视化的影响,y/n?: \n','s');

if choi==1
    maxchoi=input('\n请输入帧阈值，超过该值将会提示按不同方法进行部分抽样导出,建议值为10：\n');
end
if choi==1 || choi==2 || choi==4
    if formatout==2
        PBCchoi=input('\n请输入周期性边界条件开关选项，ON/OFF：\n','s');PBCchoi=upper(PBCchoi);
        if strcmp(PBCchoi,'ON')%写入file header信息
            PBC='PBC=ON';
            disp('周期性边界条件示例“PBC   33.4531   33.4531   33.4531   90.0000   90.0000   90.0000 (P1)”');
            PBCalpha=input('请输入周期性边界条件alpha，可由MS中lattice parameters得到，四位小数：\n');
            disp('周期性边界条件示例“PBC   33.4531   33.4531   33.4531   90.0000   90.0000   90.0000 (P1)”');
            PBCbeta=input('请输入周期性边界条件beta,可由MS中lattice parameters得到，四位小数：\n');
            disp('周期性边界条件示例“PBC   33.4531   33.4531   33.4531   90.0000   90.0000   90.0000 (P1)”');
            PBCgamma=input('请输入周期性边界条件gamma,可由MS中lattice parameters得到，四位小数：\n');
            disp('周期性边界条件示例“PBC   33.4531   33.4531   33.4531   90.0000   90.0000   90.0000 (P1)”');
            spacegroupname=input('请在括号里输入周期性边界条件group name,示例“（P1）”：\n','s');spacegroupname=upper(spacegroupname);
        elseif strcmp(PBCchoi,'OFF')
            PBC='PBC=OFF';
        else
            disp('PBCchoi周期性边界条件输入非法，请检查');
            return;
        end
    elseif formatout==3
        PBCchoi=input('\n请输入周期性边界条件开关选项，ON/OFF：\n','s');PBCchoi=upper(PBCchoi);
        if strcmp(PBCchoi,'ON')
            PBCalpha=input('\n请输入周期性边界条件alpha，可由MS中lattice parameters得到，四位小数：\n');
            PBCbeta=input('\n请输入周期性边界条件beta,可由MS中lattice parameters得到，四位小数：\n');
            PBCgamma=input('\n请输入周期性边界条件gamma,可由MS中lattice parameters得到，四位小数：\n');
            spacegroupname=input('\n请在括号里输入周期性边界条件group name,示例对于*.arc是“（P1）”，而*.pdb是“P 1”：\n','s');spacegroupname=upper(spacegroupname);
        elseif strcmp(PBCchoi,'OFF')
            warndlg('非周期性约束：坐标提取将不涉及周期性盒子尺寸，直接提取！')
        end
    end
    fprintf('\n强烈建议不要缩放坐标（使用dump_modify scale no)，容易导致计算误差，增大成图的瑕疵')
    BOXsize=input('\n请问*.lammpstrj中原子坐标标是否是缩放的（dump_modify scale yes(no)缩放(不缩放），y/n:\n','s');BOXsize=lower(BOXsize);
    if ~ismember(BOXsize,{'y','n'})
        error('非法的BOXsize参数');
    end
    
    if exist('outputdatanew','var')
        fprintf('\n警告：species_analysis_capture程序未运行，因为工作空间存在outputdatanew变量,请确保是当前考察文件的导入数据');
        msgbox('警告：species_analysis_capture程序未运行,因为工作空间存在outputdatanew变量,是否重新运行species_analysis_capture程序导入产物？');
        sperunans=input('\nspecies_analysis_capture程序未运行,因为工作空间存在outputdatanew变量,是否重新运行species_analysis_capture程序导入产物？y/n：\n','s');
        if strcmpi(sperunans,'y')
            if strcmpi(rerun_ans,'n')%不重新检查帧编号情况，继续沿用
                sperunans2=input('是否重新修改已存在的species帧编号(outputnew中)，同一文件若帧编号不一致已修改则继续研究其他物质时不需要修改！y/n：\n','s');%同一文件反复运行，存在species帧编号时不会每次修正而耗时
            else
                sperunans2='y';
            end
            species_analysis_capture %导入目标species文件
            msgbox('species_analysis_capture运行结束');
        elseif strcmpi(sperunans,'n')
            fprintf('\n继续使用工作空间已有的存储了产物信息的outputdatanew变量');
            species=strtrim(species);
            species=strsplit(species);
        end
    end
    if ~exist('outputdatanew','var')
        species_analysis_capture %导入目标species文件
    end
    %捕获目标分子式数量变化信息，用于确定分析起讫帧
    [row,~]=size(outputdatanew);
    residueseqname=species{1};%残基名确定为目标分子式
    
    tic %开始计时
    
    prompt=0;rowcopy=row-1;irow=2;%建立提示进度,找到总共需要导出的数目
    while rowcopy
        if (irow==2 && outputdatanew{irow,4}~=0) ||(irow>2 && outputdatanew{irow,4}~=0 && outputdatanew{irow,4}~=outputdatanew{irow-1,4})
            prompt=prompt+1;
        end
        rowcopy=rowcopy-1;irow=irow+1;
    end
    
    if choi==1
        if prompt>maxchoi
            fprintf('\n在outputdatanew中总共有%d帧满足导出要求,文件数目大于帧阈值，建议导出部分图像',prompt);
            promptans=input('\n是否导出部分图像,y/n：\n','s');
            if strcmpi(promptans,'y')
                promptans2=input('请选择确定导出帧方法数字代码：\n1.单调递增的帧\n2.限制帧数的等差法（含首尾两帧，相邻帧数量不同）\n3.限制单调递增的帧的等差法（含首尾两帧）在\n4.手动输入时间步的帧\n');
            end
        else
            promptans='n';
        end
    elseif choi==2 || choi==4
        promptans='y';promptans2=4;
    end
    
    
    if strcmpi(promptans,'y')%按所选方法建立到导出帧表导出较少的帧
        frame={};%存储满足所选方法的帧的时间步
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
            promptans3=input('\n请输入需要导出的帧数,大于2：\n');
            counter=0;%配合总数目prompt,计数导入的帧步数
            modcounter=mod(prompt-2,promptans3);
            moddivid=(prompt-2-modcounter)/(promptans3-2);
            if prompt<promptans3
                error('可导出的帧数少于要求的帧数');
            end
            for irow=2:row
                if (irow==2 && outputdatanew{irow,4}~=0) ||(irow>2 && outputdatanew{irow,4}~=0 && outputdatanew{irow,4}~=outputdatanew{irow-1,4})%只导出有产物数量变化的帧
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
                        counter=counter+1;%统计满足要求的单调递增的帧数
                    end
                end
            end
            countercopy=counter;
            fprintf('\n满足要求的共有%d帧',countercopy);
            promptans3=input('\n请输入需要导出的帧数：\n');
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
                promptans3=input('\n请输入想要导出的帧数的时间步，用空号隔开：\n','s');
            
            elseif choi==2 || choi==4
                fprintf('可通过species_analysis，species_capture先确定想要研究的目标产物前后反应的特定时间步');
                msgbox('species分析完毕，请输入产物所在的帧数的时间步');
                promptans3=input('\n请输入想要导出的产物所在的帧数的时间步，用空号隔开：\n','s');
            end
            
            promptans3=strtrim(promptans3);promptans3=strsplit(promptans3);promptans3=str2double(promptans3);
            frame=num2cell(promptans3);lenframe=length(frame);
            tartrajectory=promptans3;
            
            if exist('fram_num_check','var') && strcmpi(rerun_ans,'y') && rerun_ans2==1
                clear fram_num_check
                species_bonds_lammpstrj_framecheck%检查species文件帧编号(outputnew中)与bonds及files文件是否一致，不一致时修改
                rerun_ans2=rerun_ans2+1;%之后循环导出不再检查
            elseif exist('fram_num_check','var') && strcmpi(rerun_ans,'n')
                fprintf('species帧（outputnew中）编号已处理，继续沿用\n')
            elseif ~exist('fram_num_check','var')
                species_bonds_lammpstrj_framecheck
                rerun_ans2=rerun_ans2+1;%之后循环导出不再检查 
            end
            if rerun_ans2==1
                species_bonds_lammpstrj_framecheck
                rerun_ans2=rerun_ans2+1;
            end
            frame=num2cell(promptans3);%重新赋值可能修复后的帧编号
            
            for irow=1:lenframe
                if ~ismember(frame{irow},cell2mat(outputdatanew(2:end,1)))
                    error('手动输入帧至少有一帧不存在，请检查！！！');
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
            fprintf('非法的导出帧方法数字代码输入，请检查！！！');
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
    rawdatatrj=fopen(datanametrj,'r');%采用逐行读取处理lammpstrj文件
    for irow=1:lenframe%循环导出部分抽样的图
        rowcopy=rowcopy+1;
        fprintf('\n\nchemi_mechanism程序运行中，请等待...');
        fprintf('\n总共有%d帧，当前生成第%d帧car文件',lenframe,rowcopy);
        tartrajectory={};
        tartrajectory=frame{1,irow};%目标键级和坐标轨迹
            
        if exist('fram_num_check','var') && strcmpi(rerun_ans,'y') && rerun_ans2==1
            clear fram_num_check
            species_bonds_lammpstrj_framecheck %检查species文件帧编号(outputnew中)与bonds及files文件是否一致，不一致时修改
            rerun_ans2=rerun_ans2+1;%之后循环导出不再检查
        elseif exist('fram_num_check','var') && strcmpi(rerun_ans,'n')
            fprintf('species帧（outputnew中）编号已处理，继续沿用\n')
        elseif ~exist('fram_num_check','var')
            species_bonds_lammpstrj_framecheck
            rerun_ans2=rerun_ans2+1;%之后循环导出不再检查
        end
        if rerun_ans2==1
            species_bonds_lammpstrj_framecheck
            rerun_ans2=rerun_ans2+1;
        end
        
        bonds_analysis_speedup%提取目标轨迹键级信息
        bondorder_deepmining%按分子式归类各物质，便于输出目标分子式图形
        
        [~,col]=size(tarelenummatch);bondset=[];
        for j=1:col/2%组装tarelenummatch中分子式，与目标分子式比对
            molecomp='';
            for k=1:length(elementsequence)
                if tarelenummatch{k,2*j}==1
                    molecomp=strcat(molecomp,tarelenummatch{k,2*j-1});
                elseif tarelenummatch{k,2*j}>1
                    molecomp=strcat(molecomp,tarelenummatch{k,2*j-1},num2str(tarelenummatch{k,2*j}));
                end
            end
            if strcmp(molecomp,species{1})%找到目标分子所在的第i块分区，tarBOinform中'#'分割线之上
                bondsetdata=length(bondset);bondset(bondsetdata+1)=j;
            end
        end
        bondsetdata=length(bondset);
        if bondsetdata~=frame{2,irow}%发现species文件记录的分子数与bonds统计的有些数据居然对不上！！！
            if bondsetdata==0
                fprintf('\n\n警告！！！通过bonds文件未匹配到目标产物，但是species文件记录了该物质(数量为%d)，原因可能\n是in文件中fix reax/c/bonds为特定Nevery输出键级等信息，而fix reax/c/species可能设置为Nevery取样键级Nrepeat平均\n并Nfreq输出此平均键级产物。另一原因是同种原子映射为不停的原子模拟或原子type输错',frame{2,irow});
                fprintf('\n解决方法：设置in文件中fix reax/c/species中Nrepeat为1，且有Nevery(species)=Nfreq(species)=Nevery(bonds),完全匹配二者输出信息\n替换转换为后者本来的类型。\n重新检查输入原子type\n');
                msgbox('未匹配到目标产物,请检查元素顺序是否正确，请选择处理方式！');
				errorexe=input('是否忽略错误？y/n:\n','s');
				msgbox('程序检测到错误存在，请处理！');
                if strcmpi(errorexe,'y')
                    continue;%忽略本帧
                elseif strcmpi(errorexe,'n')
                    error('请参见上方错误分析，请检查！！！');
                end
            else
                fprintf('\n请注意:bonds文件中匹配到的目标分子数（%d)与species文件中数量(%d)不一致!!!\n',bondsetdata,frame{2,irow});
            end
        else
            fprintf('搜索到目标产物数量（%d)，数量与species文件中数量(%d)一致',bondsetdata,frame{2,irow});
        end
        [tarraw,~]=size(tarBOinform);ii=1;jj=1;
        while tarraw%将tarBOinform中非目标产物的键级信息删除
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
        fprintf('\n成功剔除tarBOinform中非目标产物的键级信息\n');
        lammpstrj_analysis%提取目标轨迹坐标信息
        %输出该轨迹中所有可表示为目标分子式的产物car文件，其在tarBOinform中的键级信息区块记录在bondset中，坐标信息在trjdata中
        %对坐标进行unwrap
        if strcmpi(unwrapans,'y')
            load('BondRadii.mat');
            trjdata=PBC_Unwrap(tarBOinform,trjdata,BOXsize,boxsize,element,BondRadii);
        end
        if str2num(datarep)<=tartrajectory{1}%dump.*文件帧确保是需要导出的帧（dump.*并非由0时间步开始记录导致错误）,小于最小帧（不存在）则无法写文件
            if formatout==1
                datanamecar=strcat(species{1},'-',num2str(frame{1,irow}),'-',num2str(frame{2,irow}),'-base',num2str(base),'.xyz');%以目标物质-实践步-物质个数为car文件名
            elseif formatout==2
                datanamecar=strcat(species{1},'-',num2str(frame{1,irow}),'-',num2str(frame{2,irow}),'-base',num2str(base),'.car');%以目标物质-实践步-物质个数为car文件名
            elseif formatout==3
                datanamecar=strcat(species{1},'-',num2str(frame{1,irow}),'-',num2str(frame{2,irow}),'-base',num2str(base),'.pdb');%以目标物质-实践步-物质个数为car文件名
            end
            xyz_car_pdb_filemaker%制作xyz_car文件
            seekBOinform%写出产物键级信息
        else
            fclose(rawdatatrj);%关闭lammpstrj文件
            rawdatatrj=fopen(datanametrj,'r');%采用逐行读取处理lammpstrj文件,重新打开
        end
    end
    
    
    
    
    
elseif choi==1 && strcmpi(promptans,'n')%导出所有相邻不等的帧，很多！！！
    rowcopy=0;[~,lenframe]=size(frame);%节约变量名，拼了
    rawdatatrj=fopen(datanametrj,'r');%采用逐行读取处理lammpstrj文件
    for irow=2:lenframe%本段不可改irow,循环导出
        rowcopy=rowcopy+1;
        fprintf('\n\nchemi_mechanism程序运行中，请等待...');
        fprintf('\n总共有%d帧，当前生成第%d帧car文件\n',prompt,rowcopy);
        tartrajectory=frame{1,irow};%目标键级和坐标轨迹
        
        if exist('fram_num_check','var') && strcmpi(rerun_ans,'y') && rerun_ans2==1
            clear fram_num_check
            species_bonds_lammpstrj_framecheck%检查species文件帧编号(outputnew中)与bonds及files文件是否一致，不一致时修改
            rerun_ans2=rerun_ans2+1;%之后循环导出不再检查
        elseif exist('fram_num_check','var') && strcmpi(rerun_ans,'n')
            fprintf('species帧（outputnew中）编号已处理，继续沿用\n')
        elseif ~exist('fram_num_check','var')
            species_bonds_lammpstrj_framecheck
            rerun_ans2=rerun_ans2+1;%之后循环导出不再检查
        end
        if rerun_ans2==1
            species_bonds_lammpstrj_framecheck
            rerun_ans2=rerun_ans2+1;
        end
        bonds_analysis_speedup%提取目标轨迹键级信息
        bondorder_deepmining%按分子式归类各物质，便于输出目标分子式图形
        
        [~,col]=size(tarelenummatch);bondset=[];
        for j=1:col/2%组装tarelenummatch中分子式，与目标分子式比对
            molecomp='';
            for k=1:length(elementsequence)
                if tarelenummatch{k,2*j}==1
                    molecomp=strcat(molecomp,tarelenummatch{k,2*j-1});
                elseif tarelenummatch{k,2*j}>1
                    molecomp=strcat(molecomp,tarelenummatch{k,2*j-1},num2str(tarelenummatch{k,2*j}));
                end
            end
            if strcmp(molecomp,species{1})%找到目标分子所在的第i块分区，tarBOinform中'#'分割线之上
                bondsetdata=length(bondset);bondset(bondsetdata+1)=j;
            end
        end
        bondsetdata=length(bondset);
        if bondsetdata~=outputdatanew{irow,4}%发现species文件记录的分子数与bonds统计的有些数据居然对不上！！！
            if bondsetdata==0
                fprintf('\n\n警告！！！通过bonds文件未匹配到目标产物，但是species文件记录了该物质(数量为%d)，原因可能\n是in文件中fix reax/c/bonds为特定Nevery输出键级等信息，而fix reax/c/species可能设置为Nevery取样键级Nrepeat平均\n并Nfreq输出此平均键级产物。另一原因是同种原子映射为不停的原子模拟或原子type输错',outputdatanew{irow,4});
                fprintf('\n解决方法：设置in文件中fix reax/c/species中Nrepeat为1，且有Nevery(species)=Nfreq(species)=Nevery(bonds),完全匹配二者输出信息\n替换转换为后者本来的类型\n');
                errorexe=input('是否忽略错误？y/n:\n','s');
				msgbox('程序检测到错误存在，请处理！！！');
                if strcmpi(errorexe,'y')
                    continue;%忽略本帧
                elseif strcmpi(errorexe,'n')
                    error('请参见上方错误分析，请检查！！！');
                end
            else
                fprintf('\n请注意:bonds文件中匹配到的目标分子数（%d)与species文件中数量(%d)不一致!!!\n',bondsetdata,outputdatanew{irow,4});
            end
        else
            fprintf('搜索到目标产物数量（%d)，数量与species文件中数量(%d)一致',bondsetdata,outputdatanew{irow,4});
        end
        [tarraw,~]=size(tarBOinform);ii=1;jj=1;
        while tarraw%将tarBOinform中非目标产物的键级信息删除
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
        fprintf('\n成功剔除tarBOinform中非目标产物的键级信息\n');
        lammpstrj_analysis%提取目标轨迹坐标信息
        %输出该轨迹中所有可表示为目标分子式的产物car文件，其在tarBOinform中的键级信息区块记录在bondset中，坐标信息在trjdata中
        %对坐标进行unwrap
        if strcmpi(unwrapans,'y')
            load('BondRadii.mat');
            trjdata=PBC_Unwrap(tarBOinform,trjdata,BOXsize,boxsize,element,BondRadii);
        end
        if str2num(datarep)<=tartrajectory{1}%dump.*文件帧确保是需要导出的帧（dump.*并非由0时间步开始记录导致错误）,小于最小帧（不存在）则无法写文件
            if formatout==1
                datanamecar=strcat(species{1},'-',num2str(frame{1,irow}),'-',num2str(frame{2,irow}),'-base',num2str(base),'.xyz');%以目标物质-实践步-物质个数为car文件名
            elseif formatout==2
                datanamecar=strcat(species{1},'-',num2str(frame{1,irow}),'-',num2str(frame{2,irow}),'-base',num2str(base),'.car');%以目标物质-实践步-物质个数为car文件名
            elseif formatout==3
                datanamecar=strcat(species{1},'-',num2str(frame{1,irow}),'-',num2str(frame{2,irow}),'-base',num2str(base),'.pdb');%以目标物质-实践步-物质个数为car文件名
            end
            xyz_car_pdb_filemaker%制作xyz_car文件
            seekBOinform%写出产物键级信息
        else
            fclose(rawdatatrj);%关闭lammpstrj文件
            rawdatatrj=fopen(datanametrj,'r');%采用逐行读取处理lammpstrj文件,重新打开
        end
    end
    msgbox('产物car文件导出完成');
    fclose(rawdatatrj);%关闭lammpstrj文件
    


    
    
    
elseif choi==2 || choi==4%确定该产物数量变化的连续两帧时间步数，根据单一产物，关联的多产物或反应物，某个基团或者化学键
    [~,outcol]=size(outputdatanew);%控制太多时导出文件数目
    for i=4:outcol
        if strcmpi(species{1},(outputdatanew{1,i}))
            outcol=i;
            break
        end
    end
    if choi==2
        fprintf('\n目标物质%s物质总共有%d个（初始%d个）',species{1},outputdatanew{floor((promptans3(1)-outputdatanew{2,1})/trajper+2),outcol},outputdatanew{2,outcol});
        numstop=input('\n是否限制导出反应物-产物car文件对数目，如果限制，达到此数目将强行终止搜寻任务。y/n：\n','s');
        if strcmpi(numstop,'y')
            fprintf('\n请输入不大于%d的文件限制数目',outputdatanew{floor((promptans3(1)-outputdatanew{2,1})/trajper+2),outcol})
            numstop=input(':\n');
        else
            numstop=outcol;
        end
    elseif choi==4
        fprintf('\n由于反应物%d个不一定全部消失转化（最小剩余%d个，最终剩%d个），请根据产物导出数据确定car文件对数目，达到此数目将强行终止搜寻任务',outputdatanew{floor((promptans3(1)-outputdatanew{2,1})/trajper+2),outcol},min(cell2mat(outputdatanew(2:end,outcol))),outputdatanew{end,outcol});
        numstop=input('\n请输入文件限制数目:\n');
    end
    fprintf('\n为了提高搜索效率，可以只对有数目变化的帧进行搜索，尽管不变化也可能恰好由于生成与消耗相抵消，但是概率很小');%提高效率，只搜索有数目变化的帧
    seekacc=input('\n是否限制只搜索数目变化的帧？对于2号子程序将是递减反应物的帧，对于4子程序将是递增产物的帧。y/n:\n','s');
    if ~strcmpi(seekacc,'y') && ~strcmpi(seekacc,'n')
        error('非法的限制搜索数目变化的帧选项输入，请检查！！！');
    end
    
    tartrajectory=frame{1,1};
    if exist('fram_num_check','var') && strcmpi(rerun_ans,'y') && rerun_ans2==1
        clear fram_num_check
        species_bonds_lammpstrj_framecheck%检查species文件帧编号(outputnew中)与bonds及files文件是否一致，不一致时修改
        rerun_ans2=rerun_ans2+1;%之后循环导出不再检查
    elseif exist('fram_num_check','var') && strcmpi(rerun_ans,'n')
        fprintf('species帧（outputnew中）编号已处理，继续沿用\n')
        species_bonds_lammpstrj_framecheck%检查species文件帧编号(outputnew中)与bonds及files文件是否一致，不一致时修改
        rerun_ans2=rerun_ans2+1;%之后循环导出不再检查
    elseif ~exist('fram_num_check','var')
        species_bonds_lammpstrj_framecheck
        rerun_ans2=rerun_ans2+1;%之后循环导出不再检查
    end
    if rerun_ans2==1
        species_bonds_lammpstrj_framecheck
        rerun_ans2=rerun_ans2+1;
    end
    
    fprintf('\n提取目标轨迹键级信息\n')
    bonds_analysis_speedup%提取目标轨迹键级信息
    bondorder_deepmining%按分子式归类各物质，便于输出目标分子式用于后续比对与出图
    tarBOinformfullcopy=tarBOinform;%全复制，避免剔除非产物后影响反应物反射产物搜寻相应原子所在产物
    [~,col]=size(tarelenummatch);bondset=[];
    for j=1:col/2%组装tarelenummatch中分子式，与目标分子式比对
        molecomp='';
        for k=1:length(elementsequence)
            if tarelenummatch{k,2*j}==1
                molecomp=strcat(molecomp,tarelenummatch{k,2*j-1});
            elseif tarelenummatch{k,2*j}>1
                molecomp=strcat(molecomp,tarelenummatch{k,2*j-1},num2str(tarelenummatch{k,2*j}));
            end
        end
        if strcmp(molecomp,species{1})%找到目标分子所在的第i块分区，tarBOinform中'#'分割线之上
            bondsetdata=length(bondset);bondset(bondsetdata+1)=j;
        end
    end
    bondsetdata=length(bondset);
    if bondsetdata~=frame{2,1}%发现species文件记录的分子数与bonds统计的有些数据居然对不上！！！
        if bondsetdata==0
            fprintf('\n\n通过bonds文件未匹配到目标产物，但是species文件记录了该物质(数量为%d)，原因可能\n是in文件中fix reax/c/bonds为特定Nevery输出键级等信息，而fix reax/c/species可能设置为Nevery取样键级Nrepeat平均\n并Nfreq输出此平均键级产物。另一原因是同种原子映射为不同的原子模拟或原子type顺序输入错误',frame{2,1});
            fprintf('\n解决方法：设置in文件中fix reax/c/species中Nrepeat为1，且有Nevery(species)=Nfreq(species)=Nevery(bonds),完全匹配二者输出信息\n替换转换为后者本来的类型');
            msgbox('未匹配到目标产物,请检查元素顺序是否正确，请选择处理方式！');
			errorexe=input('\n是否忽略错误？y/n:\n','s');
            if strcmpi(errorexe,'y')
%                 continue;%忽略本帧
            elseif strcmpi(errorexe,'n')
                error('请参见上方错误分析，请检查！！！');
            end
        else
            fprintf('\n请注意:bonds文件中匹配到的目标分子数（%d)与species文件中数量(%d)不一致!!!\n',bondsetdata,frame{2,1});
        end
    else
        fprintf('搜索到目标产物数量（%d)，数量与species文件中数量(%d)一致',bondsetdata,frame{2,1});
    end
    
    [tarraw,~]=size(tarBOinform);ii=1;jj=1;
    while tarraw%将tarBOinform中非目标产物的键级信息删除
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
    fprintf('\n成功剔除tarBOinform中非目标产物的键级信息\n');
    rawdatatrj=fopen(datanametrj,'r');%采用逐行读取处理lammpstrj文件
    lammpstrj_analysis%提取产物轨迹坐标信息
    %对坐标进行unwrap
    if strcmpi(unwrapans,'y')
        load('BondRadii.mat');
        trjdata=PBC_Unwrap(tarBOinform,trjdata,BOXsize,boxsize,element,BondRadii);
    end
    fclose(rawdatatrj);%关闭lammpstrj文件
    trjdatacopy=trjdata;%复制目标产物坐标信息，避免分析前端反应物被覆盖
    tarBOinformcopy=tarBOinform;%复制剔除了非目标键级信息后的目标产物键级信息，避免分析前端反应物被覆盖
    tartrajectorycopy=tartrajectory{1};
    
    m=0;productnum=[];%记录命中到反应物的产物和产物所有原子在产物帧中各产物的键级信息，对应坐标在trjdatacopy中
    [row5,~]=size(tarBOinformfullcopy);
    for i=1:row5%记录产物帧中各物质起始行和行数
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
    
    [tarraw,~]=size(tarBOinformcopy);loopnum=1;pairnum=0;%向前回溯反应物帧
    loop=1;trajectorynote=[];trajectorynote(1,1)=tartrajectorycopy(1);trajectorynote(1,2)=tartrajectorycopy(1);%存放合适帧的时间步数
    while tarraw%寻找目标帧产物的反应物帧，反应物在不同时间步的分开制作xyz_car文件
        if choi==2
            if strcmpi(seekacc,'n')
                tartrajectoryact=tartrajectorycopy(1)-loopnum*trajper;
            elseif strcmpi(seekacc,'y')
                control=1;
                while control
                    tartrajectoryact=tartrajectorycopy(1)-loop*trajper;
                    loop=loop+1;
                    if tartrajectoryact/trajper>=1 && outputdatanew{(tartrajectoryact-outputdatanew{2,1})/trajper+2,4}<outputdatanew{(trajectorynote(1,2)-outputdatanew{2,1})/trajper+2,4}
                        fprintf('\n找到新的产物帧%d，其产物物数目%d比上一帧%d搜寻的数目%d少\n',outputdatanew{(tartrajectoryact-outputdatanew{2,1})/trajper+2,1},outputdatanew{(tartrajectoryact-outputdatanew{2,1})/trajper+2,4},outputdatanew{(trajectorynote(1,2)-outputdatanew{2,1})/trajper+2,1},outputdatanew{(trajectorynote(1,2)-outputdatanew{2,1})/trajper+2,4});
                        trajectorynote(1,1)=trajectorynote(1,2);trajectorynote(1,2)=tartrajectoryact;
                        control=0;
                        break;
                    elseif tartrajectoryact/trajper<1
                        error('反应物帧%d超出第一帧，此帧不存在，请检查！！！',tartrajectoryact);
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
                        fprintf('\n找到新的产物帧%d，其产物物数目%d比上一帧%d搜寻的数目%d少\n',outputdatanew{(tartrajectoryact-outputdatanew{2,1})/trajper+2,1},outputdatanew{(tartrajectoryact-outputdatanew{2,1})/trajper+2,4},outputdatanew{(trajectorynote(1,2)-outputdatanew{2,1})/trajper+2,1},outputdatanew{(trajectorynote(1,2)-outputdatanew{2,1})/trajper+2,4});
                        trajectorynote(1,1)=trajectorynote(1,2);trajectorynote(1,2)=tartrajectoryact;
                        control=0;
                        break;
                    elseif tartrajectoryact>outputdatanew{row,1}
                        error('反应物帧%d超出最大帧，此帧不存在，请检查！！！',tartrajectoryact);
                    end
                end
            end
            [row,~]=size(outputdatanew);
            if tartrajectoryact>outputdatanew{row,1}
                fprintf('\n超出数据的最大时间步，搜寻终止')
                break;
            end
        end
        
        if tartrajectoryact/trajper>=1
            tartrajectory=tartrajectoryact;
            if exist('fram_num_check','var') && strcmpi(rerun_ans,'y') && rerun_ans2==1
                clear fram_num_check
                species_bonds_lammpstrj_framecheck%检查species文件帧编号(outputnew中)与bonds及files文件是否一致，不一致时修改
                rerun_ans2=rerun_ans2+1;%之后循环导出不再检查
            elseif exist('fram_num_check','var') && strcmpi(rerun_ans,'n')
                fprintf('species帧（outputnew中）编号已处理，继续沿用\n')
            elseif ~exist('fram_num_check','var')
                species_bonds_lammpstrj_framecheck
                rerun_ans2=rerun_ans2+1;%之后循环导出不再检查
            end
            if rerun_ans2==1
                species_bonds_lammpstrj_framecheck
                rerun_ans2=rerun_ans2+1;
            end
            
            bonds_analysis_speedup%提取目标轨迹键级信息
            bondorder_deepmining%按分子式归类各物质键级信息，获得反应帧tarBOinform
            %将tarBOinformcopy中产物在反应物帧的tarBOinform中找到键级区块
            
            frameproact=[];%存放tarBOinformcopy中产物顺序与对应反应物所在帧时间步和反应物tarBOinform各键级块区
            frameproact=react_blocklocate(tarBOinformcopy,tarBOinform,tartrajectory);%得到该时间步反应物tarBOinform各键级块区
            
            if sum(frameproact(:,5))~=0%确保有反应物
                if choi==2
                    fprintf('\n\n\n正在搜索第%d组反应物(%d)-产物(%d)匹配组,该组为命中组，已命中%d组(含该组),总共%d组(不在同一反应物帧出现多个反应时),限制导出%d组',loopnum,tartrajectory{1,1},tartrajectorycopy(1),pairnum+1,outputdatanew{floor((promptans3(1)-outputdatanew{2,1})/trajper+2),outcol},numstop);
                elseif choi==4
                    fprintf('\n\n\n正在搜索第%d组产物(%d)-反产物(%d)匹配组,该组为命中组，已命中%d组(含该组),最多有%d组(在全部转化为别的产物时),限制导出%d组',loopnum,tartrajectory{1,1},tartrajectorycopy(1),pairnum+1,outputdatanew{floor((promptans3(1)-outputdatanew{2,1})/trajper+2),outcol},numstop);
                end
                loopnum=loopnum+1;
                %写出car文件，避免搜索完所有信息单独写提取花更多时间；先写反应物再写产物
                tarBOinformcopy3={};%记录反应物键级信息，对应坐标还未提取
                
                frameproact=react_blocklocate(tarBOinformcopy,tarBOinform,tartrajectory);%得到该时间步反应物tarBOinform各键级块区/
                [lenframe,~]=size(frameproact);pronum=0;
                for i=1:lenframe
                    if frameproact(i,5)~=0
                        pronum=pronum+1;
                        tarBOinformcopy(frameproact(i,2):frameproact(i,3),:)=[];%清除找到反应物的产物键级信息
                        for j=i+1:lenframe%更新删除命中行后frameproact中巨鹿的tarBOinformcopy中删除行以下的各目标物质区块行号
                            frameproact(j,2)=frameproact(j,2)-frameproact(j,4);
                            frameproact(j,3)=frameproact(j,3)-frameproact(j,4);
                        end
                        [~,colframe]=size(frameproact);[numzero,~]=ismember(frameproact(i,:),0);%确定反应物数目
                        numreactant=(colframe-5-sum(numzero))/2;
                        if numreactant>0%找到反应物键级信息，存储到tarBOinformcopy3中
                            for j=1:numreactant
                                [row3,~]=size(tarBOinformcopy3);
                                tarBOinformcopy3(row3+1:row3+frameproact(i,2*j+5)-frameproact(i,2*j+4)+1,:)=tarBOinform(frameproact(i,2*j+4):frameproact(i,2*j+5),:);
                            end
                        else
                            if choi==2
                                error('有反应物命中，但是没有找到块区编号，请检查！！！')
                            elseif choi==4
                                error('有产物命中，但是没有找到块区编号，请检查！！！')
                            end
                        end
                    end
                end
                
                [row4,~]=size(tarBOinformcopy3);[row5,~]=size(tarBOinformfullcopy);%找到反应物其他原子在产物所在帧所在的物质，注意可能在这过程中有其他反应发生，导致原子数仍然不一样；必包含目标产物
                blockNO=[];
                for k=1:row4
                    m=1;%对于每一个新的反应物区块，开始新的计数比对
                    for ii=1:row5
                        if strcmp(tarBOinformfullcopy{ii,1},'#')
                            m=m+1;
                        end
                        if ~strcmp(tarBOinformcopy3{k,1},'#') && ~strcmp(tarBOinformfullcopy{ii,1},'#') && tarBOinformcopy3{k,1}==tarBOinformfullcopy{ii,1}%注意规避数值与符号的ASCII码相等，这里'#'==35是为真的！！！
                            if ~ismember(m,blockNO)%避免有大量的块区号重复
                                lenbloc=length(blockNO);blockNO(lenbloc+1)=m;%找到反应物所有原子在产物帧的块区号
                            end
                        end
                    end
                end
                blockNO=sort(blockNO);
                
                
                tarBOinformcopy2={};
                for i=1:length(blockNO)%重新记录产物帧键级信息，包含目标产物和反应物其他原子所在产物
                    [row6,~]=size(tarBOinformcopy2);
                    tarBOinformcopy2(row6+1:row6+productnum(blockNO(i),3),:)=tarBOinformfullcopy(productnum(blockNO(i),1):productnum(blockNO(i),2),:);
                end
                
                %tartrajectory={tartrajectory(1)};%得到目标轨迹元胞，num型转cell型!!!
                rawdatatrj=fopen(datanametrj,'r');%采用逐行读取处理lammpstrj文件
                lammpstrj_analysis%反应物帧键级信息在tarBOinformcopy3中，键级信息在trjdata中；产物键级信息在tarBOinformcopy2中，trjdata在trjdatacopy中
                 %对坐标进行unwrap
                 if strcmpi(unwrapans,'y')
                     load('BondRadii.mat');
                     trjdata=PBC_Unwrap(tarBOinform,trjdata,BOXsize,boxsize,element,BondRadii);
                 end
                fclose(rawdatatrj);%关闭lammpstrj文件
                if choi==2
                    if formatout==1
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','reactant','-base',num2str(base),'.xyz');%文件名为产物名-反应物时间步-产物时间步-产物个数-反应物-文件格式
                    elseif formatout==2
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','reactant','-base',num2str(base),'.car');%文件名为产物名-反应物时间步-产物时间步-产物个数-反应物-文件格式
                    elseif formatout==3
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','reactant','-base',num2str(base),'.pdb');%文件名为产物名-反应物时间步-产物时间步-产物个数-反应物-文件格式
                    end
                elseif choi==4
                    if formatout==1
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','product','-base',num2str(base),'.xyz');%文件名为产物名-反应物时间步-产物时间步-产物个数-反应物-文件格式
                    elseif formatout==2
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','product','-base',num2str(base),'.car');%文件名为产物名-反应物时间步-产物时间步-产物个数-反应物-文件格式
                    elseif formatout==3
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','product','-base',num2str(base),'.pdb');%文件名为产物名-反应物时间步-产物时间步-产物个数-反应物-文件格式
                    end
                end
                tarBOinform=tarBOinformcopy3;%分别写反应物和产物写car文件，因为原子数目变化，arc文件不能正常播放后面的帧，故写car文件
                xyz_car_pdb_filemaker%开始写入反应物
                
                if choi==2
                    if formatout==1
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','product','-base',num2str(base),'.xyz');%文件名为产物名-反应物时间步-产物时间步-产物个数-产物-文件格式
                    elseif formatout==2
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','product','-base',num2str(base),'.car');%文件名为产物名-反应物时间步-产物时间步-产物个数-产物-文件格式
                    elseif formatout==3
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','product','-base',num2str(base),'.pdb');%文件名为产物名-反应物时间步-产物时间步-产物个数-产物-文件格式
                    end
                elseif choi==4
                    if formatout==1
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','reactant','-base',num2str(base),'.xyz');%文件名为产物名-反应物时间步-产物时间步-产物个数-产物-文件格式
                    elseif formatout==2
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','reactant','-base',num2str(base),'.car');%文件名为产物名-反应物时间步-产物时间步-产物个数-产物-文件格式
                    elseif formatout==3
                        datanamecar=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','reactant','-base',num2str(base),'.pdb');%文件名为产物名-反应物时间步-产物时间步-产物个数-产物-文件格式
                    end
                end
                tarBOinform=tarBOinformcopy2;trjdata=trjdatacopy;
                xyz_car_pdb_filemaker%开始写入产物
                
                seekBOinform%写出tarBOinformcopy3中反应物（choi=2）或产物（choi=4）键级信息,用于查找元素坐标，与MS结合揭示机理
                
                pairnum=pairnum+1;
            else
                if choi==2
                    fprintf('\n\n\n第%d组反应物(%d)-产物(%d)匹配组未命中任何反应物，已命中%d组匹配，继续下一帧回溯反应物\n',loopnum,tartrajectory{1,1},tartrajectorycopy(1),pairnum);
                elseif choi==4
                    fprintf('\n\n\n第%d组产物(%d)-反应物(%d)匹配组未命中任何反应物，已命中%d组匹配，继续下一帧搜寻产物\n',loopnum,tartrajectory{1,1},tartrajectorycopy(1),pairnum);
                end
                loopnum=loopnum+1;
            end  
        else
            error('反应物帧超出第一帧，此帧不存在，请检查！！！');
        end
        [tarraw,~]=size(tarBOinformcopy);%控制找完目标产物帧所有产物的反应物帧
        if pairnum==numstop
            tarraw=0;
		elseif pairnum>numstop
			fprintf('警告！！！搜索到的匹配对大于终止限制，请检查文件限制数目输入值');
        end
    end
    if choi==2
        inspectBOinform%输出产物（choi=2）键接信息
        fprintf('\n目标帧%d产物索引的反应物-产物搜索完毕，共%d组反应物-产物car文件对',tartrajectorycopy(1),pairnum);
    elseif choi==4
        inspectBOinform%输出反应物（choi=4）键接信息
        fprintf('\n目标帧%d反应物索引的反应物-产物搜索完毕，共%d组反应物-产物car文件对',tartrajectorycopy(1),pairnum);
    end
    
    msgbox('反应物-产物搜索完毕，成功生成了car文件对');
   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
elseif choi==3
    fprintf('\n请输入搜索的基团、自由基或者化学键资料，请查看文件夹内输入指导：\n');
    fprintf('\n1.化学基团');
    
    
    
    
    
    
    
    
else
    disp('非法输入，未预设的子程序编号，请检查！！！')
end
fprintf('\n\nchemi_mechanism运行结束\n');
fprintf('目标基团、自由基或者化学键随时间步变化的信息储存在tarstepspecies中\n');
msgbox('chemi_mechanism运行结束');

toc %结束计时
fprintf('\n本次运行耗时：%.2f s\n',toc)

clear blockNO bondset bondsetdata BOXsize choi colframe datanamebond datanametrj date elementsequence fidBO fileheader frame frameproact
clear ii irow jj lenbloc lenframe loopnum m molecomp numreactant numstop numzero outcol pairnum PBC PBCa PBCalpha PBCb PBCbeta PBCc
clear PBCchoi PBCgamma productnum prompt promptans promptans2 promptans3 pronum rawdata rawdatatrj reactname residueseqname
clear row3 row4 row5 row6 rowcopy rowfram spacegroupname species tarelenummatch tarraw tartrajectory tartrajectoryact tartrajectorycopy
clear tartrjdata title trajper loop datanamespe i j seekacc sperunans trajectorynote producttname col maxchoi maxdata outputdata_copy
clear atomid_conv base boxsize counter datanamecar eleswap eleswapans modcounter moddivid formatout datarep errorexe fram_num_check
clear rerun_ans sperunans2 outputdatanew_frame rerun_ans2 species_frame_checkans unwrapans species_frame_check2 readline MOLE
clear iijj elemax datalinenum datacellnum control checknum check_control_origin 
