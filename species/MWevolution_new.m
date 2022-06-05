%%scrit file name MWevolution
%purpose:
%calculate the evolution of number-average moleculae weight (Mn), weight-average molacular weight(Mw) and molecular weight distribution (MWD)
%include species_analysis.m
%version 1;2018.10.5
disp('##################################################################################################################################')
disp('Welcome!--by Qiang Liu @Institute of Nuclear Physics and Chemistry, China Academy of Engineering Physics; Email: liubinqiang@163.com');
disp('Repository adress of the Source code on github: https://github.com/dadaoqiuzhi/RMD_Digging');
disp('References: 1.Fuel 287 (2021) 119484. 2.ACS Appl. Mat. Interfaces 13(34) (2021) 41287-41302. More work is coming!')
disp('##################################################################################################################################')
disp('Plwase run species_analysis first, then run MWevolution to obtain Mw and Mn data')
MDans=input('\nPlease select the option No.: 1. only MW   2.only molecular weight distribution (MWD)  3.both 1 and 2: \n');

if MDans==1 || MDans==3
    massshresholdans=input('\n Limit the Mw threshold value of species used to calculate? y/n: \n','s');
    if ~ismember(massshresholdans,{'y','n'})
        error('illegal argument for Mw threshold value, Please check it!!!');
    end
    if strcmpi(massshresholdans,'y')
        fprintf('\nPlease select the option No.: \n1.larger than a Mw threshold value\n2.less than a Mw threshold value\n3.between two Mw threshold values(closed set): \n');
        massshresholdans2=input('\nPlease select the option No.: \n');
        if massshresholdans2==1 || massshresholdans2==2
            massshreshold=input('\nPlease input the Mw threshold value: \n');
        elseif massshresholdans2==3
            massshreshold=input('\nPlease input the Mw threshold value matrix, e.g.[200 50000]: \n');
            massshreshold=sort(massshreshold,2);
        end
    end
end

if MDans==2 || MDans==3
    MWDans=input('\nIf export the MWD data? y/n: \n','s');
    if strcmpi(MWDans,'y')
        MWDfram=input('\nPlease select the option No.: 1.Specify frames manually  2.monotonically increasing frames (closed set) : \n');
        if ~ismember(MWDfram,[1 2])
            error('Illegal option No., please check it!!!');
        end
        STEP=[];
        switch MWDfram
            case 1
                MWDfram=input('\nPlease input the interested frames with white space for multiple parameters: \n','s');
                MWDfram=strtrim(MWDfram);MWDfram=strsplit(MWDfram);
                for i=1:length(MWDfram)
                    STEP(i)=str2num(MWDfram{1,i});
                end
            case 2
                minstep=input('\nPlease input the minimum timestep: \n');
                maxstep=input('\nPlease input the maximum timestep: \n');
                dumpevery=input('\nPlease input the output frequency of species file by lammps: \n');
                STEP(1)=maxstep(1);STEP(2)=minstep(1);
                MWDfram=input('\nPlease input the common difference of monotonically increasing frames, must be the integral multiple of output frequency of species file: \n');
                if mod(MWDfram,dumpevery)~=0
                    error('Illeagal common difference, please check it!!!');
                end
                while STEP(length(STEP))<maxstep(1)
                    STEP(length(STEP)+1)=STEP(length(STEP))+MWDfram;
                    if STEP(length(STEP))>maxstep(1)
                        break;
                    end
                end
        end
    end
    STEP=sort(STEP);
end

disp('MWevolution_new is running, please wait...')
datadelimiter={'C','H','O','N','He','Li','Be','B','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Pd','Ag','Cd','In','Sn','Sb','I','Xe','Cs','Ba','Pt','Au','Hg','Pb'};
matchdataMD=[];[row,col]=size(outputdata);
for i=4:col%
    [C,matches]=strsplit(outputdata{1,i},datadelimiter,'CollapseDelimiters',false);
    classmatch={};
    C=delnull(C);
    if length(matches)~=length(C)
        [C,matches]=strsplit(outputdata{1,i},datadelimiter,'CollapseDelimiters',false);
        C=C(1,2:end);
        for k=1:length(C)
            if isempty(C{k})
                C{k}='1';
            end
        end
    end
    for j=1:length(matches)
        classmatch{j,1}=matches{j};
        classmatch{j,2}=C{j};
    end
    matchdataMD(i-3)=molecuweight(classmatch);
end

if MDans==1 || MDans==3
    MD=[];Mn=0;Mw=0;molenum={};
    Mn_num=0;Mw_num=0;
    for i=2:row
        Mn=0;Mw=0;
        molenum={};
        Mn_num=0;
        Mw_num=0;
        
        molenum=outputdata(i,4:end);
        molenum=cell2mat(molenum);
        for k=1:col-3
            molenum(2,k)=molenum(1,k)*matchdataMD(k);
        end
        
        
        if strcmpi(massshresholdans,'n')
            for j=1:col-3
                Mn=Mn+molenum(1,j)*matchdataMD(j)/sum(molenum(1,:));
                Mw=Mw+molenum(2,j)*matchdataMD(j)/sum(molenum(2,:));
            end
        elseif strcmpi(massshresholdans,'y')
            if massshresholdans2==1
                for j=1:col-3
                    if matchdataMD(j)>=massshreshold
                        Mn=Mn+molenum(1,j)*matchdataMD(j);
                        Mn_num=Mn_num+molenum(1,j);
                        Mw=Mw+molenum(2,j)*matchdataMD(j);
                        Mw_num=Mw_num+molenum(2,j);
                    end
                end
            elseif massshresholdans2==2
                for j=1:col-3
                    if matchdataMD(j)<=massshreshold
                        Mn=Mn+molenum(1,j)*matchdataMD(j);
                        Mn_num=Mn_num+molenum(1,j);
                        Mw=Mw+molenum(2,j)*matchdataMD(j);
                        Mw_num=Mw_num+molenum(2,j);
                    end
                end
            elseif massshresholdans2==3
                for j=1:col-3
                    if matchdataMD(j)>=massshreshold(1,1) && matchdataMD(j)<=massshreshold(1,2)
                        Mn=Mn+molenum(1,j)*matchdataMD(j);
                        Mn_num=Mn_num+molenum(1,j);
                        Mw=Mw+molenum(2,j)*matchdataMD(j);
                        Mw_num=Mw_num+molenum(2,j);
                    end
                end
            end
        end
        
        
        if strcmpi(massshresholdans,'n')
            MD(i-1,1)=Mn;
            MD(i-1,2)=Mw;
        elseif strcmpi(massshresholdans,'y')
            if Mn_num==0
                MD(i-1,1)=0;
            else
                MD(i-1,1)=Mn/Mn_num;
            end
            if Mw_num==0
                MD(i-1,2)=0;
            else
                MD(i-1,2)=Mw/Mw_num;
            end
        end
        
    end
end

if MDans==2 || MDans==3
    k=2;MWDdata=[];i=2;
    MWDdata(1,:)=matchdataMD(1,:);
    for j=1:length(STEP)
        while k <= size(outputdata,1)
            if outputdata{k,1}==STEP(j)
                for m = 4:size(outputdata,2)
                    MWDdata(i,m-3)=outputdata{k,m};
                end
                i=i+1;
                break;
            else
                k=k+1;
                if k > size(outputdata,1)
                    break;
                end
            end
        end
    end
    MWDdatacopy=MWDdata';
    MWDdatacopy=sortrows(MWDdatacopy,1);
    MWDdata=MWDdatacopy';

    i=1;j=size(MWDdata,2);k=1;
    while k
        if i<j
            if MWDdata(1,i)==MWDdata(1,i+1)
                MWDdata(:,i)=MWDdata(:,i)+MWDdata(:,i+1);
                MWDdata(:,i+1)=[];
                j=j-1;
            else
                i=i+1;
            end
        else
            break;
        end
    end
end
MWDdata2=[];

%MD(isnan(MD)==1)=0;
if MDans==1 || MDans==3
    fprintf('\nMWevolution is successfully finished. Mn and Mw is saved in MD matrix\n')
end
if MDans==2 || MDans==3
    fprintf('\nMWD is saved in MWDdata, the first row is molecular weight, the following data is species number. The 0 value can be further treated by delzero procedure\n');
end
if MDans==1
    msgbox('Molecular weight calculation is finished!');
elseif MDans==2
    msgbox('MWD calculation is finished!');
elseif MDans==3
    msgbox('Molecular weight and MWD is successfully finished!');
end


clear C col i j k matches Mn Mw row MWDans dumpevery maxstep minstep MWDfram MWDans MWDdatacopy Mw_num Mn_num massshresholdans
clear m datadelimiter classmatch MDans molenum
