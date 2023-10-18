%scrit file name orient_entangle
%purpose:
%This program is used to analyze orientation parameter & entanglement degree
%in a specific trajectory.
%include bonds_analysis_speedup,bondorder_deepmining,lammpstrj_analysis
%reference
%Molecular dynamics simulations of deformation mechanisms of amorphous polyethylene
%半晶态聚合物拉伸变形的微观机理;碳纳米管和石墨烯诱导聚乙烯链取向的分子动力学模拟
%Molecular dynamics simulation of deformation behavior in amorphous polymer: nucleation of chain entanglements and
%network structure under uniaxial tension
%version 1;2019.3.27

disp('##################################################################################################################################')
disp('Welcome!--by Qiang Liu @Institute of Nuclear Physics and Chemistry, China Academy of Engineering Physics; Email: liubinqiang@163.com');
disp('Repository adress of the Source code on github: https://github.com/dadaoqiuzhi/RMD_Digging');
disp('References: 1.Fuel 287 (2021) 119484. 2.ACS Appl. Mat. Interfaces 13(34) (2021) 41287-41302. 3.ACS Appl. Mat. Interfaces 2022, 14.(4), 5959-5972.')
disp('4.ACS Materials Letters 2023, 2174-2188. More work is coming!')
disp('##################################################################################################################################')

fprintf('\nThis program is used to calxulated the orientation parameter and entanglement degree of polymer, H, CH3, terminal carbonyl/phenol O are not considered,\n the para- and meta- positions of benze ring are considered via its geometric center, notice: severe degradation should be avoided\n');
dataname=input('Please input the file name of bond.* file：\n','s');
trajper=input('Please input the output frequency of BO information (Positive integer):\n');


tartrajectoryset=input('\nPlease input the trajectory number of the BO information, multiple trajectories should be seperated by white space: \n','s');
tartrajectoryset=strtrim(tartrajectoryset);tartrajectoryset=strsplit(tartrajectoryset);
for i=1:length(tartrajectoryset)
    tartrajectoryset{1,i}=str2num(tartrajectoryset{1,i});
end
tartrajectoryset=cell2mat(tartrajectoryset);


atomnum=input('Please input atom number: \n');
atomnumcopy=atomnum;
elementsequence=input('Please input atom type like C,H.O,N, seperated by white space, corresponding to 1,2,3,4...n (see *.data or in.*, \nespecially for element mapping: \n','s');
datanametrj=input('\nFile name of *.lammpstrj file：\n','s');
eunitvect=input('\nPlease input the direction matrix of stress or orientation, eg. [0,0,1]：\n');
entanglength=input('\nPlease input the interval number between entangled monomer, (according to model or critical entangled length)：\n');
if entanglength<=0 || entanglength~=round(entanglength)
    error('\ninterval number between entangled monomer, please check it!');
end


fprintf('\nCriterion for molecular chain: \n1.limited by C atom number，2.limited by O atom number，3.limited by N atom number，4.limited by O atom number，5.limited by O atom number，6.limited by H atom number，7.no limitation');
anaoption=input('\nPlease select the No. of molecular chain, must remove small molecules, multiple No. should be seperated by comma in matrix ([1,2]), No.7 have mutual exclusion with other option: \n');
if ~ismember(anaoption,[1,2,3,4,5,6,7])
    error('No. for molecular chain is wrong, please check it!');
end
if ismember(1,anaoption)
    Cnumrequest=input('\nPlease input limiting condition for C atom, logical operator should be included,eg. (< 20,> 10,< = 25,= 18,str should be seperated by white space): \n','s');
    Cnumrequest=strtrim(Cnumrequest);Cnumrequest=strsplit(Cnumrequest);Cnumrequest{length(Cnumrequest)}=str2num(Cnumrequest{length(Cnumrequest)});
    if length(Cnumrequest)<2
        error('Illegal input,please input limiting condition!');
    end
else
    Cnumrequest={0};
end
if ismember(2,anaoption)
    Onumrequest=input('\nPlease input limiting condition for O atom, logical operator should be included,eg. (< 20,> 10,< = 25,= 18,str should be seperated by white space): \n','s');
    Onumrequest=strtrim(Onumrequest);Onumrequest=strsplit(Onumrequest);Onumrequest{length(Onumrequest)}=str2num(Onumrequest{length(Onumrequest)});
    if length(Onumrequest)<2
        error('Illegal input,please input limiting condition!');
    end
else
    Onumrequest={0};
end
if ismember(3,anaoption)
    Nnumrequest=input('\nPlease input limiting condition for N atom, logical operator should be included,eg. (< 20,> 10,< = 25,= 18,str should be seperated by white space): \n','s');
    Nnumrequest=strtrim(Nnumrequest);Nnumrequest=strsplit(Nnumrequest);Nnumrequest{length(Nnumrequest)}=str2num(Nnumrequest{length(Nnumrequest)});
    if length(Nnumrequest)<2
        error('Illegal input,please input limiting condition!');
    end
else
    Nnumrequest={0};
end
if ismember(4,anaoption)
    Hnumrequest=input('\nPlease input limiting condition for H atom, logical operator should be included,eg. (< 20,> 10,< = 25,= 18,str should be seperated by white space): \n','s');
    Hnumrequest=strtrim(Hnumrequest);Hnumrequest=strsplit(Hnumrequest);Hnumrequest{length(Hnumrequest)}=str2num(Hnumrequest{length(Hnumrequest)});
    if length(Hnumrequest)<2
        error('Illegal input,please input limiting condition!');
    end
else
    Hnumrequest={0};
end
if ismember(5,anaoption)
    Totalnumrequest=input('\nPlease input limiting condition for tital atom, logical operator should be included,eg. (< 20,> 10,< = 25,= 18,str should be seperated by white space): \n','s');
    Totalnumrequest=strtrim(Totalnumrequest);Totalnumrequest=strsplit(Totalnumrequest);Totalnumrequest{length(Totalnumrequest)}=str2num(Totalnumrequest{length(Totalnumrequest)});
    if length(Totalnumrequest)<2
        error('Please input limiting condition for O atom, logical operator should be included,eg. (< 20,> 10,< = 25,= 18,str should be seperated by white space): \n','s'');
    end
else
    Totalnumrequest={0};
end
if ismember(6,anaoption)
    MWrequest=input('\nPlease input limiting condition for molecular weight, logical operator should be included,eg. (< 2000,> 100,< = 250000,= 1862,str should be seperated by white space): \n','s');
    MWrequest=strtrim(MWrequest);MWrequest=strsplit(MWrequest);MWrequest{length(MWrequest)}=str2num(MWrequest{length(MWrequest)});
    if length(MWrequest)<2
        error('Illegal input,please input limiting condition!');
    end
else
    MWrequest={0};
end
if ismember(7,anaoption) && length(anaoption)>1
    error('\nOption 7 is not compatible with other option, please check it!');
end



orient_entangledata={};
for loop=1:length(tartrajectoryset)
    tartrajectory=tartrajectoryset(loop);
    len=size(orient_entangledata,1);
    orient_entangledata{len+1,1}=tartrajectory;orient_entangledata{len+1,2}='#';orient_entangledata{len+1,3}='#';
    orient_entangledata{len+2,1}='Mole';orient_entangledata{len+2,2}='orientation degree';orient_entangledata{len+2,3}='Entanglement';
    
    bonds_analysis_speedup
    bondorder_deepmining
    lammpstrj_analysis
    
    
    Porientdata=[];
    if ~ismember(7,anaoption)
        [~,col]=size(tarelenummatch);BOrderid=[];
        for ii=1:col/2
            logicvalue=anaoptlogic(anaoption,ii,tarelenummatch,Cnumrequest,Onumrequest,Nnumrequest,Hnumrequest,Totalnumrequest,MWrequest);
            if logicvalue==1
                tarBOinformcopy={};
                [row,~]=size(tarBOinform);line=0;
                while row
                    line=line+1;
                    row=row-1;
                    if strcmp(tarBOinform{line,1},'#')
                        tarBOinformcopy(1:line,1:7)=tarBOinform(1:line,1:7);
                        tarBOinformcopy(line,:)=[];
                        tarBOinform(1:line,:)=[];
                        line=0;
                        break;
                    end
                end
                endBO=chainend(tarBOinformcopy,element,tarelenummatch,ii,tartrajectory);
                [backboneBO,delatomid]=delnonbackbone(tarBOinformcopy,element);
                BOrderid=order_idtag(endBO,backboneBO,delatomid,element);
                xyz_Posi=Posi_map(BOrderid,trjdata);
                Porient=orient_cal(xyz_Posi,eunitvect);
                Porientdata(length(Porientdata)+1,1)=Porient;
                [~,degnum]=entanglement_cal(xyz_Posi,entanglength);
                
                moleform='';
                for i=1:4
                    if tarelenummatch{i,2*ii}~=0
                        moleform=strcat(moleform,tarelenummatch{i,2*ii-1},num2str(tarelenummatch{i,2*ii}));
                    end
                end
                
                len=size(orient_entangledata,1);
                orient_entangledata{len+1,1}=moleform;
                orient_entangledata{len+1,2}=Porient;
                orient_entangledata{len+1,3}=degnum;
                
                
                tarBOinformcopy={};
                
            else%
                moleform='';
                for i=1:4
                    if tarelenummatch{i,2*ii}~=0
                        moleform=strcat(moleform,tarelenummatch{i,2*ii-1},num2str(tarelenummatch{i,2*ii}));
                    end
                end
                fprintf('\nThe %d group molecular formula %s dissatisfy the orientyation condition, removed',ii,moleform);
                [row,~]=size(tarBOinform);line=0;
                while row
                    line=line+1;
                    row=row-1;
                    if strcmp(tarBOinform{line,1},'#')
                        tarBOinform(1:line,:)=[];
                        line=0;
                        break;
                    end
                end
            end
        end
        
        
    else
        [~,col]=size(tarelenummatch);
        for ii=1:col/2
            [~,col]=size(tarelenummatch);BOrderid=[];
            tarBOinformcopy={};
            [row,~]=size(tarBOinform);line=0;
            while row
                line=line+1;
                if strcmp(tarBOinform{line,1},'#')
                    tarBOinformcopy(1:line,1:7)=tarBOinform(1:line,1:7);
                    tarBOinformcopy(line,:)=[];
                    tarBOinform(1:line,:)=[];
                    row=row-line;
                    line=0;
                    break;
                end
            end
            endBO=chainend(tarBOinformcopy,element,tarelenummatch,ii,tartrajectory);
            [backboneBO,delatomid]=delnonbackbone(tarBOinformcopy,element);
            BOrderid=order_idtag(endBO,backboneBO,delatomid,element);
            xyz_Posi=Posi_map(BOrderid,trjdata);
            Porient=orient_cal(xyz_Posi,eunitvect);
            Porientdata(length(Porientdata)+1,1)=Porient;
            [entangle_deg,degnum]=entanglement_cal(xyz_Posi,entanglength);
            len=size(orient_entangledata,1);
            orient_entangledata{len+1,1}=moleform;
            orient_entangledata{len+1,2}=Porient;
            orient_entangledata{len+1,3}=degnum;
            
            tarBOinformcopy={};
            
        end
    end

    fprintf('\nThe %d th trajectory %d has been disposed\n',loop,tartrajectoryset(loop));
    if loop<length(tartrajectoryset)
        fprintf('\nFurther deal with the %d th trajectory %d\n\n',loop+1,tartrajectoryset(loop+1));
    end
end



fprintf('\nThe orientation parameter and entanglement degree of respective molecular formula is saved in orient_entangledata, seperated by # character\n');

clear ans atomnum bondnumdata control datacell datacellchar datadel dataline dataname datarep datasplit found gap i j k kk line 
clear outputans rawdata tartrajectory trajper unfound dataoutrow dataoutcol dataoutputrow dataoutcolchar dataoutputcol filename
clear alter bondrownum BOrow col datapython element elementname elementsequence elementsequence 
clear elenummatch elenumrow i j k kk  lineofbo lineofelenum numseq row rowtarBO separator speciestrjnum
clear tarbondnum tartrjnum trajectorynum readline
clear atomnum control datacell datacellchar dataline dataname datarep datasplit found gap i line rawdata tartrajectory trajper unfound ans
clear outputans dataoutrow dataoutcol dataoutputrow dataoutcolchar dataoutputcol filename atomnumcopy
clear anaoption BOinform Cnumrequest datanametrj degnum entanglength eunitvect Hnumrequest ii len logicvalue loop moleform
clear MWrequest Nnumrequest Onumrequest Porient Porientdata rawdatabond rawdatatrj tartrajectoryset Totalnumrequest


