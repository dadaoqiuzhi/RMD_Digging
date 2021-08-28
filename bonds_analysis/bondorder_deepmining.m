%scrit file name bondorder_deepmining
%purpose:
%This program is used to analyze  bond order information of complex molecule
%in a specific trajectory.
%version 1;2018.6.29
disp('Welcome!--by Qiang Liu @Institute of Nuclear Physics and Chemistry, China Academy of Engineering Physics; Email: liubinqiang@163.com');
disp('Repository adress of the Source code on github: https://github.com/dadaoqiuzhi/RMD_Digging');

fprintf('This program is will analyze the BO information in bondoutdata, classifiec by molecular formula')
speciestrjnum=input('\nPlease input the trajectory timestep. It can be obtained from the analysis results of species files: \n');
elementsequence=input('Please input the atom type, e.g.C H O N etc with space interval, should be in line with the in.* or data file, especially for these with elements mapping to different elements.: \n','s');
disp('bondorder_deepmining is running, please wait...')
element=upper(elementsequence);
element=strtrim(element);element=strsplit(element);
eleswapans=input('\nIs there elements mapping to different elements? y/n: \n','s');eleswapans=lower(eleswapans);
if strcmp(eleswapans,'y')
    fprintf('\nPlease input the cell array of mapped elements, e.g.{''Si'',''C'';''S'',''O''},single quotes in practical, meaning Si-->C and S-->O:\n');
    eleswap=input('');
    eleswap=upper(eleswap);
else
    eleswap={'Nan'};
end
numseq={};
for i=1:length(element)
    numseq{1,i}=i;
end
[row,~]=size(bondoutdata);%find the trajectory in bondoutdata
bondrownum=0;
for i=1:row
    if strcmp(bondoutdata{i,1},'Timestep')
        if bondoutdata{i,2}==speciestrjnum
            tartrjnum=i;%Timestep num row
            break
        end
    end
end
tarbondnum=[];
for i=tartrjnum+1:row
    if strcmp(bondoutdata{i,1},'Timestep')
        tarbondnum=i;
        break
    end
end
if ~isempty(tarbondnum)
    tarbondnum=tarbondnum-tartrjnum-1;%obtain the BO information of the target trajectory
else
    tarbondnum=row-tartrjnum;
end

separator={'#','#','#','#','#','#','#','#','#','#','#','#','#','#','#'};%to seperate the BO of different species
bondoutdata(row+1,:)=separator(1,:);
tarelenummatch={};tarBOinform={};lineofelenum=1;
while ~ischar(bondoutdata{tartrjnum+1,1})
    datapython={};datapython{1,1}=bondoutdata{tartrjnum+1,1};
    for i=1:bondoutdata{tartrjnum+1,3}
        datapython{1,i+1}=bondoutdata{tartrjnum+1,3+i};
    end
    elenummatch={};elenummatch(:,1)=element';
    for i=1:length(element)
        elenummatch{i,2}=0;
    end
    [elenumrow,~]=size(elenummatch);
    elementname=charnum_match(element,numseq,bondoutdata{tartrjnum+1,2});
    if ismember(elementname,eleswap(:,1)) 
        [~,lib]=ismember(elementname,eleswap(:,1));
        elementname=eleswap{lib,2};
    end
    elenummatch=eleme_molecule(elenummatch,elementname);
    BOinform={};BOinform(1,:)=bondoutdata(tartrjnum+1,:);
    lineofbo=2;
    bondoutdata{tartrjnum+1,1}='NaN';
    datapython{1,1}='NaN';
    bondoutdata=cellrowcol_del(bondoutdata,'delrow','NaN');
    tarbondnum=tarbondnum-1;
    datapython=cellrowcol_del(datapython,'delcol','NaN');
    
    
    while ~isempty(datapython) && ~ischar(bondoutdata{tartrjnum+1,1})
        alter=datapython{1,1};kk=0;
        for i=tartrjnum+1:tartrjnum+tarbondnum
            if alter==datapython{1,1}
                kk=kk+1;
                if datapython{1,1}==bondoutdata{i,1}
                    j=length(datapython)+1;
                    for k=1:bondoutdata{i,3}
                        datapython{1,j}=bondoutdata{i,3+k};
                        j=j+1;
                    end
                    
                    datapython{1,1}='NaN';
                    datapython=cellrowcol_del(datapython,'delcol','NaN');
                    elementname=charnum_match(element,numseq,bondoutdata{i,2});
                    if ismember(elementname,eleswap(:,1)) 
                        [~,lib]=ismember(elementname,eleswap(:,1));
                        elementname=eleswap{lib,2};
                    end
                    elenummatch=eleme_molecule(elenummatch,elementname);
                    BOinform(lineofbo,:)=bondoutdata(i,:);
                    lineofbo=lineofbo+1;
                    bondoutdata{i,1}='NaN';
                    bondoutdata=cellrowcol_del(bondoutdata,'delrow','NaN');
                    tarbondnum=tarbondnum-1;
                    break
                end
                if alter==datapython{1,1} && kk==tarbondnum 
                    datapython{1,1}='NaN';
                    datapython=cellrowcol_del(datapython,'delcol','NaN');
                    break
                end
            end
        end
    end
    tarelenummatch(1:elenumrow,2*lineofelenum-1:2*lineofelenum)=elenummatch(1:elenumrow,1:2);
    lineofelenum=lineofelenum+1;
    [BOrow,~]=size(BOinform);
    BOinform(BOrow+1,:)=separator(1,:);
    [rowtarBO,~]=size(tarBOinform);
    tarBOinform(rowtarBO+1:rowtarBO+BOrow+1,:)=BOinform(1:BOrow+1,:);
    continue
end

fprintf('\nbondorder_deepmining is successfully finished\n')
disp('Molecular formula is saved in tarelenummatch with corresponding BO information in tarBOinform')

clear alter bondrownum BOrow col datapython element elementname elementsequence elementsequence
clear elenummatch elenumrow i j k kk  lineofbo lineofelenum numseq row rowtarBO separator speciestrjnum
clear tarbondnum tartrjnum trajectorynum

