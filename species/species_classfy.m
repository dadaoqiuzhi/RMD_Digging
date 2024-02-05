%scrit file name species_classfy
%purpose:
%This program is used to analysis species file
%(1)C20 means species with 20 C, C42+ denotes species with C number larger
%than 42, C100- is species with C less than 100
%(2)M100 indicates species with Mw of 100, M125+ denotes Mw larger than 125, M400- is species less than 400
%(3)eleC are species have C, eleCO are species have C and O
%(4)eleonlyCH are species only have C and H, eleonlyCO are species only have C and O
%version 1;2018.6.23
disp('##################################################################################################################################')
disp('Welcome!--by Qiang Liu @Institute of Nuclear Physics and Chemistry, China Academy of Engineering Physics; Email: liubinqiang@163.com');
disp('Repository adress of the Source code on github: https://github.com/dadaoqiuzhi/RMD_Digging');
disp('References: 1.Fuel 287 (2021) 119484. 2.ACS Appl. Mat. Interfaces 13(34) (2021) 41287-41302. 3.ACS Appl. Mat. Interfaces 2022, 14.(4), 5959-5972.')
disp('4.ACS Materials Letters 2023, 2174-2188. More work is coming!')
disp('##################################################################################################################################')

disp('When species_analysis is executed, this procedure can obtain the interested species')

fprintf('\n(1)C20 means species with 20 C,C42+ denotes species with C number larger than 42, C100- is species with C less than 100')
fprintf('\n(2)M100 indicates species with Mw of 100, M125+ denotes Mw larger than 125, M400- is species less than 400')
fprintf('\n(3)eleC are species have C, eleCO are species have C and O')
fprintf('\n(4)eleonlyCH are species only have C and H, eleonlyCO are species only have C and O')
fprintf('\n\nMethods to filter out the interested species: \na:C1,C20,C42+,C100-,+ means >=,- means <\n')
fprintf('b:M100,M125+,M400-\n')
fprintf('c:eleC,eleCO\n')
fprintf('d:eleonlyC,eleonlyCO\n\n\n')
tarclass=input('Please select the option (a, b, c or d): \n','s');
tarclass=lower(tarclass);
sumans=input('Sum the data? y/n:\n','s');
sumans=lower(sumans);
%long characters in datadelimiter should be list first to avoid find such case: (1) target Cl but find/match C, (2) target Na but find/match N
datadelimiter={'eleonly','ele','Li','Be','He','Ne','Na','Mg','Cl','Ar','Ca','Sc','Ti','Al','Si','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Pd','Ag','Cd','In','Sn','Sb','Xe','Cs','Ba','Pt','Au','Hg','Pb','M','C','H','O','N','+','-','B','F','P','S','K','V','I'};
outputdatast=outputdata(1,:);
fprintf('\nif the exported data only have three column, not hit the interested species, please delet the irrelevant data in the work space\n')
if tarclass=='a'||tarclass=='b'
    classid=input('Please input the specific requirements according to the selected methods a or b, e.g. C100-, M100: \n','s');
    classid=upper(classid);
    clear dataexport matchdatacol sumdata
    [C,matches]=strsplit(classid,datadelimiter,'CollapseDelimiters',false);
    if isempty(matches)
        error('Not match the option, please check input or complete the parameter in datadelimiter')
    end
    Clength=length(C);
    if Clength==2
        classidcell={};
        classidcell{1}=matches{1};
        classidcell{2}=C{2};
    else
        classidcell={};
        classidcell{1}=matches{1};
        classidcell{2}=C{2};
        classidcell{3}=matches{2};
    end
    fprintf('\nspecies_classfy is running, please wait...\n')
end


if tarclass=='c' || tarclass=='d'
    classid=input('Please input the specific requirements according to the selected methods c or d, e.g. eleonlyCH, eleCO: \n','s');
    [C,matches]=strsplit(classid,datadelimiter,'CollapseDelimiters',false);
    classidcell={};
    for i=1:length(matches)
        classidcell{i}=matches{i};
    end
    fprintf('\nspecies_classfy is running, please wait...\n')
end

matchdatacol=[]; kk=1;[row,col]=size(outputdata);
for i=4:col
    [C,matches]=strsplit(outputdatast{i},datadelimiter,'CollapseDelimiters',false);
    classmatch={};
    C=delnull(C);
    if length(matches)~=length(C)
        [C,matches]=strsplit(outputdatast{i},datadelimiter,'CollapseDelimiters',false);
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
    
    
    if tarclass=='a'
        memcheck=ismember(classmatch(:,1),classidcell{1});
        if sum(memcheck)>=1
            [row,~]=size(memcheck);matchrow=[];matchnum=[];
            for j=1:row
                if memcheck(j,1)==1
                    matchrow(length(matchrow)+1,1)=j;
                    matchnum(length(matchnum)+1,1)=str2num(classmatch{j,2});
                end
            end
            if length(classidcell)==2 && sum(matchnum)==str2num(classidcell{2})
                matchdatacol(kk)=i;
                kk=kk+1;
            elseif length(classidcell)==3
                memcheck=ismember(classidcell,{'+'});
                if sum(memcheck)==1
                    if sum(matchnum)>=str2num(classidcell{2})
                        matchdatacol(kk)=i;
                        kk=kk+1;
                    end
                end
                memcheck=ismember(classidcell,{'-'});
                if sum(memcheck)==1
                    if sum(matchnum)<str2num(classidcell{2})
                        matchdatacol(kk)=i;
                        kk=kk+1;
                    end
                end
            end
        end
    end
    
    
    
    if tarclass=='b'
        mw=str2num(classidcell{2});
        datamw=molecuweight(classmatch);
        if length(classidcell)==2
            if mw==round(datamw)
                matchdatacol(kk)=i;
                kk=kk+1;
            end
        end
        if length(classidcell)==3
            memcheck=ismember(classidcell,{'+'});
            if sum(memcheck)==1
                if datamw>=mw
                    matchdatacol(kk)=i;
                    kk=kk+1;
                end
            end
            memcheck=ismember(classidcell,{'-'});
            if sum(memcheck)==1
                if datamw<mw
                    matchdatacol(kk)=i;
                    kk=kk+1;
                end
            end
        end
    end
    
    if tarclass=='c'
        sumcheck=0;
        for j=2:length(classidcell)
            memcheck=ismember(classmatch,classidcell{j});
            sumcheck=sumcheck+sum(sum(memcheck));
        end
        if sumcheck==length(classidcell)-1
            matchdatacol(kk)=i;
            kk=kk+1;
        end
    end
    
    if tarclass=='d'
        if length(classidcell)-1==length(matches)
            sumcheck=0;
            for j=2:length(classidcell)
                memcheck=ismember(classmatch,classidcell{j});
                sumcheck=sumcheck+sum(sum(memcheck));
                if sumcheck==length(matches)
                    matchdatacol(kk)=i;
                    kk=kk+1;
                end
            end
        end
    end  
end



dataexport={};
for j=1:length(matchdatacol)
    dataexport(:,j+3)=outputdata(:,matchdatacol(j));
end
dataexport(:,1:3)=outputdata(:,1:3);
disp('Results of nspecies_classfy are saved in dataexport')
 
if isempty(matchdatacol)
    fprintf('\n\nDo not find the species information, please check it!\n\n')
end

if sumans=='y'
    fprintf('\nSum calculation is running...\n')
    [a,b]=size(dataexport);
    sumdata=[];
    for i=2:a
        sumsum=0;
        for j=4:b
            sumsum=sumsum+dataexport{i,j};
        end
        sumdata(i-1,1)=sumsum;
    end
end
fprintf('\nspecies_classfy is successfully finished\n')
disp('Sum of data are saved in sumdata')



fprintf('\nResults of species_classfy are saved in dataexport, sum of data are saved in sumdata')
fprintf('\nMore complex data abstraction can be performed by copying data in dataexport to outputdata and go on!\n')
msgbox('species_classfy is successfully finished!');

% clear C classid classidcell classmatch Clength col datadelimiter dataoutcol dataoutcolchar dataoutputcol dataoutputrow sumsum
% clear dataoutrow filename i j k kk matches matchnum matchrow memcheck outputdatast row sumcheck tarclass saveans sumans a b mw
