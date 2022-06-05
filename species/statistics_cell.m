%scrit file name statistics
%purpose:
%This is the main program used to analysis species file
%version 1;2018.6.24
disp('##################################################################################################################################')
disp('Welcome!--by Qiang Liu @Institute of Nuclear Physics and Chemistry, China Academy of Engineering Physics; Email: liubinqiang@163.com');
disp('Repository adress of the Source code on github: https://github.com/dadaoqiuzhi/RMD_Digging');
disp('References: 1.Fuel 287 (2021) 119484. 2.ACS Appl. Mat. Interfaces 13(34) (2021) 41287-41302. More work is coming!')
disp('##################################################################################################################################')
fprintf('\nCell data is required\n');
statans=input('Average data? y/n:\n','s');
statans=lower(statans);
if statans=='y'
    datastat=input('data name in work space, in general outputdata, outputdatanew or dataexport: \n');
    timestep=input('Please input timestep, unit:fs: \n');
    avesteps=input('Please input frame span for data average, must be an integer: \n');
    thermoper=datastat{3,1}-datastat{2,1};
    thermomax=datastat{size(datastat,1),1};
    fprintf('\n\nstatistics is running, please wait...\n\n')
    imax=thermomax/thermoper+1;re=mod(imax,avesteps);
    imaxend=(imax-re)/avesteps;
    [~,datacol]=size(datastat);outdatastat={};
    outdatastat(1,:)=datastat(1,:);
    if strcmp(datastat{1,1},'Timestep')
        for i=2:imaxend
            outdatastat{i,1}=datastat{(i-2)*avesteps+2,1};
        end
        
        for i=2:imaxend
            for j=2:datacol
                sumnum=0;
                for k=avesteps*(i-2)+2:avesteps*(i-1)+1
                    sumnum= sumnum+datastat{k,j};
                end
                outdatastat{i,j}=sumnum/avesteps;
            end
        end
    else
        for i=2:imaxend
            for j =1:datacol
                sumnum=0;
                for k=avesteps*(i-2)+2:avesteps*(i-1)+1
                    sumnum= sumnum+datastat{k,j};
                end
                outdatastat{i,j}=sumnum/avesteps;
            end
        end
    end
end
for i=2:size(outdatastat,1)
    outdatastat{i,1}=outdatastat{i,1}*timestep/1000;
end
disp('Results of statistics is saved in outdatastat')


outputans=input('Export results to Excel? y/n: n','s');
outputans=lower(outputans);
if outputans=='y'
    [dataoutrow,dataoutcol]=size(outdatastat);
    dataoutputrow=strcat('A','1');
    if dataoutcol<=26
        dataoutcolchar=char(65+dataoutcol-1);
        dataoutputcol=strcat(dataoutcolchar,num2str(dataoutrow));
    else
          dataoutcol=dec2base(dataoutcol,26);
          dataoutcol=num2str(dataoutcol);
          datadelimiter={'1','2','3','4','5','6','7','8','9','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','0'};
          [C,matches]=strsplit(dataoutcol,datadelimiter,'CollapseDelimiters',false);
          charcor=char26cor(matches);
          dataoutputcol=strcat(charcor,num2str(dataoutrow));
    end    
    filename='output_mydata.xlsx';
    xlswrite(filename,outdatastat,dataoutputrow:dataoutputcol)
    disp('Averaged results are saved in outdatastat and exported to excel:output_mydata')
end
fprintf('\n\nstatistics is successfully finished\n\n')
clear avesteps datacol datarow datastat i imax imaxend j k outputans re statans sumnum thermomax thermoper timestep statans
