%scrit file name  statisave 
%purpose:
%This program is used to  statistically average a given file
%version 1;2018.7.15
disp('##################################################################################################################################')
disp('Welcome!--by Qiang Liu @Institute of Nuclear Physics and Chemistry, China Academy of Engineering Physics; Email: liubinqiang@163.com');
disp('Repository adress of the Source code on github: https://github.com/dadaoqiuzhi/RMD_Digging');
disp('References: 1.Fuel 287 (2021) 119484. 2.ACS Appl. Mat. Interfaces 13(34) (2021) 41287-41302. More work is coming!')
disp('##################################################################################################################################')
fprintf('\nMatrix data is required\n');
datasource=input('Where is the source data from: 1.Excel, 2.matrix in workspace. Input the No.: \n');
if datasource==1
   datasource=xlsread('input_data.xlsx');
elseif datasource==2
    datasource=input('Input the name of the matrix in the workspace: \n');
else
    disp('Illegal input! Please check it.')
end
avesteps=input('Step/Span number used to average the data, should be a positive integer: \n');
disp('statisave is running...')
if iscell(datasource)
    if ischar(datasource{1,1})
        datasource(1,:)=[];
        datasource=cell2mat(datasource);
    end
end
[datarow,datacol]=size(datasource);
imax=(datarow-mod(datarow,avesteps))/avesteps;
datastatis=[];
for i=1:imax
    for j=1:datacol
        datastatis{i,j}=mean(datasource(avesteps*(i-1)+1:avesteps*i,j));
    end
end
disp('statisave is successfully accomplished')
fprintf('\nResult is saved in datastatis matrix\n')

clear avesteps datacol datarow i j imax
