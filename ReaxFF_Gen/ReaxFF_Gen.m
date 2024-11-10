%scrit file name ReaxFF_Gen
%purpose:
%Generation of standard file with ReaxFF force field parameters from literatures
disp('##################################################################################################################################')
disp('Welcome!--by Qiang Liu @Institute of Nuclear Physics and Chemistry, China Academy of Engineering Physics; Email: liubinqiang@163.com');
disp('Repository adress of the Source code on github: https://github.com/dadaoqiuzhi/RMD_Digging');
disp('References: 1.Fuel 287 (2021) 119484. 2.ACS Appl. Mat. Interfaces 13(34) (2021) 41287-41302. More work is coming!')
disp('##################################################################################################################################')

fprintf('Paste data from literatures into the input.txt, one blank Line for each block without the initio lines of pure character string\n');
fidin=fopen('input.txt','r');
fidout=fopen('output.txt','w');
charset=input('\ncell of elements involved, e.g.{"C","H","O"}:\n');
example= '{''C'', ''H'', ''O''}';
fprintf('Please input cell of elements involved, eg.:%s\n',example)
charset=input('\nPlease input：\n');

while ~feof(fidin)%Number of general parameters
    str=fgetl(fidin);
    str=strtrim(str);
    str=strsplit(str);
    if strcmp('',str)
        fprintf(fidout,'\n');
        break;
    end
    fprintf(fidout,'%10s ',str{1});
    fprintf(fidout,'%-s',str{2});
    if length(str)==2
        fprintf(fidout,'\n');
    elseif length(str)>3
        for i=4:length(str)-1
            fprintf(fidout,' %-s',str{i});
        end
        fprintf(fidout,' %-s\n',str{length(str)});
    elseif length(str)==3
        fprintf(fidout,' %-s\n',str{3});
    end
end

while ~feof(fidin)%Nr of atoms
    str=fgetl(fidin);
    str=strtrim(str);
    str=strsplit(str);
    if strcmp('',str)
        fprintf(fidout,'\n');
        break;
    end
    len=length(str);
    if ismember(str{1},charset)
        for i=1:len
            if i==1
                fprintf(fidout,'%1s',str{i});
            elseif i>1 && i<len
                fprintf(fidout,'%10s',str{i});
            elseif i==len
                fprintf(fidout,'%10s\n',str{i});
            end
        end
    else
        for j=1:len
            if j==1
                fprintf(fidout,'%11s',str{j});
            elseif  j>1 && j<len
                fprintf(fidout,'%10s',str{j});
            elseif j==len
                fprintf(fidout,'%10s\n',str{j});
            end
        end
    end
end

while ~feof(fidin)%Nr of bonds
    str=fgetl(fidin);
    str=strtrim(str);
    str=strsplit(str);
    if strcmp('',str)
        fprintf(fidout,'\n');
        break;
    end
    len=length(str);
    if len==10
        for i=1:len
            if i==1
                fprintf(fidout,'%3s',str{i});
            elseif i==2
                fprintf(fidout,'%3s',str{i});
            elseif i>2 && i<10
                fprintf(fidout,'%10s',str{i});
            elseif i==10
                fprintf(fidout,'%10s\n',str{i});
            end
        end
    elseif len==8
        for i=1:len
            if i==1
                fprintf(fidout,'%16s',str{i});
            elseif i>1 && i<8
                fprintf(fidout,'%10s',str{i});
            elseif i==8
                fprintf(fidout,'%10s\n',str{i});
            end
        end
    end
end

while ~feof(fidin)%! Nr of off-diagonal terms
    str=fgetl(fidin);
    str=strtrim(str);
    str=strsplit(str);
    if strcmp('',str)
        fprintf(fidout,'\n');
        break;
    end
    len=length(str);
    for i=1:len
        if i==1
            fprintf(fidout,'%3s',str{i});
        elseif i==2
            fprintf(fidout,'%3s',str{i});
        elseif i>2 && i<8
            fprintf(fidout,'%10s',str{i});
        elseif i==8
            fprintf(fidout,'%10s\n',str{i});
        end
    end
end

while ~feof(fidin)%!  Nr of angles
    str=fgetl(fidin);
    str=strtrim(str);
    str=strsplit(str);
    if strcmp('',str)
        fprintf(fidout,'\n');
        break;
    end
    len=length(str);
    for i=1:len
        if i==1
            fprintf(fidout,'%3s',str{i});
        elseif i==2
            fprintf(fidout,'%3s',str{i});
        elseif i==3
            fprintf(fidout,'%3s',str{i});
        elseif i>3 && i<10
            fprintf(fidout,'%11s',str{i});
        elseif i==10
            fprintf(fidout,'%11s\n',str{i});
        end
    end
end

while ~feof(fidin)%!  Nr of torsions
    str=fgetl(fidin);
    str=strtrim(str);
    str=strsplit(str);
    if strcmp('',str)
        fprintf(fidout,'\n');
        break;
    end
    len=length(str);
    for i=1:len
        if i==1
            fprintf(fidout,'%3s',str{i});
        elseif i==2
            fprintf(fidout,'%3s',str{i});
        elseif i==3
            fprintf(fidout,'%3s',str{i});
        elseif i==4
            fprintf(fidout,'%3s',str{i});
        elseif i>4 && i<11
            fprintf(fidout,'%11s',str{i});
        elseif i==11
            fprintf(fidout,'%11s\n',str{i});
        end
    end
end

while ~feof(fidin)%!  Nr of hydrogen bonds
    str=fgetl(fidin);
    str=strtrim(str);
    str=strsplit(str);
    if strcmp('',str)
        fprintf(fidout,'\n');
        break;
    end
    len=length(str);
    for i=1:len
        if i==1
            fprintf(fidout,'%3s',str{i});
        elseif i==2
            fprintf(fidout,'%3s',str{i});
        elseif i==3
            fprintf(fidout,'%3s',str{i});
        elseif i>3 && i<7
            fprintf(fidout,'%11s',str{i});
        elseif i==7
            fprintf(fidout,'%11s\n',str{i});
        end
    end
end

fclose(fidin);
fclose(fidout);
fprintf('\nElements with two letters should be manually modified.\n'),
fprintf('\nSuccessfully finished!!!Please check the output.txt.\n');

clear charset fidin fidout i j len str