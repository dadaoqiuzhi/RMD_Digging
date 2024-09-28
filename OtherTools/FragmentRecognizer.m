%FragmentRecognizer
%purpose: to recongnize the species fragments in the car/arc file seperated
%by end keyword.

fprintf('\nThis program distinguishes the fragment line numbers of each species based on the "end" keyword in the car and arc files')
fprintf('and obtains the atomic numbers after converting them into LAMMPS data files');
FileName = input('\nPlease enter the file name for processing: \n','s');
tic 
disp('FragmentRecognizer program is running, please wait')
rawdata=fopen(FileName,'r');
block_start = 0;
AtomNo = [];
line = 0;
while ~feof(rawdata)
    dataline=fgetl(rawdata);
    datacell=textscan(dataline,'%s','delimiter','\n');
    datacellchar=char(datacell{1});
    datarep=strtrim(datacellchar);
    datasplit=strsplit(datarep);
    if strcmpi(datasplit(1),'PBC')
        block_start = 1;
    end
    
    if block_start == 1
        while ~feof(rawdata)
            dataline=fgetl(rawdata);
            datacell=textscan(dataline,'%s','delimiter','\n');
            datacellchar=char(datacell{1});
            datarep=strtrim(datacellchar);
            datasplit=strsplit(datarep);
            if isscalar(datasplit)
                if strcmpi(datasplit(1),'end')
                    AtomNo(size(AtomNo,1)+1,1) = line;
                end
            else
                line = line + 1;
            end
        end
    end
end

for i = 1:size(AtomNo,1)
    if i == 1
        AtomNo(i,2) = AtomNo(i,1);
        AtomNo(i,1) = 1;
    else
        AtomNo(i,2) = AtomNo(i,1);
        AtomNo(i,1) = AtomNo(i-1,2) + 1;
    end
end
fprintf('\n\nFragmentRecognizer program is finished\n')
fprintf('\nData are stored in AtomNo\n\n')
Elapsedtime = toc; 
fprintf('\nDuration of this run: %.2f s\n',Elapsedtime)

clear block_start datacell datacellchar dataline datarep datasplit FileName i line rawdata Elapsedtime