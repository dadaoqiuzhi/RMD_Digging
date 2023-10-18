%scrit file name bondorder_capture
%purpose:
%This program is used to analyze  bond order information of simple molecule
%in a specific trajectory.
%version 1;2018.6.26
disp('##################################################################################################################################')
disp('Welcome!--by Qiang Liu @Institute of Nuclear Physics and Chemistry, China Academy of Engineering Physics; Email: liubinqiang@163.com');
disp('Repository adress of the Source code on github: https://github.com/dadaoqiuzhi/RMD_Digging');
disp('References: 1.Fuel 287 (2021) 119484. 2.ACS Appl. Mat. Interfaces 13(34) (2021) 41287-41302. 3.ACS Appl. Mat. Interfaces 2022, 14.(4), 5959-5972.')
disp('4.ACS Materials Letters 2023, 2174-2188. More work is coming!')
disp('##################################################################################################################################')
fprintf('This program will analysis the BO information in bondoutdata of specified trajectory, helping to obtain the structure information')
fprintf('\nRepeated refine should considered for some complex case and the obtained last BO information should overwrite these in bondoutdata in advance')

speciestrjnum=input('\nPlease input the trajectory timestep. It can be obtained from the analysis results of species files: \n');
disp('atomid,bondnum and lpnum can be used to anasylyze the species with known atom id, especially from the preliminary analysis.')
atomid=input('\nPlease input the atom id of the characteristic species (positive integer). Other input can be given accordingly. \n"ignore" to skip this input: \n','s');
atomid=lower(atomid);
fprintf('\nCharacteristic atom type of the species: \n1.C-type of CH4, C-type of C2H4, H-type of H2, O-type of H2O, C-type or O-typ of CO')
fprintf('2.C-type or O-type of CH3OH, C2-type, O1-type and O2-type of CH3COOH, O1-type of CH3COOCH2CH3')
atomtype=input('\nPlease input the characteristic atom type options, a positive integer must be given: \n');
atomtype=lower(atomtype);
bondnum=input('\nPlease input the chemical bond number (positive integer). "ignore" to skip this input\n','s');
bondnum=lower(bondnum);
lpnum=input('\nPlease input the lone pair electrons number of the characteristic atom(e.g.O or N),decimals is possiple, "ignore" to skip this input\n','s');
lpnum=lower(lpnum);

[row,col]=size(bondoutdata);
bondrownum=0;
for i=1:row
    if strcmp(bondoutdata{i,1},'Timestep') 
        if bondoutdata{i,2}==speciestrjnum
            tartrjnum=i;%Timestep num line
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
    tarbondnum=tarbondnum-tartrjnum-1;%start-stop line of bond order
else
    tarbondnum=row-tartrjnum;%
end
        
bondcapturedata={};bondorderdata={};
if ~strcmp(atomid,'ignore')%with explicit tomid,only one result
    taratomidrow=atomidfind(atomid,tartrjnum,tarbondnum,bondoutdata);
    bonddatacapture(1,:)=bondoutdata(taratomidrow,:);
    bondorderdata(1,:)=bonddatacapture(1,:);
    
    
    
else%without atomid, must know atomtype, many results in general
    taratomtyperow=atomtypefind(atomtype,tartrjnum,tarbondnum,bondoutdata);
    for i=1:length(taratomtyperow)
        for j=1:col
            bonddatacapture{i,j}=bondoutdata{taratomtyperow(i),j};
        end
    end%cell data for known atom type, generally based on the known atom type to start work
    if ~strcmp(bondnum,'ignore')%known bond number, refine the results after atomtype search
        tarabondnum=bondnumfind(bondnum,bonddatacapture);%cell for atom type data with known bond number
        if ~strcmp(lpnum,'ignore')%known lone pair electrons number
            tarlpnum=lpnumfind(lpnum,tarabondnum);%final result,unknown atom id, known atom type, known bondnum, known lpnum
            bondorderdata(:,:)=tarlpnum(:,:); 
        else%unknown lone pair electrons number
            bondorderdata(:,:)=tarabondnum(:,:);
        end
    else%unknown bond number
        if ~strcmp(lpnum,'ignore')%known lone pair electrons number, refine the result after bondnumsearch
            tarlpnum=lpnumfind(lpnum,bonddatacapture);
            bondorderdata(:,:)=tarlpnum(:,:);
        else%nuknown lone pair electrons number
            bondorderdata(:,:)=bonddatacapture(:,:);
        end
    end
end
fprintf('\nsigma bond is formed with bond length less than 1.5, sigma bond is broken with bond length larger than 2.5\n(2) the first pi bond is formed with bond length less than 1.2, the first pi bond is broken with bond length larger than 1.75\n')
fprintf('(3)the second pi bond is formed with bond length less than 1.0, the second pi bond is broken with bond length larger than 1.4\n(4)maximum BO of C-C is 3, maximum BO of C-H and H-H is 1 (only sigma bond contribution)\n')
disp('bondorder_capture is successfully finished, and results is saved in bondorderdata')