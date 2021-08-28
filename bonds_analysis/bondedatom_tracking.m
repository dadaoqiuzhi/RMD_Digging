%scrit file name bondedatom_tracking
%purpose:
%This program is used to tarck the specified atom with atom id during
%degradation of macromolecule.
%include bonds_analysis, bondorder_deepmining and bond_classify program
%version 1;2018.6.30
disp('Welcome!--by Qiang Liu @Institute of Nuclear Physics and Chemistry, China Academy of Engineering Physics; Email: liubinqiang@163.com');
disp('Repository adress of the Source code on github: https://github.com/dadaoqiuzhi/RMD_Digging');

fprintf('\nThis program is intended to track the specific atom in a species, the BO information is returned\n')
atomid=input('\nPlease input the interested atom id: \n');
fprintf('\nbondedatom_tracking is running, please wait...')

[tarbonewrow,~]=size(tarBOinformnew);
for i=1:tarbonewrow%find the atom id in tarBOinformnew
    if atomid==tarBOinformnew{i,1}
        tarrow=i;
        break
    end
end
[dataclassrow,~]=size(dataclass);
for i=1:dataclassrow%confirm the block of the species (in dataclass) with the interested atom id in tarBOinformnew 
    if dataclass{i,3}<=tarrow && dataclass{i,4}>=tarrow
        tarrowtwo=i;
        inirow=dataclass{i,3};endrow=dataclass{i,4};
    end
end

atombondmolecule={};%export the data
atombondmolecule{1,1}=dataclass{tarrowtwo,1};%export the molecular formula 
atomidrow=tarrow-inirow+2;
atombondmolecule{1,2}=atomidrow;
for j=3:15
    atombondmolecule{1,j}=[];
end
rownum=endrow-inirow+1;%
atombondmolecule(2:rownum+1,:)=tarBOinformnew(inirow:endrow,:);
fprintf('\n\nbondedatom_tracking is successfully finished\n')
fprintf('\nthe molecular formula and BO information with the interested atom id is saved in atombondmolecule, the start-stop line of BO is shown here\n')

clear atomid atomidrow dataclassrow endrow i j inirow rownum tarbonewrow tarrow tarrowtwo