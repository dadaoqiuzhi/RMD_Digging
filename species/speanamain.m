%scrit file name speanamain
%purpose:This is the main program used to analysis species file
%include functions:species_analysis.species_capture, species_classfy and
%statistics£¬some orther invoked auxiliary function like char26cor, delnull,
%membercheck and molecuweight.
%version 1;2018.6.25
disp('##################################################################################################################################')
disp('Welcome!--by Qiang Liu @Institute of Nuclear Physics and Chemistry, China Academy of Engineering Physics; Email: liubinqiang@163.com');
disp('Repository adress of the Source code on github: https://github.com/dadaoqiuzhi/RMD_Digging');
disp('References: 1.Fuel 287 (2021) 119484. 2.ACS Appl. Mat. Interfaces 13(34) (2021) 41287-41302. More work is coming!')
disp('##################################################################################################################################')
fprintf('\nThis procedure integrate the following programs: \n1.species_analysis\n2.species_capture\n')
fprintf('3.species_classfy\n4.statistics\n\n')

while true
    progchoi=input('Please select the program No.: 1,2,3 or 4: \n');
    switch progchoi
        case 1
            species_analysis
        case 2
            species_capture
        case 3
            species_classfy
        case 4
            statistics
        otherwise
            disp('Illegal program No, please input the correct program No.: 1,2,3 or 4')
    end
    fprintf('\nThis procedure integrate the following programs: \n1.species_analysis\n2.species_capture\n')
    fprintf('3.species_classfy\n4.statistics\n\n')
    fprintf('\n\nCtrl+c to exit\n')
    clear progchoi
end


