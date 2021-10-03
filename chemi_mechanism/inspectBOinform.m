if choi==2
    reactname=strcat(species{1},'-','any','-',num2str(tartrajectorycopy(1)),'-','1','-','productBO','.txt');
elseif choi==4
    reactname=strcat(species{1},'-','any','-',num2str(tartrajectorycopy(1)),'-','1','-','reactantBO','.txt');
end
[row6,~]=size(tarBOinformcopy2);
fidBO=fopen(reactname,'wt');
fprintf(fidBO,'%-12s%-12s%-12s%-12s%-12s%-12s%-12s\n','id','atom type','bonds','1bo','2bo','3bo','4bo');
for i=1:row6
    for j=1:7
        if ~ischar(tarBOinformcopy2{i,j})
            fprintf(fidBO,'%-12d',tarBOinformcopy2{i,j});
        elseif ischar(tarBOinformcopy2{i,j})
            fprintf(fidBO,'%-12s',tarBOinformcopy2{i,j});
        end
        if j==7
            fprintf(fidBO,'\n');
        end
    end
end
fclose(fidBO);
if choi==2
    fprintf('\nBO information of products in group %d is exported',tartrajectorycopy(1));
elseif choi==4
    fprintf('\nBO information of reactants in group %d is exported',tartrajectorycopy(1));
end