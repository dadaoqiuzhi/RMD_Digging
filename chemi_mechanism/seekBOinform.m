if choi==2
    reactname=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','reactantBO','.txt');
elseif choi==4
    reactname=strcat(species{1},'-',num2str(tartrajectory{1,1}),'-',num2str(tartrajectorycopy(1)),'-',num2str(pronum),'-','productBO','.txt');
elseif choi==1
    producttname=strcat(species{1},'-',num2str(frame{1,irow}),'-',num2str(frame{2,irow}),'-','productBO','.txt');
    tarBOinformcopy3=tarBOinform;
end
[row3,~]=size(tarBOinformcopy3);
if choi==2 || choi==4
    fidBO=fopen(reactname,'wt');
elseif choi==1
    fidBO=fopen(producttname,'wt');
end

fprintf(fidBO,'%-12s%-12s%-12s%-12s%-12s%-12s%-12s\n','id','atom type','bonds','1bo','2bo','3bo','4bo');
for i=1:row3
    for j=1:7
        if ~ischar(tarBOinformcopy3{i,j})
            fprintf(fidBO,'%-12d',tarBOinformcopy3{i,j});
        elseif ischar(tarBOinformcopy3{i,j})
            fprintf(fidBO,'%-12s',tarBOinformcopy3{i,j});
        end
        if j==7
            fprintf(fidBO,'\n');
        end
    end
end
fclose(fidBO);
if choi==2
    fprintf('\nBO information of Group %d reactants is successfully exported.',tartrajectory{1,1});
elseif choi==4
    fprintf('\nBO information of Group %d products is successfully exported.',tartrajectory{1,1});
end
if choi==1
    fprintf('\nBO information is successfully exported.');
end
