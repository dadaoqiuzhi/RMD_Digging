%scrit file name anaoptlogic
%purpose:
%This program is used to analyze the logic value (true or false) of a given
%moleculae at the predetermined request
function logicvalue=anaoptlogic(anaoption,ii,tarelenummatch,Cnumrequest,Onumrequest,Nnumrequest,Hnumrequest,Totalnumrequest,MWrequest)
count=0;logicvalue=1;
if ismember(1,anaoption)
    count=count+1;
    for i=1:4
        if strcmp(tarelenummatch{i,2*ii-1},'C')
            Cnumlogic=logicinspect(Cnumrequest,tarelenummatch,i,ii);
            logicvalue=logicvalue*Cnumlogic;
        end
    end
else
    Cnumlogic=1;
    logicvalue=logicvalue*Cnumlogic;
end

if ismember(2,anaoption)
    count=count+1;
    for i=1:4
        if strcmp(tarelenummatch{i,2*ii-1},'O')
            Onumlogic=logicinspect(Onumrequest,tarelenummatch,i,ii);
            logicvalue=logicvalue*Onumlogic;
        end
    end
else
    Onumlogic=1;
    logicvalue=logicvalue*Onumlogic;
end

if ismember(3,anaoption)
    count=count+1;
    for i=1:4
        if strcmp(tarelenummatch{i,2*ii-1},'N')
            Nnumlogic=logicinspect(Nnumrequest,tarelenummatch,i,ii);
            logicvalue=logicvalue*Nnumlogic;
        end
    end
else
    Nnumlogic=1;
    logicvalue=logicvalue*Nnumlogic;
end

if ismember(4,anaoption)
    count=count+1;
    for i=1:4
        if strcmp(tarelenummatch{i,2*ii-1},'H')
            Hnumlogic=logicinspect(Hnumrequest,tarelenummatch,i,ii);
            logicvalue=logicvalue*Hnumlogic;
        end
    end
else
    Hnumlogic=1;
    logicvalue=logicvalue*Hnumlogic;
end

if ismember(5,anaoption)%Ãÿ ‚
    count=count+1;
    sumatom=cell2mat(tarelenummatch(:,ii*2));
    sumatom=sum(sumatom);sumatomcell{1}=0;
    sumatomcell{2}=sumatom(1);
    Totalnumlogic=logicinspect(Totalnumrequest,sumatomcell,1,1);
    logicvalue=logicvalue*Totalnumlogic;
else
    Totalnumlogic=1;
    logicvalue=logicvalue*Totalnumlogic;
end

if ismember(6,anaoption)
    count=count+1;
    classmatch=tarelenummatch(:,ii*2-1:ii*2);
    MW=molecuweight(classmatch);MWcell{1}=0;
    MWcell{2}=MW(1);
    MWlogic=logicinspect(MWrequest,MWcell,1,1);
    logicvalue=logicvalue*MWlogic;
else
    MWlogic=1;
    logicvalue=logicvalue*MWlogic;
end

if count ~= length(anaoption)
    error('\nNo. for molecular chain is wrong, please check it!')
end

%logicvalue=Cnumlogic*Onumlogic*Nnumlogic*Hnumlogic*Totalnumlogic*MWlogic;

