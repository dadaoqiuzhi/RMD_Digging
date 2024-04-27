%scrit file name logicinspect
%purpose:
%This program is used to help analyze the logic value (true or false) of a given
%moleculae at the predetermined request
function logicvalue=logicinspect(numrequest,tarelenummatch,i,ii)
if length(numrequest)==2
    if strcmp(numrequest{1},'<')
        if numrequest{2}>tarelenummatch{i,ii*2}
            logicvalue=1;
        else
            logicvalue=0;
        end
    elseif strcmp(numrequest{1},'>')
        if numrequest{2}<tarelenummatch{i,ii*2}
            logicvalue=1;
        else
            logicvalue=0;
        end
    elseif strcmp(numrequest{1},'=')
        if numrequest{2}==tarelenummatch{i,ii*2}
            logicvalue=1;
        else
            logicvalue=0;
        end
    else
        logicvalue=0;
        error('\nIllegal logic operator, possibly caused by lossed white space!');
    end
elseif length(numrequest)==3
    if strcmp(numrequest{2},'=')
        if strcmp(numrequest{1},'<')
            if numrequest{3}>=tarelenummatch{i,ii*2}
                logicvalue=1;
            else
                logicvalue=0;
            end
        elseif strcmp(numrequest{1},'>')
            if numrequest{3}<=tarelenummatch{i,ii*2}
                logicvalue=1;
            else
                logicvalue=0;
            end
        else
            logicvalue=0;
            error('\nCombined operator? But the first operator is not a < or > operator, please check it!');
        end  
    else
        logicvalue=0;
        error('\nError triggered for criterion of molecular chain limitation, = operator is not found, please check it!');
    end
elseif length(numrequest)==1
    if numrequest{1}==1;
        logicvalue=1;
    else
        logicvalue=0;
        error('\nNo limited criterion for molecular chain, but initialization does not set cell array with one or loss white space, please check it');
    end
else
    logicvalue=0;
    error('\nError triggered for criterion of molecular chain limitation, please check it!');
end
