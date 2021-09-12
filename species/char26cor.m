%scrit file name char26cor
%purpose:
%mapping number systems to column numbers in excel
function charcor=char26cor(matches)
charcor=char();
datadelimiter={'1','2','3','4','5','6','7','8','9','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','0'};
char26={'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};
for i=length(matches):-1:1
    if matches{i}==0
        charcor(i)='Z';
        for j=1:length(datadelimiter)
            if datadelimiter(j)==matches{i-1}
                charcor(i-1)=char26{j-1};
            end
        end
    else
        for j=1:length(datadelimiter)
            if matches{i}==datadelimiter{j}
                charcor(i)=char26{j};
            end
        end
    end
end
        