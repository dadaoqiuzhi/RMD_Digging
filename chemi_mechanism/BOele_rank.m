%scrit file BOele_rank
%purpose:
%sort the found elements in BO information and used to obtain BO from tarBOinform
function id=BOele_rank(temporary,tarBOinform,elepriority)
temporarycopy={};
temporarycopy{1}=temporary{2};temporarycopy{2:5}=temporary{4:7};
temporary2=BOid_element(temporarycopy,tarBOinform);
temporary3=temporary2(1,2:5);
temporary3(2,:)=temporary(1,4:7);
for i=1:4
    if ischar(temporary3{1,i})
        temporary3(:,i)=[];
    end
end
[~,col]=size(elepriority);
[~,col2]=size(temporary3);
if ~isempty(temporary3)
    for i=1:col2
        for j=1:col
            if temporary3{1,i}==elepriority{3,j}
                temporary3{3,i}=elepriority{2,j};
            end
        end
    end
end
temporary3=cell2mat(temporary3);
temporary3=temporary3';
temporary3=sortrows(temporary3,3);

for i=1:4
    elesame=ismember(temporary3(1,:),i);
    sumelesame=sum(elesame);
    if sumelesame>1
        sumelesame=0;
        for j=1:length(elesame)                
            if elesame(j)~=0
                sumelesame=sumelesame+1;
            end
            if sumelesame==1
                startpo=j;
            elseif sumelesame==sumelesame
                endpo=j;
                break
            end
        end
        
    end
end

lenid=length(id);
if lenid<5
    for i=lenid+1:5
        id{lenid+i}='Nan';
    end
end
    



