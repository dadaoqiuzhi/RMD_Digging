%scrit file name delnonbackbone
%purpose:
%This program is used to rebuild backbone of PC with nonsignificant atom
%delated and coordinates of o- & m-position C atoms in bezene simplified
function BOrderid=order_idtag(endBO,backboneBO,delatomid,element)
BOrderid=[];bondid=[];bondidtemp=[];
row=size(backboneBO,1);
BOrderid(1,1)=1;
BOrderid(1,2)=endBO{1,1};
for i=1:row
    if backboneBO{i,1}==BOrderid(1,2)
        bondid(1:backboneBO{i,3})=cell2mat(backboneBO(i,4:3+backboneBO{i,3}));
    end
end
len=length(bondid);
while len
    if ismember(bondid(len),delatomid)
        bondid(len)=[];
        len=len-1;
    else
        len=len-1;
    end
end

num=2;
while ~isempty(bondid)
    if length(bondid)>2
        error('\nThe atom number directly bonded to the same level atom on the main chain is larger than 2, maybe caused by small molecules bonded to the main chain, please check it!')
    end
    BOrderid(num,1)=num;
    len=length(bondid);
    BOrderid(num,2:1+len)=bondid(1,:);
    bondidcopy=[];bondidcopy=bondid;
    bondid=[];
    for i=1:len
        for j=1:row
            if bondidcopy(i)==backboneBO{j,1}
                for k=1:backboneBO{j,3}
                    if ~ismember(backboneBO{j,k+3},delatomid)
                        zeroruleout=BOrderid(num-1,2:end);
                        zeroruleout(find(zeroruleout==0))=[];
                        if ~ismember(backboneBO{j,k+3},zeroruleout)
                            if ~ismember(backboneBO{j,k+3},bondidtemp)
                                bondidtemp(length(bondidtemp)+1)=backboneBO{j,3+k};
                            end
                        end
                    end
                end
            end
        end
    end
    num=num+1;
    bondid=bondidtemp;
    bondidtemp=[];
end



elenum={'C';'H';'O';'N'};
for i=1:4
    for j=1:4
        if strcmp(element(j),elenum{i,1})
            elenum{i,2}=j;
        end
    end
end

for i=1:size(BOrderid,1)
    if BOrderid(i,3)~=0
        for j=2:3
            for k=1:size(backboneBO,1)
                if BOrderid(i,j)==backboneBO{k,1} && backboneBO{k,2}==elenum{3,2}
                    BOrderid(i,j)=0;
                end
            end
        end
        if BOrderid(i,2)*BOrderid(i,3)==0
            BOrderid(i,2)=max(BOrderid(i,2:3));
            BOrderid(i,3)=0;
        end
    end
    
end

if size(BOrderid,2)~=3
    error('\nThe ordered id for main chain is not equal to 3, please check it!')
end
