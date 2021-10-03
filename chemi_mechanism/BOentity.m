%scrit file BOentity
%purpose:
%This program is used to generate the instance of all possible BOtable rows
%instantiation of BOtable
function BOtablerowset=BOentity(BOtablerow)
BOtablerowset={};starnum=0;starlocate=[];
for i=2:5
    if strcmp(BOtablerow{i},'*')
        starnum=starnum+1;
        lenstarlocate=length(starlocate);
        starlocate(lenstarlocate+1)=i; 
    end
end
Nannum=0;
for j=2:5
    if strcmp(BOtablerow(j),'Nan')
        Nannum=Nannum+1;
    end
end
if starnum==0
    sum=0;product=1;
    BOtablerowset=BOtablerow;
    for i=2:5
        if ~isletter(BOtablerowset{i})
            sum=sum+BOtablerowset{i};
            product=product*BOtablerowset{i};
        end
    end
    BOtablerowset{6}=Nannum;
    BOtablerowset{7}=sum;
    BOtablerowset{8}=product;
    
else
    A=nchoosek({1,2,3,4},starnum);
    %B=reshape(A(:,perms(1:4)),[],starnum);
    [row,~]=size(A);
    for j=1:row
        BOtablerowset(j,:)=BOtablerow;
    end
    lenstarlocate=length(starlocate);
    for j=1:row
        for k=1:lenstarlocate
            BOtablerowset{j,starlocate(k)}=A{j,k};
        end
    end
    
    for i=1:row
        sum=0;product=1;
        for j=2:5
            if ~isletter(BOtablerowset{i,j})
                sum=sum+BOtablerowset{i,j};
                product=product*BOtablerowset{i,j};
            end
        end
        BOtablerowset{i,6}=Nannum;
        BOtablerowset{i,7}=sum;
        BOtablerowset{i,8}=product;
    end
end

