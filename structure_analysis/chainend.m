%scrit file name chainend
%purpose:
%This program is used to seek for the BO of chain end of polycarbonate and
%corresponding chaintype
%0 for gem-dimethyl benzene radical,1 for benzene end C，2 for gem-dimethyl end C，3 for peracid end O，4 for carbonate end C，5 for phenol end O,
%6 for peracid、carbonate end O，7 for phenol/ether end C
function endBO=chainend(tarBOinformcopy,element,tarelenummatch,ii,tartrajectory)
moleform='';
for i=1:4
    if tarelenummatch{i,2*ii}~=0
        moleform=strcat(moleform,tarelenummatch{i,2*ii-1},num2str(tarelenummatch{i,2*ii}));
    end
end
fprintf('\nThe No. %d in tarelenummatch of trajectory %d, the molecular formula is: %s',tartrajectory{1},ii,moleform);

endBO={};chainendtype=[];
elenum={'C';'H';'O';'N'};
for i=1:4
    for j=1:4
        if strcmp(element(j),elenum{i,1})
            elenum{i,2}=j;
        end
    end
end


[row,~]=size(tarBOinformcopy);
for i=1:row
    if tarBOinformcopy{i,2}==elenum{1,2} && tarBOinformcopy{i,3}==2
        Cbondatom=[0,0];
        found=0;
        for j=1:2
            for k=1:row
                if tarBOinformcopy{i,j+3}==tarBOinformcopy{k,1}
                    Cbondatom(j)=tarBOinformcopy{k,2};found=found+1;
                end
            end
            if found ~=j
                error('原子丢失，其id为%d',tarBOinformcopy{i,j+3});
            end
        end
        if sum(Cbondatom)==2*elenum{1,2}
            [rowend,~]=size(endBO);endBO(rowend+1,:)=tarBOinformcopy(i,:);
            chainendtype(length(chainendtype)+1)=0;
        elseif sum(Cbondatom)==2*elenum{3,2}
            [rowend,~]=size(endBO);endBO(rowend+1,:)=tarBOinformcopy(i,:);
            chainendtype(length(chainendtype)+1)=4;
        end
        
    elseif tarBOinformcopy{i,2}==elenum{1,2} && tarBOinformcopy{i,3}==3
        Cbondatom=[0,0,0];
        found=0;
        for j=1:3
            for k=1:row
                if tarBOinformcopy{i,j+3}==tarBOinformcopy{k,1}
                    Cbondatom(j)=tarBOinformcopy{k,2};found=found+1;
                end
            end
            if found ~=j
                error('原子丢失，其id为%d',tarBOinformcopy{i,j+3});
            end
        end
        
        if sum(Cbondatom)==2*elenum{1,2}+elenum{2,2}
            
            
            Cprimary=[];
            found=0;
            for j=1:3
                for k=1:row
                    if tarBOinformcopy{i,j+3}==tarBOinformcopy{k,1}
                        found=found+1;
                        if tarBOinformcopy{k,2}==elenum{1,2}
                            Cprimary(length(Cprimary)+1)=tarBOinformcopy{k,1};
                        end
                    end
                end
                if found ~=j
                    error('The primary C is loss, whose id is %d',tarBOinformcopy{i,j+3});
                end
            end
            lenCprimary=length(Cprimary);
            found=0;
            for j=1:lenCprimary
                for k=1:row
                    if tarBOinformcopy{k,1}==Cprimary(1,j)
                        Cnum=0;found=found+1;
                        found2=0;
                        for jj=1:tarBOinformcopy{k,3}
                            for kk=1:row
                                if tarBOinformcopy{k,jj+3}==tarBOinformcopy{kk,1}
                                    found2=found2+1;
                                    if tarBOinformcopy{kk,2}==elenum{1,2}
                                        Cnum=Cnum+1;
                                    end
                                end
                            end
                            if found2~=jj
                                error('Atom is not found, whose id is %d',tarBOinformcopy{k,jj+3})
                            end
                        end
                        Cprimary(2,j)=Cnum;
                    end
                end
                if found~=j
                    error('The primary C is loss, whose id is %d',Cprimary(1,j))
                end
            end
            
            if ~ismember(3,Cprimary(2,:))
                Csecondary=[];[~,col]=size(Cprimary);
                found=0;
                for j=1:col
                    for k=1:row
                        if tarBOinformcopy{k,1}==Cprimary(1,j);
                            found=found+1;
                            found2=0;
                            for jj=1:tarBOinformcopy{k,3}
                                if tarBOinformcopy{k,jj+3}~=tarBOinformcopy{i,1}
                                    for kk=1:row
                                        if tarBOinformcopy{k,jj+3}==tarBOinformcopy{kk,1}
                                            found2=found2+1;
                                            if tarBOinformcopy{kk,2}==elenum{1,2}
                                                Csecondary(length(Csecondary)+1)=tarBOinformcopy{kk,1};
                                            end
                                        end
                                    end
                                    if found2~=jj
                                        error('Atom is not found, whose id is %d',tarBOinformcopy{k,jj+3})
                                    end
                                else
                                    found2=found2+1;
                                end
                            end
                        end
                    end
                    if found~=j
                        error('The primary C is loss, whose id is %d',Cprimary(1,j))
                    end
                end
                lenCsecondary=length(Csecondary);
                found=0;
                for j=1:lenCsecondary
                    for k=1:row
                        if tarBOinformcopy{k,1}==Csecondary(1,j)
                            Cnum=0;found=found+1;
                            for jj=1:tarBOinformcopy{k,3}
                                for kk=1:row
                                    if tarBOinformcopy{k,jj+3}==tarBOinformcopy{kk,1} && tarBOinformcopy{kk,2}==elenum{1,2}
                                        Cnum=Cnum+1;
                                    end
                                end
                            end
                            Csecondary(2,j)=Cnum;
                        end
                    end
                    if found~=j
                        error('The secondary C is loss, whose id is %d',Csecondary(1,j));
                    end
                end
                if ~ismember(3,Csecondary(2,:))
                    [rowend,~]=size(endBO);endBO(rowend+1,:)=tarBOinformcopy(i,:);
                    chainendtype(length(chainendtype)+1)=1;
                end
            end
            
        elseif sum(Cbondatom)==3*elenum{1,2}%
            atomid=[];found=0;
            for j=1:3%
                for k=1:row
                    if tarBOinformcopy{k,1}==tarBOinformcopy{i,j+3}
                        found=found+1;
                        for jj=1:tarBOinformcopy{k,3}
                            atomid(length(atomid)+1)=tarBOinformcopy{k,jj+3};
                        end
                    end
                end
                if found~=j
                    error('The primary C is loss, whose id is %d',tarBOinformcopy{i,j+3});
                end
            end
            Hnum=0;
            for j=1:length(atomid)
                for k=1:row
                    if tarBOinformcopy{k,1}==atomid(j) && tarBOinformcopy{k,2}==elenum{2,2}
                        Hnum=Hnum+1;
                    end
                end
            end
            if Hnum>4
                [rowend,~]=size(endBO);endBO(rowend+1,:)=tarBOinformcopy(i,:);
                chainendtype(length(chainendtype)+1)=2;
            end
        end
        
    elseif tarBOinformcopy{i,2}==elenum{3,2} 
        if tarBOinformcopy{i,3}==2
            OHnum=0;
            for j=1:2
                for k=1:row
                    if tarBOinformcopy{k,1}==tarBOinformcopy{i,j+3} && tarBOinformcopy{k,2}==elenum{1,2}
                        OHnum=OHnum+elenum{1,2};
                    elseif tarBOinformcopy{k,1}==tarBOinformcopy{i,j+3} && tarBOinformcopy{k,2}==elenum{2,2}
                        OHnum=OHnum+elenum{2,2};
                    end
                end
            end
            if OHnum==elenum{1,2}+elenum{2,2}
                Onum=0;
                for j=1:2
                    for k=1:row
                        if tarBOinformcopy{k,1}==tarBOinformcopy{i,j+3} && tarBOinformcopy{k,2}==elenum{1,2}
                            for jj=1:tarBOinformcopy{k,3}
                                for kk=1:row
                                    if tarBOinformcopy{kk,1}==tarBOinformcopy{k,jj+3} && tarBOinformcopy{kk,2}==elenum{3,2}
                                        Onum=Onum+1;
                                    end
                                end
                            end
                        end
                    end
                end
                if Onum==3
                    [rowend,~]=size(endBO);endBO(rowend+1,:)=tarBOinformcopy(i,:);
                    chainendtype(length(chainendtype)+1)=3;
                elseif Onum==1
                    [rowend,~]=size(endBO);endBO(rowend+1,:)=tarBOinformcopy(i,:);
                    chainendtype(length(chainendtype)+1)=5;
                end
            end
            
        elseif tarBOinformcopy{i,3}==1
            for j=1:1
                for k=1:row
                    if tarBOinformcopy{k,1}==tarBOinformcopy{i,j+3} && tarBOinformcopy{k,2}==elenum{1,2}
                        Onum=0;
                        for jj=1:tarBOinformcopy{k,3}
                            for kk=1:row
                                if tarBOinformcopy{kk,1}==tarBOinformcopy{k,jj+3} && tarBOinformcopy{kk,2}==elenum{3,2}
                                    Onum=Onum+1;
                                end
                            end
                        end
                    end
                end
                if Onum==1
                    for k=1:row
                        if tarBOinformcopy{k,1}==tarBOinformcopy{i,4} && tarBOinformcopy{k,2}==elenum{1,2}
                            [rowend,~]=size(endBO);endBO(rowend+1,:)=tarBOinformcopy(k,:);
                            chainendtype(length(chainendtype)+1)=7;
                        end
                    end
                elseif Onum==3
                    Oatomid=[];
                    for k=1:row
                        if tarBOinformcopy{k,1}==tarBOinformcopy{i,4} && tarBOinformcopy{k,2}==elenum{1,2}
                            endBOCO=tarBOinformcopy(k,:);
                            for jj=1:tarBOinformcopy{k,3}
                                if tarBOinformcopy{k,jj+3} ~=tarBOinformcopy{i,1} 
                                    for kk=1:row
                                        if tarBOinformcopy{kk,1}==tarBOinformcopy{k,jj+3} && tarBOinformcopy{kk,2}==elenum{3,2}
                                            Oatomid(length(Oatomid)+1)=tarBOinformcopy{kk,1};%
                                        end
                                    end
                                end
                            end
                        end
                    end
                    
                    OCbondnum=0;
                    for k=1:length(Oatomid)
                        for jj=1:row
                            if tarBOinformcopy{jj,1}==Oatomid(k)
                                OCbondnum=OCbondnum+tarBOinformcopy{jj,3};
                            end
                        end
                    end
                    if OCbondnum==3
                        [rowend,~]=size(endBO);endBO(rowend+1,:)=endBOCO;
                        chainendtype(length(chainendtype)+1)=6;
                    else
                        endBOCO={};
                    end
                end
            end
        end
        
    end
end

len=size(endBO,1);
if len>2
    i=1;
    while len
        len=size(endBO,1);
        if len==2
            len=0;
            break;
        end
        if i>size(endBO,1)
            break;
        end
        if endBO{i,2}==elenum{3,2} && endBO{i,3}==3
            endBO(i,:)=[];
            chainendtype(i)=[];
            len=len-1;
        elseif endBO{i,2}==elenum{3,2} && endBO{i,3}==2
            for j=1:size(tarBOinformcopy,1)
                for k=1:2
                    if tarBOinformcopy{j,1}==endBO{i,k+3}
                        if tarBOinformcopy{j,2}==elenum{1,2} && tarBOinformcopy{j,3}==4
                            endBO(i,:)=[];
                            chainendtype(i)=[];
                            len=len-1;
                            break;
                        elseif tarBOinformcopy{j,2}==elenum{1,2} && tarBOinformcopy{j,3}==3
                            for kk=1:size(tarBOinformcopy,1)
                                for jj=1:3
                                    if tarBOinformcopy{j,3+jj}==tarBOinformcopy{kk,1} && tarBOinformcopy{kk,2}==elenum{1,2}
                                        endBO(i,:)=[];
                                        chainendtype(i)=[];
                                        len=len-1;
                                        break;
                                    end
                                end
                            end
                        else
                            i=i+1;
                            if i>size(endBO,1)
                                break;
                            end
                            len=len-1;
                        end
                    end
                end
            end
        else
            i=i+1;
            if i>size(endBO,1)
                break;
            end
            len=len-1;
        end
        if i>size(endBO,1)
            break;
        end
    end
end




for i=1:length(chainendtype)
    endBO{i,8}=chainendtype(i);
end

if length(chainendtype)~=2
    endBO
    fprintf('\nPlease try to remove small molecules, like water, oxygen');
    error('\nSee above, the found end atom is %d, is not 2, please check it!',length(chainendtype));
else
    fprintf('\nSuccessfully find the end atoms, whose type is saved in chainendtype, bond orser si saved in tarBOinformcopy');
    fprintf('\n0 for gem-dimethyl benzene radical,\n1 for benzene end C，\n2 for gem-dimethyl end C，\n3 for peracid end O，\n4 for carbonate end C，\n5 for phenol end O, \n6 for peracid、carbonate end O，\n7 for phenol/ether end C\n');
end







