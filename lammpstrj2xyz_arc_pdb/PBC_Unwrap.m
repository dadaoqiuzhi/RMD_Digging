%scrit file name PBC_Unwrap
%purpose:
%This program is used to unwrap the coordinates of atoms sue to PBC
%condition
%version 1;2021.10.22

function trjdata=PBC_Unwrap(unwrapans,tarBOinform,trjdata,BOXsize,boxsize,element)
if strcmpi(unwrapans,'y')
    load('BondRadii.mat');
end
fprintf('\nPerform coordinate unwrap, aiming to treat the ghost atom out of the box for a continuous display, please wait ...')
tarBOinform_copy=tarBOinform;
trjdata_copy=trjdata;
if isempty(tarBOinform_copy)
    error('No BO information can be used to be processed');
end

% count1=size(tarBOinform_copy,1);
% count2=1;
% while count1
%     if strcmpi(tarBOinform_copy{count2,1},'#')
%         tarBOinform_copy(count2,:)=[];
%         count1=count1-1;
%     else
%         count2=count2+1;
%         count1=count1-1;  
%     end
% end


RadiiData=[];
fprintf('\nReference literature for bond length: ')
fprintf('\n1.Pauling, L.; Pauling, P. Chemistry. San Francisco: W. H. Freeman Company, 1975')
fprintf('\n2.Bondi, A. J. Phys. Chem., 1964, 68: 441')
fprintf('\n3.Bokii, G. B. Kristallokhimiya ( Crystal Chemistry ). Moscow: Nauka, 1971')
fprintf('\n4.Allinger, N. L.; Zhou, X.; Bergsma, J. J. Mol. Struct., Theochem., 1994, 312:69')
fprintf('\n5.Zefirov, Yu V. Russ. J. Inorg. Chem., 2000, 45: 1552')
fprintf('\n6.Batsanov, S. S. Russ. J. Inorg. Chem., 1991, 36: 1694; Inorg. Mater., 2001, 37: 871')
for i=1:length(element)
    for j=1:length(element)
        if i<=j
            count5=size(RadiiData,1);
            RadiiData(count5+1,1)=i;
            RadiiData(count5+1,2)=j;
            ii=3;
            for k=1:size(BondRadii,2)
                if strcmpi(element{i},BondRadii{1,k})
                    RadiiData(count5+1,ii)=BondRadii{2,k};
                    ii=ii+1;
                elseif strcmpi(element{j},BondRadii{1,k})
                    RadiiData(count5+1,ii)=BondRadii{2,k};
                    ii=ii+1;
                end
            end
            if i==j
                RadiiData(count5+1,4)=RadiiData(count5+1,3);
            end
        end
    end
end
RadiiData(:,5)=(RadiiData(:,3)+RadiiData(:,4))*1.15/2;

%
species_table=[];count1=0;
for i=1:size(tarBOinform_copy,1)
    if strcmpi(tarBOinform_copy(i,1),'#')
        count1=count1+1;
        species_table(count1,1)=count1;
        if count1==1
            species_table(count1,2)=1;
        else
            species_table(count1,2)=species_table(count1-1,3)+1;
        end
        species_table(count1,3)=i;
    end
end



if strcmpi(BOXsize,'y')
    fprintf('\nScaled coordinate exists, real coordinate will be calculated before performing unwrap')
    count1=size(species_table,1);
    count2=1;
    while count2<=count1
        atomPBC=[];atomPBCmember=[];tarBOinform_copyi={};trjdata_temp=[];
        tarBOinform_copyi(:,:)=tarBOinform_copy(species_table(count2,2):species_table(count2,3)-1,:);
        
        atomPBC(1,1)=tarBOinform_copyi{1,1};
        atomPBC(1,2)=tarBOinform_copyi{1,2};
        atomPBC(1,3)=tarBOinform_copyi{1,3};
        atomPBCmember(length(atomPBCmember)+1,1)=atomPBC(1,1);
        if atomPBC(1,3)>0
            for i=1:atomPBC(1,3)
                atomPBC(1,i+3)=tarBOinform_copyi{1,i+3};
            end
            [atomPBC,atomPBCmember]=BondLink(atomPBC,atomPBCmember,tarBOinform_copyi); 
        end
        
        while ~isempty(atomPBC)
            %
            xyzdata=[];
            for i=1:size(trjdata_copy,1)%
                if trjdata_copy(i,1)==atomPBC(1,1)
                    xyzdata(1,1)=atomPBC(1,1);
                    xyzdata(1,2)=atomPBC(1,2);%
                    xyzdata(1,3)=boxsize(1,1)+trjdata_copy(i,3)*(boxsize(1,2)-boxsize(1,1));%X
                    xyzdata(1,4)=boxsize(2,1)+trjdata_copy(i,4)*(boxsize(2,2)-boxsize(2,1));%Y
                    xyzdata(1,5)=boxsize(3,1)+trjdata_copy(i,5)*(boxsize(3,2)-boxsize(3,1));%Z
                end
            end
            if atomPBC(1,3)>0
                for i=1:atomPBC(1,3)
                    for j=1:size(trjdata_copy,1)
                        if trjdata_copy(j,1)==atomPBC(1,i+3)
                            count4=size(xyzdata,1)+1;
                            xyzdata(count4,1)=atomPBC(1,i+3);
                            xyzdata(count4,2)=trjdata_copy(j,2);
                            xyzdata(count4,3)=boxsize(1,1)+trjdata_copy(j,3)*(boxsize(1,2)-boxsize(1,1));%X
                            xyzdata(count4,4)=boxsize(2,1)+trjdata_copy(j,4)*(boxsize(2,2)-boxsize(2,1));%Y
                            xyzdata(count4,5)=boxsize(3,1)+trjdata_copy(j,5)*(boxsize(3,2)-boxsize(3,1));%Z
                        end
                    end
                end
            end
            %
            if size(xyzdata,1)>1
                xyzdata=BondForm(xyzdata,RadiiData,boxsize);%
            end
            %             trjdata_temp(size(trjdata_temp,1)+1:size(trjdata_temp,1)+size(xyzdata,1),:)=xyzdata(:,:);
            for i=1:size(xyzdata,1)%
                if isempty(trjdata_temp)
                    trjdata_temp(size(trjdata_temp,1)+1,:)=xyzdata(i,:);
                else
                    if ~ismember(xyzdata(i,1),trjdata_temp(:,1))
                        trjdata_temp(size(trjdata_temp,1)+1,:)=xyzdata(i,:);
                    end
                    
                end
            end
            %
            atomPBC(1,:)=[];
            if ~isempty(atomPBC)
                [atomPBC,atomPBCmember]=BondLink(atomPBC,atomPBCmember,tarBOinform_copyi); %
            end
        end
        %
        xmean=mean(trjdata_temp(:,3));
        if xmean<-(boxsize(1,2)-boxsize(1,1))/2
            trjdata_temp(:,3)=trjdata_temp(:,3)+(boxsize(1,2)-boxsize(1,1));
            fprintf('\nPBC transformation according to the X coordinate deviation');
        elseif xmean>(boxsize(1,2)-boxsize(1,1))*1.5
            trjdata_temp(:,3)=trjdata_temp(:,3)-(boxsize(1,2)-boxsize(1,1));
            fprintf('\nPBC transformation according to the X coordinate deviation');
        end
        ymean=mean(trjdata_temp(:,4));
        if ymean<-(boxsize(2,2)-boxsize(2,1))/2
            trjdata_temp(:,4)=trjdata_temp(:,4)+(boxsize(2,2)-boxsize(2,1));
            fprintf('\nPBC transformation according to the Y coordinate deviation');
        elseif ymean>(boxsize(2,2)-boxsize(2,1))*1.5
            trjdata_temp(:,4)=trjdata_temp(:,4)-(boxsize(2,2)-boxsize(2,1));
            fprintf('\nPBC transformation according to the Y coordinate deviation');
        end
        zmean=mean(trjdata_temp(:,5));
        if zmean<-(boxsize(3,2)-boxsize(3,1))/2
            trjdata_temp(:,5)=trjdata_temp(:,5)+(boxsize(3,2)-boxsize(3,1));
            fprintf('\nPBC transformation according to the Z coordinate deviation');
        elseif zmean>(boxsize(2,2)-boxsize(2,1))*1.5
            trjdata_temp(:,5)=trjdata_temp(:,5)-(boxsize(3,2)-boxsize(3,1));
            fprintf('\nPBC transformation according to the Z coordinate deviation');
        end
        for i=1:size(trjdata_temp,1)
            for j=1:size(trjdata,1)
                if trjdata_temp(i,1)==trjdata(j,1)
                    trjdata(j,3:5)=trjdata_temp(i,3:5);
                end
            end
        end
        
        fprintf('\nSuccessfully unwrap the %d th species, %d species need unwrap',count2,count1)
        count2=count2+1;%
    end
    fprintf('\nGood job! The coordnate of %d species in total is unwraped',size(species_table,1))

    
    
    
elseif strcmpi(BOXsize,'n')
    fprintf('\nThe coordinate is not scaled, directly to perform unwrap')
    count1=size(species_table,1);
    count2=1;
    while count2<=count1
        atomPBC=[];atomPBCmember=[];tarBOinform_copyi={};trjdata_temp=[];
        tarBOinform_copyi(:,:)=tarBOinform_copy(species_table(count2,2):species_table(count2,3)-1,:);%
        %
        atomPBC(1,1)=tarBOinform_copyi{1,1};%
        atomPBC(1,2)=tarBOinform_copyi{1,2};%
        atomPBC(1,3)=tarBOinform_copyi{1,3};%
        atomPBCmember(length(atomPBCmember)+1,1)=atomPBC(1,1);%
        if atomPBC(1,3)>0
            for i=1:atomPBC(1,3)
                atomPBC(1,i+3)=tarBOinform_copyi{1,i+3};
            end
            [atomPBC,atomPBCmember]=BondLink(atomPBC,atomPBCmember,tarBOinform_copyi); %
        end
        
        while ~isempty(atomPBC)
            %
            xyzdata=[];
            for i=1:size(trjdata_copy,1)
                if trjdata_copy(i,1)==atomPBC(1,1)
                    xyzdata(1,1)=atomPBC(1,1);
                    xyzdata(1,2)=atomPBC(1,2);
                    xyzdata(1,3)=trjdata_copy(i,3);%X
                    xyzdata(1,4)=trjdata_copy(i,4);%Y
                    xyzdata(1,5)=trjdata_copy(i,5);%Z
                end
            end
            if atomPBC(1,3)>0%
                for i=1:atomPBC(1,3)
                    for j=1:size(trjdata_copy,1)
                        if trjdata_copy(j,1)==atomPBC(1,i+3)
                            count4=size(xyzdata,1)+1;
                            xyzdata(count4,1)=atomPBC(1,i+3);
                            xyzdata(count4,2)=trjdata_copy(j,2);
                            xyzdata(count4,3)=trjdata_copy(j,3);%X
                            xyzdata(count4,4)=trjdata_copy(j,4);%Y
                            xyzdata(count4,5)=trjdata_copy(j,5);%Z
                        end
                    end
                end
            end
            %
            if size(xyzdata,1)>1
                xyzdata=BondForm(xyzdata,RadiiData,boxsize);
            end
            %             trjdata_temp(size(trjdata_temp,1)+1:size(trjdata_temp,1)+size(xyzdata,1),:)=xyzdata(:,:);%
            for i=1:size(xyzdata,1)
                if isempty(trjdata_temp)
                    trjdata_temp(size(trjdata_temp,1)+1,:)=xyzdata(i,:);
                else
                    if ~ismember(xyzdata(i,1),trjdata_temp(:,1))
                        trjdata_temp(size(trjdata_temp,1)+1,:)=xyzdata(i,:);
                    end
                    
                end
            end
            
            atomPBC(1,:)=[];
            if ~isempty(atomPBC)
                [atomPBC,atomPBCmember]=BondLink(atomPBC,atomPBCmember,tarBOinform_copyi); 
            end
        end
        %
        xmean=mean(trjdata_temp(:,3));
        if xmean<-(boxsize(1,2)-boxsize(1,1))/2
            trjdata_temp(:,3)=trjdata_temp(:,3)+(boxsize(1,2)-boxsize(1,1));
            fprintf('\nPBC transformation according to the X coordinate deviation');
        elseif xmean>(boxsize(1,2)-boxsize(1,1))*1.5
            trjdata_temp(:,3)=trjdata_temp(:,3)-(boxsize(1,2)-boxsize(1,1));
            fprintf('\nPBC transformation according to the X coordinate deviation');
        end
        ymean=mean(trjdata_temp(:,4));
        if ymean<-(boxsize(2,2)-boxsize(2,1))/2
            trjdata_temp(:,4)=trjdata_temp(:,4)+(boxsize(2,2)-boxsize(2,1));
            fprintf('\nPBC transformation according to the Y coordinate deviation');
        elseif ymean>(boxsize(2,2)-boxsize(2,1))*1.5
            trjdata_temp(:,4)=trjdata_temp(:,4)-(boxsize(2,2)-boxsize(2,1));
            fprintf('\nPBC transformation according to the Y coordinate deviation');
        end
        zmean=mean(trjdata_temp(:,5));
        if zmean<-(boxsize(3,2)-boxsize(3,1))/2
            trjdata_temp(:,5)=trjdata_temp(:,5)+(boxsize(3,2)-boxsize(3,1));
            fprintf('\nPBC transformation according to the Z coordinate deviation');
        elseif zmean>(boxsize(2,2)-boxsize(2,1))*1.5
            trjdata_temp(:,5)=trjdata_temp(:,5)-(boxsize(3,2)-boxsize(3,1));
            fprintf('\nPBC transformation according to the Z coordinate deviation');
        end
        for i=1:size(trjdata_temp,1)
            for j=1:size(trjdata,1)
                if trjdata_temp(i,1)==trjdata(j,1)
                    trjdata(j,3:5)=trjdata_temp(i,3:5);
                end
            end
        end
        
        fprintf('\nSuccessfully unwrap the %d th species, %d species need unwrap',count2,count1)
        count2=count2+1;%
    end
    fprintf('\nGood job! The coordnate of %d species in total is unwraped',size(species_table,1))

else
    error('Illegal parameters for BOXsize!')
end

if strcmpi(BOXsize,'y')
    trjdata(:,3)=(trjdata(:,3)-boxsize(1,1))/(boxsize(1,2)-boxsize(1,1));
    trjdata(:,4)=(trjdata(:,4)-boxsize(2,1))/(boxsize(2,2)-boxsize(2,1));
    trjdata(:,5)=(trjdata(:,5)-boxsize(3,1))/(boxsize(3,2)-boxsize(3,1));
end
end