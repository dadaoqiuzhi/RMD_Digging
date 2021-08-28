%scrit file name speciesbond_classify
%purpose:
%This program is used to rearrange the molecule and corresponding bond boder information processed and generated
%the by bondorder_deepmining program
%version 1;2018.6.29
disp('Welcome!--by Qiang Liu @Institute of Nuclear Physics and Chemistry, China Academy of Engineering Physics; Email: liubinqiang@163.com');
disp('Repository adress of the Source code on github: https://github.com/dadaoqiuzhi/RMD_Digging');

fprintf('\nThis program takes the following preferential options into consideration: elements (C,H,O, etc), total atom number (atsum). multi-step refinement can be performed\n')
elementsort=input('\nPlease input the ranking priority, seperated by white space\n','s');
elementsort=upper(elementsort);
fprintf('\nspeciesbond_classify is running, please wait...\n')
tarBOinformcopy=tarBOinform;tarelenummatchcopy=tarelenummatch;
[tarelenumrow,tarelenumcol]=size(tarelenummatch);
speciescol=tarelenumcol/2;
for i=1:speciescol
    block={};block(:,:)=tarelenummatch(:,2*i-1:2*i);
    atomnummolecule=atomnummolecule_strcat(block);
    tarelenummatch{tarelenumrow+1,2*i-1}='#';
    tarelenummatch{tarelenumrow+2,2*i-1}=atomnummolecule;
    tarelenummatch{tarelenumrow+1,2*i}='#';
    tarelenummatch{tarelenumrow+2,2*i}=i;
end

eleidmat=[];
for i=1:speciescol
    for j=1:tarelenumrow
    eleidmat(i,j)=tarelenummatch{j,2*i};
    end
    eleidmat(i,j+1)=0;
    for k=1:tarelenumrow
    eleidmat(i,j+1)=eleidmat(i,j+1)+eleidmat(i,k);
    end
    eleidmat(i,j+2)=tarelenummatch{j+2,2*i};
end

elementsort=strtrim(elementsort);
elementsort=strsplit(elementsort);
elementstd=tarelenummatch(1:tarelenumrow,1)';
elementstd{1,tarelenumrow+1}='atsum';
stdsort=[];sortrows(matrixA,[3,2,1])
for i=1:length(elementsort)
    loglocat=ismember(elementstd,elementsort{i});
    if sum(loglocat(1:tarelenumrow))==0
        disp('the input element has no legal molecular formula or is not assigned to sort, please check it')
    else
        index=find(loglocat==1);
        stdsort(i)=index;
    end
end
eleidmat=sortrows(eleidmat,stdsort);

dataclass={};
for i=1:speciescol
    dataclass{i,1}=tarelenummatch{tarelenumrow+2,2*eleidmat(i,5)-1};
    dataclass{i,2}=eleidmat(i,5);
end

%resort tarBOinform block
[tarBOrow,tarBOcol]=size(tarBOinform);
seperatorrow=[];k=1;
for i=1:tarBOrow
    if tarBOinform{i,1}=='#'
      seperatorrow(k)=i;
      k=k+1;
    end
end

moleculebond=[];
if length(seperatorrow)==1
    moleculebond=[1,1,seperatorrow(1)];
else
    moleculebond(1,:)=[1,1,seperatorrow(1)];
    for i=2:speciescol
        moleculebond(i,1)=i;
        moleculebond(i,2)=seperatorrow(i-1)+1;
        moleculebond(i,3)=seperatorrow(i);
    end
end

tarBOinformnew={}; 
for i=1:speciescol
    for j=1:speciescol
        if moleculebond(j,1)==dataclass{i,2}
            rownum=moleculebond(j,3)-moleculebond(j,2);
            [tarrow,tarcol]=size(tarBOinformnew);
            tarBOinformnew(tarrow+1:tarrow+rownum+1,:)=tarBOinform(moleculebond(j,2):moleculebond(j,3),:);
            dataclass{i,3}=tarrow+1;dataclass{i,4}=tarrow+rownum+1;
        end
    end
end
tarBOinform=tarBOinformcopy;tarelenummatch=tarelenummatchcopy;

tarmoleculebond={};
for i=1:speciescol
    [tarmoborow,tarmobocol]=size(tarmoleculebond);
    tarmoleculebond{tarmoborow+1,1}=dataclass{i,1};
    tarmoleculebond{tarmoborow+1,2}=dataclass{i,2};
    for j=3:tarBOcol
        tarmoleculebond{tarmoborow+1,j}=[];
    end
    rownum=dataclass{i,4}-dataclass{i,3}+1;
    tarmoleculebond(tarmoborow+2:tarmoborow+rownum+1,:)=tarBOinformnew(dataclass{i,3}:dataclass{i,4},:);
end
disp('speciesbond_classify is successfully finished')
fprintf('\nResorted molecular formula, previous number, BO in tarBOinformnew is saved in dataclass\n')
fprintf('Resorted BO corresponding to dataclass is saved in tarBOinformnew\n')
fprintf('Resorted molecular formula and BO informaton is saved in tarmoleculebond\n')

clear atomnummolecule block eleidmat elementsort i index j k loglocat moleculebond rownum seperatorrow speciescol
clear stdsort tarBOcol tarBOinformcopy tarBOrow tarcol tarelenumcol tarelenummatchcopy tarelenumrow
clear tarrow tarmoborow tarmobocol elementstd

