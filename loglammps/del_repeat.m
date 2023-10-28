function rawdata = del_repeat(rawdata,data)
if ismember(data(1),rawdata(:,1))
    for i = 1:size(rawdata,1)
        if data(1) == rawdata(i,1)
            break
        end
    end
    rawdata(i:size(rawdata,1),:)=[];
end
end