function data = num_obtain(datasplit,coloum_num)
data = [];
num_length_maybe = length(datasplit);
num_length_true = 0;

for i = 1:size(datasplit,2)
    TF = isstrprop(datasplit{i},'digit'); 
    if ismember('.',datasplit{i})
        dot_num = 1;
    else
        dot_num = 0;
    end
    if ismember('-',datasplit{i})
        minus_num = 1;
    else
        minus_num = 0;
    end
    if sum(TF)+dot_num+minus_num == length(datasplit{i})
        num_length_true = num_length_true + 1;
    end
end

if num_length_maybe == num_length_true 
    if num_length_true == coloum_num
        for i = 1:size(datasplit,2)
            data(length(data)+1) = str2num(datasplit{i});
        end
    end
else
    data='not number line';
end

end