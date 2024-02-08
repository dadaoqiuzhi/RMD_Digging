%This function is used to transform a number to decimal system.
function num_10base = Deconvolu_SysConvert(given_num_char,base)
num_table = {'0','1','2','3','4','5','6','7','8','9',...
    'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z',...
    'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};
[~,matches] = strsplit(given_num_char,num_table,'CollapseDelimiters',false);
num_10base = 0;
for i = 1:length(matches)
    [~,num_true] = ismember(matches{i},num_table);
    num_10base = num_10base + (num_true-1)*base^(length(matches)-i);
end
end