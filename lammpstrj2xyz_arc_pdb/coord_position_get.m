function trjdata=coord_position_get(coord_position,datasplit)
trjdata=[];
for i=1:5
    trjdata(i)=str2num(datasplit{coord_position(i)});
end
end