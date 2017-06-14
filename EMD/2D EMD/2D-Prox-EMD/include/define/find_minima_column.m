function columnmax=find_minima_column(column)

[sizec,~]=size(column);
columnmax=zeros(sizec,1);
for c=2:sizec-1
    if (column(c)<column(c-1) & column(c)<column(c+1))
        columnmax(c)=1;
    end
end
end