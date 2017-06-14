function linemax=find_minima_line(line)

[~,sizel]=size(line);
linemax=zeros(1,sizel);
for l=2:sizel-1
    if (line(l)<line(l-1) & line(l)<line(l+1))
        linemax(l)=1;
    end
end
end