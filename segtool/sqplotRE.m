function sqplotRE(Tx,dt,fs,Cs,clwin,hpar)
% sqplotRE displays the SQ transqform Tx together with the ridges contained in
%   Cs. Internal function

if isempty(Tx) || isempty(dt) || isempty(fs) || isempty(Cs)
    warndlg('This data has not been computed yet');
    return;
end

N = size(Tx,2);
na = size(Tx,1);
nr = size(Cs,1);

cdat = Tx;
vmax = max(cdat(:))*1.5;
vmin = min(cdat(:))*0.8;
    
for b=1:N
    for j=1:nr
        %cdat(max(1,Cs(j,b)-clwin):min(na,Cs(j,b)+clwin),b) = vmin;
        cdat(Cs(j,b),b) = vmax;
    end
end

imagesc(dt*(0:N-1),(fs),log(1+abs(cdat)),'Parent',hpar);
set(hpar,'YDir','normal');
set(hpar,'YTick',(linspace(fs(1),fs(end),10)),'YTickLabel',2*fs(floor(1:(length(fs)/10):end)));

end

