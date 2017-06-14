function sqplot(Tx,dt,fs,hpar)
% sqplot : displays a synchrosqueezed or wavelet transform.

cla(hpar);

if isempty(Tx) || isempty(dt) || isempty(fs)
    warndlg('This data has not been computed yet');
    return;
end

N = size(Tx,2);
h = imagesc(dt*(0:N-1),fs,log(1+abs(Tx)),'Parent',hpar);
set(hpar,'YDir','normal');
set(hpar,'YTick',(linspace(fs(1),fs(end),10)),'YTickLabel',fs(floor(linspace(1,length(fs),10))));

end

