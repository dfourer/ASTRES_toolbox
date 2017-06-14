function slidercall(hobj,hfig,haxes)
% slidercall : internal function for the GUI

Wx = getappdata(hfig,'Wx');
as = getappdata(hfig,'as');
if isempty(Wx) || isempty(as)
	warndlg('No Wavelet transform to plot');
    return
end
val = get(hobj,'Value');
N = size(Wx,2);
b = floor(N*val+1-0.001);
if b<1
    b=1;
elseif b>N
    b=N;
end
cla(haxes);
plot(log(as)/log(2),abs(Wx(:,b)),'Parent',haxes);
     