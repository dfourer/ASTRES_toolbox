function changepopupv(hobj,hfig,haxes,hsli)
% change the graphics axes3 or axes4 of figure segtool.
%   Internal function.

switch(get(hobj,'Value'))
    case 1 % signal
        s = getappdata(hfig,'s');
        dt = getappdata(hfig,'dt');
        if ~isempty(s) && ~isempty(dt)
            cla(haxes);
            plot(dt*(0:length(s)-1),s,'Parent',haxes);
            set(hsli,'Visible','off');
        else
            warndlg('No signal loaded');
        end
    case 2 % Wx section
        Wx = getappdata(hfig,'Wx');
        as = getappdata(hfig,'as');
        if isempty(Wx) || isempty(as)
            warndlg('No Wavelet transform to plot');
            return
        end
        set(hsli,'Visible','on','Value',0.5);
        b = floor(size(Wx,2)/2);
        cla(haxes);
        plot(log(as)/log(2),abs(Wx(:,b)),'Parent',haxes);
    case 3 % IMFs
        imf = getappdata(hfig,'imf');
        dt = getappdata(hfig,'dt');
        if isempty(imf) || isempty(dt)
            warndlg('No IMFs to plot');
            return
        end
        set(hsli,'Visible','off');
        cla(haxes);
        imfplot(imf,dt,haxes);
    otherwise
        imf = getappdata(hfig,'imf');
        dt = getappdata(hfig,'dt');
        if isempty(imf) || isempty(dt)
            warndlg('No IMFs to plot');
            return
        end
        set(hsli,'Visible','off');
        [nr N] = size(imf);
        if get(hobj,'Value')>nr+3
            warndlg('Bad nummer of IMF');
        end
        cla(haxes);
        plot(dt*(0:N-1),imf(get(hobj,'Value')-3,:),'Parent',haxes);
        legend(['IMF ' num2str(get(hobj,'Value')-3)]);
end
        