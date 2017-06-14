function upplot(hf,ha,val)

switch(val)
    case 1 % clear graphic
        cla(ha);
    case 2 % Wx
        dt = getappdata(hf,'dt');
        Wx = getappdata(hf,'Wx');
        as = getappdata(hf,'as');
        sqplot(Wx,dt,as,ha);
    case 3 % Tx
        dt = getappdata(hf,'dt');
        Tx = getappdata(hf,'Tx');
        fs = getappdata(hf,'fs');
        sqplot(Tx,dt,fs,ha);
    case 4 % Tx and ridge extraction
        dt = getappdata(hf,'dt');
        Tx = getappdata(hf,'Tx');
        fs = getappdata(hf,'fs');
        Cs = getappdata(hf,'Cs');
        clwin = getappdata(hf,'clwin');
        sqplotRE(Tx,dt,fs,Cs,clwin,ha);
end
