function update_param(hObject,handles)
% update_param : initialize the parameters for ridge extraction. Called
% by popupmenu8 callback.

switch(get(hObject,'Value'))
    case 1 % Brevdo
        set(handles.edit4,'String','1'); % nr
        set(handles.edit5,'String','0.01'); % lambda
        set(handles.edit6,'String','2'); % clwin
        set(handles.text3,'String','Num Modes');
        set(handles.text4,'String','lambda');
        set(handles.text5,'String','clwin');
    case 2 % Brevdo
        set(handles.edit4,'String','1'); % nr
        set(handles.edit5,'String','0.01'); % lambda
        set(handles.edit6,'String','2'); % clwin
        set(handles.text3,'String','Num Modes');
        set(handles.text4,'String','lambda');
        set(handles.text5,'String','clwin');
    case 3 % Meignen
        set(handles.edit4,'String','auto'); % nr
        set(handles.edit5,'String','auto'); 
        set(handles.edit6,'String','auto');
        set(handles.text3,'String','Num Modes');
        set(handles.text4,'String','rrec');
        set(handles.text5,'String','gamma min');
end
