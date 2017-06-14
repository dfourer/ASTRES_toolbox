function Shuffle_Figures(Delta_X, Delta_Y, Location_Reset)
%   Figure_Shuflle(Delta_X, Delta_Y, Location_Reset)
%       Delta_X     = Horizontal shift between figures (Left to Right)
%       Delta_Y     = Vertical shift between figures (Bottom Up)
%       Location_Reset  = true -> Figures clustered at the middle up of the screen
%                         false -> No Location shift
%
%
%   Ph. Depalle 
%   2015, July 2nd


List_of_Figures = get(0, 'Children');
Screen_Size     = get(0, 'screensize');
Screen_Length   = Screen_Size(3);
Screen_Height   = Screen_Size(4);

if(Location_Reset)
    Initial_Position    = List_of_Figures(end).Position;
    Initial_Position(1) = Screen_Length/2;
    Initial_Position(2) = Screen_Height - Initial_Position(4) - 96;
    fprintf('Reset Position');
else
    Initial_Position = List_of_Figures(end).Position;
    fprintf('Current Position');
end


for Figure_Index = length(List_of_Figures):-1:1
    this_Figure = List_of_Figures(Figure_Index);
    [this_Figure.Position] = Initial_Position;
    Initial_Position = Initial_Position + [Delta_X -Delta_Y 0 0];
    figure(Figure_Index);
end
end
