function triZ_out = Select_Triangles(triZ_in, xZ_in, yZ_in, x_min, x_max, y_min, y_max, first_time)
% triZ_out = Select_Triangles(triZ_in, xZ_in, yZ_in, x_min, x_max, y_min, y_max, first_time)
%   Returns a subset of triangles, which vertices coordinates values lie
%   within a rectangle defined by its vertices
%
%   triZ_out = Select_Triangles(triZ_in, xZ_in, yZ_in, x_min, x_max, y_min, y_max, first_time)
%       triZ_in     - Input Array of triangles to be selected
%       xZ_in       - Vector of X vertices' coordinates
%       yZ_in       - Vector of Y vertices' coordinates
%       x_min       - lower limit on X vertices' coordinates
%       x_max       - Upper limit on X vertices' coordinates
%       y_min       - lower limit on Y vertices' coordinates
%       y_max       - Upper limit on Y vertices' coordinates
%       first_time  - Use 0 instead of y_min for the 2nd condition

% Calls
%	Get_Tri_Max_Coord.m
%	Get_Tri_Min_Coord.m
%
% Ph. Depalle 
% 2015, June, 10th

Min_XZl           = Get_Tri_Min_Coord(triZ_in, xZ_in);
Max_XZl           = Get_Tri_Max_Coord(triZ_in, xZ_in);

Min_YZl           = Get_Tri_Min_Coord(triZ_in, yZ_in);
Max_YZl           = Get_Tri_Max_Coord(triZ_in, yZ_in);

if(first_time == true)
    y_min_2 = 0;
else
    y_min_2 = y_min;
end

Selected_Indices  = find(  (Max_YZl >= y_min) & (Min_YZl <= y_max)...
                         & (Min_YZl >= y_min_2)...
                         & (Min_XZl >= x_min) & (Max_XZl <= x_max));
triZ_out          = triZ_in(Selected_Indices, :);

end