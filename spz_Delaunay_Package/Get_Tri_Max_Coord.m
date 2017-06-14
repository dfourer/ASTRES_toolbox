function Max_Coord = Get_Tri_Max_Coord(triZ, This_Coord)
% Get_Tri_Max_Coord(triZ, This_Coord)
%   Returns max coordinate value out of the three vertices' coordinate for
%   each triangle element of triZ.
%
%   Max_Coord = Get_Tri_Max_Coord(triZ, This_Coord)
%       tri1Z       - Array of triangles
%       This_Coord  - Vector of vertices' coordinate
%       Max_Coord   - Vector of max coordinate values for each triangle
%
% Ph. Depalle 
% 2015, June, 10th

Taille = length(triZ);
Max_Coord = zeros(Taille, 1);

for k = 1:Taille
    Max_Coord(k) = max(This_Coord(triZ(k, :)));  
end
