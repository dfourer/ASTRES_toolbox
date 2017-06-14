function Min_Coord = Get_Tri_Min_Coord(triZ, This_Coord)
% Get_Tri_Min_Coord(triZ, This_Coord)
%   Returns min coordinate value out of the three vertices' coordinate for
%   each triangle element of triZ.
%
%   Min_Coord = Get_Tri_Min_Coord(triZ, This_Coord)
%       tri1Z       - Array of triangles
%       This_Coord  - Vector of vertices' coordinate
%       Min_Coord   - Vector of min coordinate values for each triangle
%
%
% Ph. Depalle 
% 2015, June, 10th

Taille = length(triZ);
Min_Coord = zeros(Taille, 1);

for k = 1:Taille
    Min_Coord(k) = min(This_Coord(triZ(k, :)));  
end
