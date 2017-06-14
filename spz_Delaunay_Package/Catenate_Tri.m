function [tri_out, x_out, y_out, i_min_out] = Catenate_Tri(tri1, x1, y1, tri2, x2, y2, min_index)
% Catenate_Tri(tri1, x1, y1, tri2, x2, y2, min_index)
%   Smart catenatation of two sets of triangles, as well as their corresponding vectors
%   of vertices' coordinates.
%
%   [tri_out, x_out, y_out, i_min_out] = Catenate_Tri(tri1, x1, y1, tri2, x2, y2, min_index)
%       tri1      - Current set of triangles  
%       x1        - Vector of vertices' coordinate x for current set of triangles
%       y1        - Vector of vertices' coordinate y for current set of triangles
%       tri2      - New set of triangles
%       x2        - Vector of vertices' coordinate x for new set of triangles
%       y2        - Vector of vertices' coordinate y for new set of triangles
%       min_index - Index for starting exploration of current vector of vert.
%       tri_out   - Resulting set of triangles
%       x_out     - Vector of vertices' coordinate x for resulting set of triangles
%       y_out     - Vector of vertices' coordinate y for resulting set of triangles
%       i_min_out - Index for starting exploration of resulting vector of vert.
    

Taille_in1  = length(tri1);
Taille_in2  = length(tri2);

tri_out     = [tri1; tri2];     % Assumes there is no duplicate in the two sets
x_out       = x1;
y_out       = y1;
max_index   = length(x_out);    % Useful when adding a new vertex
statut      = zeros(size(x2));  % Store vertex index pointing to updated 
                                % coordinates or equals 0 when vertex never 
                                % met before

if(max_index == 0)
    i_min_out = 1;
else
    i_min_out = max_index;
end
                              
for tri_index = Taille_in1+1:Taille_in1+Taille_in2
    for this_tri_sommet = 1:3
       this_sommet = tri_out(tri_index, this_tri_sommet);
       % Already identified in the resulting vertices' set ?
       if(statut(this_sommet) ~= 0)
           tri_out(tri_index, this_tri_sommet) = statut(this_sommet);
       else
           is_new = 0; % Not identified, but might be there
           for sommet_index = min_index:max_index
               if (x2(this_sommet) == x_out(sommet_index)...
                   &&...
                   y2(this_sommet) == y_out(sommet_index))
                       is_new = 1;
                       found_sommet = sommet_index;
                end
           end
           if(is_new == 1) % Vertex already in the resulting vector
               tri_out(tri_index, this_tri_sommet) = found_sommet;
           else % New vertex. Add it to the resulting vector
               x_out = [x_out; x2(this_sommet)];
               y_out = [y_out; y2(this_sommet)]; 
               tri_out(tri_index, this_tri_sommet) = max_index+1; 
               max_index = max_index+1;
           end
       end      % Is this vertex already registered in the resulting set
    end         % and for each of its vertices 
end             % For each triangle of the new set
