function spz_delaunay_save(Triangles, x_Coord, y_Coord, Domains, Block_Size, Root_File_Name)

% Calls:
%   spz_delaunay_domains_save
% 	spz_delaunay_triangle_coordinate_save
% 	spz_delaunay_triangle_save
    
    
Tri_File_Name = sprintf('%s%s', Root_File_Name, '_tri.spz');
spz_delaunay_triangle_save(Triangles, Block_Size, Tri_File_Name);

x_Coord_File_Name = sprintf('%s%s', Root_File_Name, '_x.spz');
spz_delaunay_triangle_coordinate_save(x_Coord, x_Coord_File_Name);

y_Coord_File_Name = sprintf('%s%s', Root_File_Name, '_y.spz');
spz_delaunay_triangle_coordinate_save(y_Coord, y_Coord_File_Name);

if(max(size(Domains)) ~= 0)
    Domains_File_Name = sprintf('%s%s', Root_File_Name, '_dom.spz');
    spz_delaunay_domains_save(Domains, Domains_File_Name);
else
    disp('Pas de domaine');
end

end
 

