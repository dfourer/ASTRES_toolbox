function spz_delaunay_triangle_save(Triangles, Block_size, File_Name)


[Nb_Triangles, Nb_Vertices] = size(Triangles);

File_Descriptor = fopen(File_Name, 'wb');

% Write a header with descriptive values
fwrite(File_Descriptor, Nb_Triangles, 'int32');
fwrite(File_Descriptor, Nb_Vertices, 'int32');
fwrite(File_Descriptor, Block_size, 'int32');

for this_Vertex = 1:Nb_Vertices
    fwrite(File_Descriptor, Triangles(:, this_Vertex), 'int32');
end

fclose(File_Descriptor);

end