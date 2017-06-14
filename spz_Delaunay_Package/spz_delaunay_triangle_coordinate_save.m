function spz_delaunay_triangle_coordinate_save(this_Coord, File_Name)

File_Descriptor = fopen(File_Name, 'wb');

fwrite(File_Descriptor, this_Coord, 'int32');

fclose(File_Descriptor);

end