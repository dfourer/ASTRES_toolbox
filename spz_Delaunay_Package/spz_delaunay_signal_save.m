function spz_delaunay_signal_save(this_Signal, File_Name)

File_Descriptor = fopen(File_Name, 'wb');

fwrite(File_Descriptor, this_Signal, 'double');

fclose(File_Descriptor);

end