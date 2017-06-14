function Read_Signal = spz_delaunay_signal_read(File_Name)

File_Descriptor = fopen(File_Name, 'rb');
Read_Signal     = fread(File_Descriptor, 'double');
fclose(File_Descriptor);

end