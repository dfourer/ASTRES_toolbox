function spz_delaunay_domains_image_save(this_Image, N1, N2, File_Name)
%   spz_delaunay_domains_image_save(this_Image, N1, N2, File_Name)
%   Save a time-frequency image into a file - Stored as binary 16-bit integers.
%   N1 and N2 are saved first into the file, prior to the image content.
%
%   spz_delaunay_domains_image_save(this_Image, N1, N2, File_Name)
%       this_Image  - N1xN2 Matrix
%       N1          - Size 1 (Vertical, e.g. Frequency)
%       N2          - Size 2 (Horizontal, e.g. Time)
%       File_Name   - Name of the file where the image is to be stored.
%
% Ph. Depalle 
% 2015, June, 14th


File_Descriptor = fopen(File_Name, 'wb');

fwrite(File_Descriptor, N1, 'int16');
fwrite(File_Descriptor, N2, 'int16');

fwrite(File_Descriptor, this_Image, 'int16');

fclose(File_Descriptor);
end