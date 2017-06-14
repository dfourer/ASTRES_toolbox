function spz_delaunay_domains_save(Domains, File_Name)

Nb_Domains = length(Domains);
for this_Domain = 1:Nb_Domains
    Length_Domain(this_Domain) = length(Domains{this_Domain});
end

File_Descriptor = fopen(File_Name, 'wb');

%Write a header with descriptive values
fwrite(File_Descriptor, Nb_Domains, 'int32');
fwrite(File_Descriptor, Length_Domain, 'int32');

for this_Domain = 1:Nb_Domains
    fwrite(File_Descriptor, Domains{this_Domain}, 'int32');
end

fclose(File_Descriptor);
end