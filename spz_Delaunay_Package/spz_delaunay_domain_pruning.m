function [new_Domains] = spz_delaunay_domain_pruning(Domains, Threshold)

N_Domain = length(Domains);

new_Domains_index = 1;
new_Domains = {};

for this_domain = 1:N_Domain
    if(length(Domains{this_domain}) > Threshold)
        length(Domains{this_domain});
        new_Domains{new_Domains_index} = Domains{this_domain};
        new_Domains_index = new_Domains_index+1;
    end
end
