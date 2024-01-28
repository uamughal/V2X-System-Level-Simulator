function [ CoMP_sites ] = initialize_CoMP_sites(LTE_config, eNodeBs )
% (c) Thomas Blazek, 2014, thomas.blazek@nt.tuwien.ac.at
% This script initializes the CoMP sites. Each site has eNodeBs assigned. 
% Right now, this is just a showcase. Possible configurations are: 
% 'trivial': Each eNodeB is a seperate site, no coordination
% 'global': All eNodeBs are connected, cooperation across all eNodeBs is
% possible

% One has to differentiate between the CoMP scheduler and the LTE
% scheduler. 

  CoMP_sites = CoMP.CoMP_site(LTE_config, eNodeBs, 1);
 
  CoMP_sites.scheduler = CoMP.schedulers.roundRobinCoMP(CoMP_sites);
  CoMP_sites = [];
end

