function [eNodeBs_sites,eNodeBs,networkMacroscopicPathlossMap] = generated_network(LTE_config)
% Generate an hexagonal network with a free space pathloss model
% (c) Josep Colom Ikuno, INTHFT, 2008
% www.nt.tuwien.ac.at
% output:   eNodeBs            ... contains info reagarding the BTSs and its
%                                  sectors
%           pathloss_data      ... [heightxwidthx3xnBTS] double
%                                  Pathloss data for each sector (including
%                                  antenna gain).
%                                  [y,x,sector_num,brts_num]

%% Necessary data loaded in the LTE_init_config config
data_res               = LTE_config.map_resolution;           % Meters/pixel, resolution of the map
eNodeB_sector_tx_power = LTE_config.eNodeB_tx_power;          % eNodeB tx power (Watts/sector)

%% Some initialisation params
ROI_increase_factor = 0.1;

%% Create the eNodeBs
if LTE_config.debug_level>=1
    fprintf(['Creating eNodeBs, (geometry: "' LTE_config.network_geometry '")\n']);
end

if ~isfield(LTE_config,'network_geometry') % in case an older config file is used - hexagonal grid by default
   LTE_config.network_geometry = 'regular_hexagonal_grid'; 
end

switch LTE_config.network_geometry
    case 'regular_hexagonal_grid'
        eNodeBs_sites = LTE_init_create_hexagonal_eNodeB_grid(LTE_config);
        inter_eNodeB_distance = LTE_config.inter_eNodeB_distance;
end

% Add the Antennas to the eNodeBs
s_idx   = 1;
eNodeBs = network_elements.eNodeB_sector; % Initialization
for b_ = 1:length(eNodeBs_sites)
    % Create the eNodeB_sector objects
    % Writing eNodeBs(b_).sectors(1) gave me an error. Maybe a bug??
    eNodeBs_sites(b_).sectors    = network_elements.eNodeB_sector;
    for s_ = 1:length(LTE_config.sector_azimuths)
        eNodeBs_sites(b_).sectors(s_)               = network_elements.eNodeB_sector;
        eNodeBs_sites(b_).sectors(s_).parent_eNodeB = eNodeBs_sites(b_);
        eNodeBs_sites(b_).sectors(s_).id            = s_;
        eNodeBs_sites(b_).sectors(s_).azimuth       = utils.miscUtils.wrapTo359(LTE_config.antenna_azimuth_offsett + LTE_config.sector_azimuths(s_));
        eNodeBs_sites(b_).sectors(s_).max_power     = eNodeB_sector_tx_power;
        eNodeBs_sites(b_).sectors(s_).antenna_type  = LTE_config.antenna.antenna_gain_pattern;
        eNodeBs_sites(b_).sectors(s_).nTX           = LTE_config.nTX;
        eNodeBs_sites(b_).sectors(s_).tx_height     = LTE_config.tx_height;
        
        eNodeBs_sites(b_).sectors(s_).eNodeB_id     = s_idx;
        eNodeBs(s_idx)                              = eNodeBs_sites(b_).sectors(s_);
        
        % Attach the correct antenna to the eNodeB
        antennas.antenna.attach_antenna_to_eNodeB(eNodeBs_sites(b_).sectors(s_),LTE_config);
        
        % Create the macroscopic pahloss model that will be used
        if LTE_config.debug_level>=1
            fprintf('Site %d, eNodeB %d: ',b_,s_);
        end
        
        eNodeBs_sites(b_).sectors(s_).macroscopic_pathloss_model = macroscopic_pathloss_models.generalPathlossModel.generateMacroscopicPathlossModel(...
            LTE_config,...
            LTE_config.macroscopic_pathloss_model,...
            LTE_config.frequency1,...
            LTE_config.macroscopic_pathloss_model_settings);

        s_idx = s_idx + 1;
    end
end

%% Calculate ROI
number_of_eNodeB_sites = length(eNodeBs_sites);  % number of eNodeB sites
if (number_of_eNodeB_sites > 1)
    tx_pos = reshape([eNodeBs_sites.pos],2,[])';
    % Calculate ROI border points in ABSOLUTE coordinates
    roi_x = [min(tx_pos(:,1)),max(tx_pos(:,1))];
    roi_y = [min(tx_pos(:,2)),max(tx_pos(:,2))];
else % set the ROI equal to eNodeB distance if just one base station is used
    roi_x = [-inter_eNodeB_distance,inter_eNodeB_distance];
    roi_y = [-inter_eNodeB_distance,inter_eNodeB_distance];
end


%% Define an area of the ROI to map
% roi_reduction_factor times smaller and draw it. ABSOLUTE COORDINATES
roi_x = roi_x + ROI_increase_factor*abs(roi_x(2)-roi_x(1))*[-3,3];
roi_y = roi_y + ROI_increase_factor*abs(roi_y(2)-roi_y(1))*[-5,5];

%% Create pathlossMap
networkMacroscopicPathlossMap                        = channel_gain_wrappers.macroscopicPathlossMap;
networkMacroscopicPathlossMap.data_res               = data_res;
networkMacroscopicPathlossMap.roi_x                  = roi_x;
networkMacroscopicPathlossMap.roi_y                  = roi_y;

%% Calculate final pathloss. With and without minimum coupling loss
if LTE_config.debug_level>=1
    fprintf('Creating cell pathloss map\n');
end
macroscopic_pathloss_models.generalPathlossModel.calculate_pathloss_maps(LTE_config,eNodeBs_sites,networkMacroscopicPathlossMap);

%% Apply Minimum Coupling Loss (MCL). After the antenna gain, as TS.942-900
%  states: RX_PWR = TX_PWR – Max (pathloss – G_TX – G_RX, MCL)
%  KNOWN ISSUE: this assumes that G_RX is 0dB, which will normally be the
%  case for a mobile terminal. This would have to be moved to the link
%  level model (UE) if G_RX is to be taken into account


 if LTE_config.debug_level>=1
     fprintf('Applying Minimum Coupling Loss of %d dB\n',LTE_config.minimum_coupling_loss);
 end
 networkMacroscopicPathlossMap.apply_MCL(LTE_config.minimum_coupling_loss);

all_sectors = [eNodeBs_sites.sectors];
all_pathloss_models = [all_sectors.macroscopic_pathloss_model];
all_pathloss_model_names = {all_pathloss_models.name};
networkMacroscopicPathlossMap.name = all_pathloss_model_names;

