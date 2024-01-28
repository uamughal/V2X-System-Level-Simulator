function [ UEs, eNodeBs_, eNodeBs_sectors_ , networkPathlossMap] = add_Vehicle_eNodeBs( LTE_config,UEs,sites,eNodeBs,networkPathlossMap)
%ADD_VEHICLE_ENODEBS Summary of this function goes here
%   Detailed explanation goes here

for u_=1:length(UEs)   
    LTE_V_pos(u_,:)= UEs(u_).pos;
    % match the pos with the map (resolusion)
    roi_min=networkPathlossMap.coordinate_origin;
    data_res=networkPathlossMap.data_res;
    [pos_pixel,pos_pixel_exact] = LTE_common_pos_to_pixel(LTE_V_pos, roi_min, data_res);
    LTE_V_pos=LTE_common_pixel_to_pos(pos_pixel, roi_min, data_res);
end

s_idx = 1+length(eNodeBs);
for r_=1:size(LTE_V_pos,1)
    % Create the LTE-V sites
    UEs(r_).sites     = network_elements.eNodeB;
    LTE_V_sites(r_)           = network_elements.eNodeB;
    LTE_V_sites(r_).id        = length(sites)+r_;   
    LTE_V_sites(r_).site_type = 'LTE-V';   %need modification
   
    % Create the LTE-V sectors
    LTE_V_sites(r_).sectors              =network_elements.eNodeB_sector;
    UEs(r_).sectors                      =network_elements.eNodeB_sector;
    for s_=1:length(LTE_config.LTE_V_config.sector_azimuths)
    UEs(r_).sectors(s_)                   =network_elements.eNodeB_sector;
    LTE_V_sites(r_).sectors(s_)           = network_elements.eNodeB_sector;
    LTE_V_sites(r_).sectors(s_).parent_eNodeB = LTE_V_sites(r_);
    LTE_V_sites(r_).sectors(s_).id            = s_;
    
    if strcmp(LTE_config.LTE_V_config.eNB_distribution,'four_side_of_lane')
    LTE_V_sites(r_).sectors(s_).azimuth       =  utils.miscUtils.wrapTo359(LTE_config.antenna_azimuth_offsett + LTE_config.LTE_V_config.sector_azimuths(s_)); % 3sector antenna
    LTE_V_sites(r_).pos       = LTE_V_pos(r_,:);  
    end
    
    
    LTE_V_sites(r_).sectors(s_).max_power     = LTE_config.LTE_V_config.tx_power_W;
    LTE_V_sites(r_).sectors(s_).antenna_type  = LTE_config.antenna.antenna_gain_pattern;
    LTE_V_sites(r_).sectors(s_).nTX           = LTE_config.LTE_V_config.nTX;
    LTE_V_sites(r_).sectors(s_).eNodeB_id     = s_idx;
    LTE_V_sites(r_).sectors(s_).tx_height      =LTE_config.LTE_V_config.tx_height;
    LTE_V_sectors(s_idx-length(eNodeBs)) =  LTE_V_sites(r_).sectors(s_);
    % Attach the correct antenna to the eNodeB
    antennas.antenna.attach_antenna_to_eNodeB(LTE_V_sites(r_).sectors(s_),LTE_config.LTE_V_config);
         
    % Create the macroscopic pahloss model that will be used. Use the same one as in the macro sites (may change in the future)
    if LTE_config.debug_level>=1
        fprintf('LTE-V Site %d: eNodeB %d ',length(sites)+r_,s_);
    end
    %LTE_V eNodeB sectors
   LTE_V_sectors(s_idx-length(eNodeBs)).macroscopic_pathloss_model= macroscopic_pathloss_models.generalPathlossModel.generateMacroscopicPathlossModel(...
            LTE_config,...
            LTE_config.LTE_V_config.macroscopic_pathloss_model,...
            LTE_config.frequency,...
            LTE_config.LTE_V_config.macroscopic_pathloss_model_settings);  
        
    UEs(r_).sectors(s_)= LTE_V_sites(r_).sectors(s_);
     s_idx=1+s_idx;
    end
    UEs(r_).sites=LTE_V_sites(r_);
end


if size(LTE_V_pos,1)>0
    % There are LTE-V cells
    eNodeBs_         = [ sites LTE_V_sites ];
    eNodeBs_sectors_ = [ eNodeBs  LTE_V_sectors ];
else
    % No LTE-V cells. Return the same values
    eNodeBs_           = sites;
    eNodeBs_sectors_    = eNodeBs;
    return
end

%% Calculate final pathloss. With and without minimum coupling loss
%{
if LTE_config.debug_level>=1
    fprintf('Creating LTE-V eNodB pathloss map\n');
end

% Preallocate to make it faster for the following function (but not 100% necessary)
    networkPathlossMap.pathloss(:,:,(length(eNodeBs)+1):(length(eNodeBs)+length(LTE_V_sites))) = 0;
    macroscopic_pathloss_models.generalPathlossModel.calculate_pathloss_maps(LTE_config,LTE_V_sites,networkPathlossMap);

%% Apply Minimum Coupling Loss (MCL). After the antenna gain, as TS.942-900
%  states: RX_PWR = TX_PWR ¨C Max (pathloss ¨C G_TX ¨C G_RX, MCL)
%  KNOWN ISSUE: this assumes that G_RX is 0dB, which will normally be the
%  case for a mobile terminal. This would have to be moved to the link
%  level model (UE) if G_RX is to be taken into account

% if LTE_config.debug_level>=1
%    fprintf('Applying Minimum Coupling Loss of %d dB to the LTE-V eNodeBs\n',LTE_config.LTE_V_config.minimum_coupling_loss);
% end
% networkPathlossMap.apply_MCL(LTE_config.LTE_V_config.minimum_coupling_loss);

%% Fill in pathloss data in the pathlossMap
all_pathloss_models      = [LTE_V_sectors.macroscopic_pathloss_model];
all_pathloss_model_names = {all_pathloss_models.name};

networkPathlossMap.name = {networkPathlossMap.name{1:end} all_pathloss_model_names{:}};

%}
end

