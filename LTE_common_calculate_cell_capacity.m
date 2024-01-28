function [...
    sector_capacity,...
    max_SINR_dB_all,...
    sector_assignment,...
    maxSINR_assignment,...
    diff_SINR_dB_all,...
    cell_sizes,...
    cell_centers...
    ] = LTE_common_calculate_cell_capacity(LTE_config,networkPathlossMap,eNodeBs,eNodeBs_sectors,varargin)
% Calculates average cell capacity based on the system settings and returns
% the cdf of the SINR for the target sector and cell (the same, just that smoother)
if ~isempty(varargin)
    shadow_fading_used = true;
    networkShadowFadingMap = varargin{1};
else
    shadow_fading_used = false;
end

if shadow_fading_used && LTE_config.debug_level>=1
    fprintf('Calculating average sector capacity (macroscopic and shadow fading)\n');
else
    fprintf('Calculating average sector capacity (macroscopic fading)\n');
end

%% Preallocate for the SINR matrices
num_eNodeBs = length(eNodeBs);
% For hexagonal grid:
num_macro_eNodeBs = sum(strcmp({eNodeBs.site_type},'macro'));
num_macro_sectors = num_macro_eNodeBs * length(LTE_config.sector_azimuths);


% Look for the eNodeB closer to the center (0,0) and store the eNodeBs' pixel position (for plotting purposes)
eNodeB_pixel_pos = zeros(num_eNodeBs,2);
RX_powers_W = zeros(size(networkPathlossMap.pathloss));
for b_ = 1:num_eNodeBs
    eNodeB_pixel_pos(b_,:) = LTE_common_pos_to_pixel(eNodeBs(b_).pos,networkPathlossMap.coordinate_origin,networkPathlossMap.data_res);
    
    % Matrix containing the received power
    if shadow_fading_used
        shadow_fading_new_map_size   = [size(networkPathlossMap.pathloss,1) size(networkPathlossMap.pathloss,2)];
        if isa(networkShadowFadingMap,'channel_gain_wrappers.shadowFadingDummyMap')
            shadow_fading_current_eNodeB = 10.^(imresize(0,shadow_fading_new_map_size)/10);
        else
            if networkShadowFadingMap.oneMapPerSite
                shadow_fading_current_eNodeB = 10.^(imresize(networkShadowFadingMap.pathloss(:,:,b_),shadow_fading_new_map_size)/10);
            end
        end
    else
        shadow_fading_current_eNodeB = 1;
    end
    for s_ = 1:length(eNodeBs(b_).sectors)
        cellIdx = eNodeBs(b_).sectors(s_).eNodeB_id;
        if shadow_fading_used&&~isa(networkShadowFadingMap,'channel_gain_wrappers.shadowFadingDummyMap')&&~networkShadowFadingMap.oneMapPerSite
            shadow_fading_current_eNodeB = 10.^(imresize(networkShadowFadingMap.pathloss(:,:,cellIdx),shadow_fading_new_map_size)/10);
        end
        s_idx = eNodeBs(b_).sectors(s_).eNodeB_id;
        RX_powers_W(:,:,s_idx) = eNodeBs(b_).sectors(s_).max_power./10.^(networkPathlossMap.pathloss(:,:,s_idx)/10) ./ shadow_fading_current_eNodeB;
    end
end

%% Calculate SINR map for all sectors
SINR_linear_all = zeros(size(RX_powers_W));
SNR_linear_all  = zeros(size(RX_powers_W));
thermal_noise_W = 10^(LTE_config.UE.thermal_noise_density/10) / 1000 * LTE_config.bandwidth * 10^(LTE_config.UE.receiver_noise_figure/10);

tot_eNodeBs = size(RX_powers_W,3);
for s_=1:tot_eNodeBs
    SNR_linear_all(:,:,s_)  = RX_powers_W(:,:,s_) ./ (thermal_noise_W);
    SINR_linear_all(:,:,s_) = RX_powers_W(:,:,s_) ./ (sum(RX_powers_W,3) + thermal_noise_W - RX_powers_W(:,:,s_));
end
SINR_dB_all = 10*log10(SINR_linear_all);

% Calculate the matrix needed to show the SINR difference map
[SINR_dB_all_sorted, SINR_IX] = sort(SINR_dB_all,3);

max_SINR_dB_all    = SINR_dB_all_sorted(:,:,end);
if tot_eNodeBs>1
    diff_SINR_dB_all   = SINR_dB_all_sorted(:,:,end)-SINR_dB_all_sorted(:,:,end-1);
else
    diff_SINR_dB_all = nan(size(SINR_dB_all_sorted(:,:,1)));
end

maxSINR_assignment = SINR_IX(:,:,end);
sector_assignment  = maxSINR_assignment;


%% Calculate sector sizes
cell_sizes = zeros(1,length(eNodeBs_sectors));
cell_centers_pixel = zeros(length(eNodeBs_sectors),2);
   for s_idx = 1:length(eNodeBs_sectors)
       cell_sizes(s_idx) = sum(sector_assignment(:)==s_idx);
       [row,col] = find(sector_assignment==s_idx);
       cell_centers_pixel(s_idx,:) = [mean(col) mean(row)];
   end
 
cell_centers = LTE_common_pixel_to_pos(cell_centers_pixel,networkPathlossMap.coordinate_origin,networkPathlossMap.data_res);

SINR_dB_ROI  = max_SINR_dB_all(:);
SINR_lin_ROI = 10.^(SINR_dB_ROI/10);



%% Calculate average capacity (finally!!!)
bandwidth       = LTE_config.N_RB*LTE_config.RB_bandwidth;
CP_length_s     = LTE_config.CP_length_samples/LTE_config.fs;
symbol_length_s = LTE_config.TTI_length/(LTE_config.N_sym*2);
CP_ratio        = 1-(CP_length_s/symbol_length_s);

nTXantennas        = 1;
subcarriers_per_RB = 12;
switch nTXantennas
    case 1
        nRef_sym = 4;
    case 2
        nRef_sym = 8;
    case 4
        nRef_sym = 12;
end
subframe_size_Sym = LTE_config.N_sym*subcarriers_per_RB*2*LTE_config.N_RB;       % 2 for 2 slots (2x0.5 ms)

RefSym_ratio  = 1-(nRef_sym / (LTE_config.N_sym*subcarriers_per_RB*nTXantennas)); % Ratio of reference_symbols/total_subframe_symbols
SyncSym_ratio = 1-(72 / (subframe_size_Sym*5)); % 72 symbols used for sync every 5 subframes

% Integrate over all of the ROI (sum). Apply correction factors for used bandwidth, Cyclic Prefix and reference/sync symbols.
sector_capacity_vec     = bandwidth*CP_ratio*RefSym_ratio*SyncSym_ratio*log2(1+SINR_lin_ROI);

sector_avg_capacity_mbps = mean(sector_capacity_vec) / 1e6;
sector_min_capacity_mbps = min(sector_capacity_vec) / 1e6;
sector_max_capacity_mbps = max(sector_capacity_vec) / 1e6;

sector_capacity.avg_mbps = sector_avg_capacity_mbps;
sector_capacity.min_mbps = sector_min_capacity_mbps;
sector_capacity.max_mbps = sector_max_capacity_mbps;

end
