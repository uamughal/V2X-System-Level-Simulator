function [eNodeBs, eNodeBs_sectors, networkPathlossMap, networkShadowFadingMap] = LTE_init_network_generation(LTE_config,varargin)
% Network generation. Either from a file (cache) or calling the necessary function.
% (c) Josep Colom Ikuno, INTHFT, 2009
% www.nt.tuwien.ac.at

if ~LTE_config.trace_simulation_mode
    vararginIsEmpy = isempty(varargin);
    if ~vararginIsEmpy
        pathlossDataIsEmpty = isempty(varargin{1}{2});
    else
        pathlossDataIsEmpty = true;
    end
    
    if LTE_config.usePathlossMapFromMemory&&pathlossDataIsEmpty
        error('When setting LTE_config.usePathlossMapFromMemory to true, pathloss data must be input.\n');
    end
%% Generate network (eNodeBs and macroscopic pathloss)
        % generate networkPathlossMap
        if LTE_config.debug_level>=1
            sprintf('Generating network\n');
        end
        switch LTE_config.network_source
            case 'generated'
                [eNodeBs, eNodeBs_sectors, networkPathlossMap] = network_generation.generated_network(LTE_config);
        end
                 
        % Add the power separation. X% to signaling/pilots (always on) and the rest for data
        set_signaling_power(LTE_config,eNodeBs);
        
        % Store the other eNodeBs as (potential) interferers
        for s_ = 1:length(eNodeBs_sectors)
            eNodeBs_sectors(s_).neighbors_eNodeB = eNodeBs_sectors([1:(s_-1) (s_+1):length(eNodeBs_sectors)]);
        end
        
        % Generate shadow fading
        if LTE_config.macroscopic_pathloss_is_model
            if LTE_config.debug_level>=1
                fprintf('Generating shadow fading\n');
            end
            switch LTE_config.shadow_fading_type
                case 'claussen'
                    [LTE_config.roi_x,LTE_config.roi_y] = networkPathlossMap.valid_range;
                    if LTE_config.debug_level>=1
                        fprintf('Generating Claussen space-correlated shadow fading map, ');
                    end
                    if ~LTE_config.decouple_site_shadow_fading_maps
                        fprintf('one map per site\n');
                        networkShadowFadingMap = channel_gain_wrappers.shadowFadingMapClaussen(...
                            LTE_config.shadow_fading_map_resolution,...
                            LTE_config.roi_x,...
                            LTE_config.roi_y,...
                            LTE_config.shadow_fading_n_neighbors,...
                            length(eNodeBs),...
                            LTE_config.shadow_fading_mean,...
                            LTE_config.shadow_fading_sd,...
                            LTE_config.r_eNodeBs,...
                            LTE_config.deactivate_claussen_spatial_correlation);
                        networkShadowFadingMap.oneMapPerSite = true;
                    else
                        fprintf('one map per cell\n');
                        networkShadowFadingMap = channel_gain_wrappers.shadowFadingMapClaussen(...
                            LTE_config.shadow_fading_map_resolution,...
                            LTE_config.roi_x,...
                            LTE_config.roi_y,...
                            LTE_config.shadow_fading_n_neighbors,...
                            length(eNodeBs_sectors),...
                            LTE_config.shadow_fading_mean,...
                            LTE_config.shadow_fading_sd,...
                            LTE_config.r_eNodeBs,...
                            LTE_config.deactivate_claussen_spatial_correlation);
                        networkShadowFadingMap.oneMapPerSite = false;
                    end
            end
        end
        
        %% Calculate the SINR for each sector based on the pathloss and maximum TX power alone
        if LTE_config.macroscopic_pathloss_is_model
            [networkPathlossMap.capacity,...
                networkPathlossMap.SINR,...
                networkPathlossMap.sector_assignment,...
                networkPathlossMap.maxSINR_assignment,...
                networkPathlossMap.diff_SINR_dB,...
                networkPathlossMap.sector_sizes,...
                networkPathlossMap.sector_centers] = LTE_common_calculate_cell_capacity(LTE_config,networkPathlossMap,eNodeBs,eNodeBs_sectors,networkShadowFadingMap);
            [networkPathlossMap.capacity2,...
                networkPathlossMap.SINR2,...
                networkPathlossMap.sector_assignment2,...
                networkPathlossMap.maxSINR_assignment,...
                networkPathlossMap.diff_SINR_dB2,...
                networkPathlossMap.sector_sizes2,...
                networkPathlossMap.sector_centers2] = LTE_common_calculate_cell_capacity(LTE_config,networkPathlossMap,eNodeBs,eNodeBs_sectors);
        else
            [networkPathlossMap.capacity,...
                networkPathlossMap.SINR,...
                networkPathlossMap.sector_assignment,...
                networkPathlossMap.maxSINR_assignment,...
                networkPathlossMap.diff_SINR_dB,...
                networkPathlossMap.sector_sizes,...
                networkPathlossMap.sector_centers] = LTE_common_calculate_cell_capacity(LTE_config,networkPathlossMap,eNodeBs,eNodeBs_sectors);
        end
        
        % Save network
        if LTE_config.cache_network
            try
                [pathstr, name, ext]     = fileparts(LTE_config.network_cache);
                hashtag_network_cache                 = '';       
                LTE_config.network_cache = fullfile(pathstr,[name hashtag_network_cache ext]);
                
                if LTE_config.debug_level>=1
                    fprintf('Saving network to file: %s\n',LTE_config.network_cache);
                end
                
                if exist(LTE_config.network_cache,'file')
                    throw(MException('LTEsim:cacheExists', 'The cache file was concurrently generated during another simulation run'));
                end
                
                if exist('networkShadowFadingMap','var')
                    save(LTE_config.network_cache,'eNodeBs','networkPathlossMap','networkShadowFadingMap','eNodeBs_sectors');
                else
                    save(LTE_config.network_cache,'eNodeBs','networkPathlossMap','eNodeBs_sectors');
                end
            catch err
                fprintf('Network cache could not be saved. If needed, it will be generated again in the next run (%s).\n',err.message);
            end
        end
    end


% To avoid error of missing return argument
if ~exist('networkShadowFadingMap','var')
    networkShadowFadingMap = [];
end

function set_signaling_power(LTE_config,eNodeBs)
% Set signaling power
for b_=1:length(eNodeBs)
    for s_=1:length(eNodeBs(b_).sectors)
        data_power      = eNodeBs(b_).sectors(s_).max_power * (1-LTE_config.signaling_ratio);
        signaling_power = eNodeBs(b_).sectors(s_).max_power * LTE_config.signaling_ratio;
        eNodeBs(b_).sectors(s_).max_power       = data_power;
        eNodeBs(b_).sectors(s_).signaling_power = signaling_power;
        LTE_config.scheduler_params.max_power   = data_power; % max data transmit power in Watts
    end
end

