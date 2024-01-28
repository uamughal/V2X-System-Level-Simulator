classdef LTE_V_ < utils.simulatorConfig
    methods (Static)
        function LTE_config = apply_parameters(LTE_config)
            % Load the LTE System Level Simulator config parameters
            % (c) Josep Colom Ikuno, INTHFT, 2008
            % www.nt.tuwien.ac.at
            
            %% Debug options
            LTE_config.debug_level = 1;  % 0=no output
            % 1=basic output
            % 2=extended output
            
            %% Plotting options
            LTE_config.show_network = 1; % 0= show no plots
            % 1= show some plots
            % 2= show ALL plots (moving UEs)
            % 3=plot even all of the pregenerated fast fading
            
            %% General options
            LTE_config.frequency       = 2e9;         % Frequency in Hz
            LTE_config.frequency1       = 6e9;         % Frequency in Hz         
            LTE_config.bandwidth       = 10e6;           % Frequency in Hz
            
            LTE_config.nTX           = 1;
            LTE_config.nRX           = 1;
            LTE_config.tx_mode       = 1; 
            
            LTE_config.always_on     = true; % Controls whether the eNodeB is always radiating power (default and worse-case scenario) or no power is used when no UEs are attached
            
            %% Random number generation options
            LTE_config.seedRandStream = false;
            LTE_config.RandStreamSeed = 0;      % Only used if the latter is set to 'true'
            
            %% Simulation time
            LTE_config.simulation_time_tti =1; % Simulation time in TTIs
            
            %% Cache options. Saves the generated eNodeBs, Pathloss map and Shadow fading map to a .mat file
            LTE_config.cache_network = true;
            LTE_config.network_cache = 'auto';
            LTE_config.delete_ff_trace_at_end = true; % Reduces the amount needed to store the traces by deleting the fading parameters trace from the results file
            LTE_config.delete_pathloss_at_end = true; % Further reduces the amount of space needed to store the traces by deleting the pathloss maps from the results file
            LTE_config.UE_cache      = false;  % Option to save the user position to a file. This works in the following way:
            %   - cache=true and file exists: read position from file
            %   - cache=true and file does not exist: create UEs and save to cache
            %   - cache=false: do not use cache at all
            LTE_config.UE_cache_file = 'auto';
            % LTE_config.UE_cache_file = fullfile('./data_files/UE_cache_2rings_2UEs_sector_20100824_174824.mat');
            
            %% How to generate the network. If the map is loaded, this parameters will be overridden by the loaded map
            % Use network planning tool Capesso by placing 'capesso'
            LTE_config.network_source = 'generated';
             LTE_config.network_geometry = 'regular_hexagonal_grid'; % Geometry refers to spatial distribution of macro sites
            
             % Configure the network source. Overridden if a pregenerated network pathloss map is used.
            switch LTE_config.network_source
                case 'generated'
                    % Network size
                    LTE_config.inter_eNodeB_distance = 500; % In meters. When the network is generated, this determines the
                    % distance between the eNodeBs.
                    
                    LTE_config.average_eNodeB_distance = LTE_config.inter_eNodeB_distance; % what??????
                    LTE_config.map_resolution = 5; % In meters/pixel. Also the resolution used for initial user creation
                    LTE_config.nr_eNodeB_rings = 1; % Number of eNodeB rings
                    LTE_config.minimum_coupling_loss = 70; % Minimum Coupling Loss: the parameter describing the minimum
                    % loss in signal [dB] between BS and UE or UE and UE in the worst case and is defined as the minimum distance loss including
                    % antenna gains measured between antenna connectors.
                    % Recommended in TS 36.942 are 70 dB for urban areas, 80 dB for rural.
                    %% 
                    LTE_config.sector_azimuths         = 60:360/3:359;  % for PS_LTE need modification % what is meaning?
                    LTE_config.macroscopic_pathloss_model = 'Urban macro-cell'; % I need to change to winner II pathloss model
                    LTE_config.macroscopic_pathloss_model_settings.environment = 'LOS';
                    
                    % eNodeB TX power settings
                    LTE_config.eNodeB_tx_power = 10^(46/10)*1/1000; % macro eNodeB's transmit power, in Watts.
        
                     %% Options for adding eNodeB Proporties to each Vehicle UEs
                   LTE_config.add_eNodeB_to_Vehicle                            = true;
                   
                   LTE_config.LTE_V_config.eNB_distribution                     = 'four_side_of_lane';   %'four_side_of_road_line'
                   LTE_config.LTE_V_config.inter_eNodeB_distance                    = 100;  % Femtocell density in macro cell sector
                   LTE_config.LTE_V_config.tx_power_W                               = 10^(23/10)*1/1000;  % Power in Watts 23dBm
                   LTE_config.LTE_V_config.sector_azimuths                          = [0 120 240];
                   LTE_config.LTE_V_config.tx_height                                = 1.5; % 1.5m
                   LTE_config.LTE_V_config.nTX                                      = 1; % default
                   
                   LTE_config.LTE_V_config.macroscopic_pathloss_model               = 'Urban macro-cell';
                   LTE_config.LTE_V_config.macroscopic_pathloss_model_settings.environment = 'LOS';
                   LTE_config.LTE_V_config.minimum_coupling_loss     =70;
                                                                                                                                                        
                %% Generation of the shadow fading
                 LTE_config.shadow_fading_type = 'claussen'; % Right now only 2D space-correlated shadow fading maps implemented
            
                 % Configure the network source
                 switch LTE_config.shadow_fading_type
                   case 'claussen'
                    LTE_config.shadow_fading_map_resolution = 5;
                    LTE_config.shadow_fading_n_neighbors    = 12;
                    LTE_config.shadow_fading_mean           = 0;
                    LTE_config.shadow_fading_sd             = 4;
                    LTE_config.r_eNodeBs                    = 0.5; % inter-site shadow fading correlation
                end
            
            %% Microscale Fading Generation config
            % Microscale fading trace to be used between the eNodeB and its attached UEs.
            LTE_config.channel_model.type = 'winner+'; % 'winner+' 'PedB' 'extPedB' 'TU' --> the PDP to use
            LTE_config.channel_model.type = 'winner+'; % 'winner+' 'PedB' 'extPedB' 'TU' --> the PDP to use
            LTE_config.channel_model.trace_length = 5; % Length of the trace in seconds. Be wary of the size you choose, as it will be loaded in memory.
            LTE_config.channel_model.correlated_fading = true;
            LTE_config.pregenerated_ff_file            = 'auto';
            % With this option set to 'true', even if cache is present, the channel trace will be recalculated
            LTE_config.recalculate_fast_fading = false;
            
            %% UE (users) settings
            % note that for reducing trace sizes, the UE_id is stored as a uint16, so
            % up to 65535 users in total are supported. To change that, modify the scheduler class.
            
            LTE_config.UE.receiver_noise_figure = 9;    % Receiver noise figure in dB
            LTE_config.UE.thermal_noise_density = -174;  % Thermal noise density in dBm/Hz (-174 dBm/Hz is the typical value)
            LTE_config.UE_distribution = 'LTE_V UEs';
            %LTE_config.UE_per_eNodeB   = 2;     % number of users per eNodeB sector (calculates it for the center sector and applies this user density to the other sectors)
            LTE_config.UE_speed        = 60/3.6; % Speed at which the UEs move. In meters/second: 5 Km/h = 1.38 m/s
            LTE_config.keep_UEs_still  = false;
            %% eNodeB options
            LTE_config.antenna.antenna_gain_pattern = 'TS 36.942 3D'; % As defined in TS 36.942. Identical to Berger, but with a 65¡ã 3dB lobe       
            LTE_config.antenna.max_antenna_gain = 15; % LTE antenna, urban area (2000 MHz)
              
            %% Additional parameters needed when using 'kathreinTSAntenna' and partially by the 'TS 36.942 3D' anntennas
            switch LTE_config.antenna.antenna_gain_pattern
                case {'TS 36.942 3D'}
                    % Default values used for each site
                    LTE_config.site_altitude               = 0;    % Altiude of site [m] (doesn't make sense if you use an elevation map, thought)
                    LTE_config.tx_height                   = 18;   % Height of transmitter  [m]
                    LTE_config.rx_height                   = 1.5;  % Receiver height [m]
                    LTE_config.antenna.mechanical_downtilt = 0; % [°]
                    LTE_config.antenna.electrical_downtilt = 8; % [°]   
                    
                    %%Vehicle UEs
                    % Default values used for each site
                    LTE_config.LTE_V_config.antenna.antenna_gain_pattern='TS 36.942 3D';
                    LTE_config.LTE_V_config.antenna.max_antenna_gain = 3; %3
                    LTE_config.LTE_V_config.site_altitude               = 0;    % Altiude of site [m] (doesn't make sense if you use an elevation map, thought)
                    LTE_config.LTE_V_config.tx_height                   = 1.5;   % Height of transmitter 1.5 [m]
                    LTE_config.LTE_V_config.rx_height                   = 1.5;  % Receiver height1.5 [m]
                    LTE_config.LTE_V_config.antenna.mechanical_downtilt = 0; % [°]
                    LTE_config.LTE_V_config.antenna.electrical_downtilt = 1; % [°]
                                                              
            end 
            %% Scheduler options
            LTE_config.scheduler        = 'prop fair traffic';  %'broadcast' 'round robin', 'best cqi', 'max min', 'max TP', 'resource fair', 'proportional fair' or 'prop fair Sun'
           % LTE_config.scheduler        = 'prop fair traffic';  
            LTE_config.scheduler_params.fairness = 0.5; % fairness for the variable fairness scheduler
            LTE_config.power_allocation = 'homogeneous;'; % 'right now no power loading is implemented, so just leave it as 'homogeneous'
            
            %% CQI mapper options
            LTE_config.CQI_mapper.CQI2SNR_method = 1; % 1 to use less conservative method to map from CQI back to SNR
            % uses the value in the middle of the SNR interval corresponding to a CQI instead of the lower boarder value
            % this value will just be used in connection with quantized CQI feedback  
            %% Uplink channel options
            LTE_config.feedback_channel_delay   = 1; % In TTIs
            LTE_config.unquantized_CQI_feedback = false;
            
            %% Where to save the results
            LTE_config.results_folder         = './results';
            LTE_config.results_file           = 'auto'; % NOTE: 'auto' assigns a filename automatically
            
            %% Values that should not be changed
            LTE_config.antenna_azimuth_offsett = 30;    % This controls the antenna layout that will be generated. 0 degrees generates hexagonal cells,
            % while 30 degrees hexagonal sectors.
            end
        end
    end
end

