function [ UEs, extra_info] = LTE_init_generate_users_and_add_schedulers(LTE_config,sites,eNodeBs,networkPathlossMap,CQI_mapper,BLER_curves,networkClock)

%function [ UEs, extra_info ] = LTE_init_generate_users_and_add_schedulers(LTE_config,sites,eNodeBs,networkPathlossMap,CQI_mapper,BLER_curves,networkClock)
%% Data needed also for plotting: generate every time
fprintf('Generating UEs: %s\n',LTE_config.UE_distribution);
switch LTE_config.UE_distribution   
    case 'LTE_V UEs'
        UE_spatial_distribution = spatial_distributions.LTE_V_UEsSpatialDistribution(networkPathlossMap, LTE_config.vehicle_UEs,...
                                                                                            LTE_config.inter_eNodeB_distance,...
                                                                                            LTE_config.LTE_V_config.inter_eNodeB_distance,...                                                                         
                                                                                            LTE_config.lane_width);
    otherwise
        error('UE distribution %s not supported.',LTE_config.UE_distribution);
end

UE_positions = UE_spatial_distribution.generate_positions;

      
%% Creating UE position, depending on the configuration
% Create UEs according to the previously generated positions
UEs = network_elements.UE;
    for u_ = 1:size(UE_positions,1)
        % General UE settings that can be saved and re-used
        UEs(u_)     = network_elements.UE;
        UEs(u_).id  = u_;
        UEs(u_).pos = [UE_positions(u_,1) UE_positions(u_,2)];
        
        % Add UE move model. UEs are move with fixed speed and angle
        UEs(u_).walking_model = walking_models.straightWalkingModel(LTE_config.UE_speed*LTE_config.TTI_length, LTE_config.UE_direction);
        % Add UEs moving direction
        UEs(u_).walking_model.direction= UE_positions(u_,3);
        UEs(u_).trace_UE = false;                             
    end
  
              
       %% add eNodeBs Functions to Vehicles
        if LTE_config.add_eNodeB_to_Vehicle
            [UEs, VeNodeBs, VeNodeBs_sectors, networkPathlossMap]= network_generation.add_Vehicle_eNodeBs(LTE_config,UEs,sites,eNodeBs,networkPathlossMap);
        end
        
%% Steps needed for FFR
if LTE_config.FFR_active
    switch LTE_config.shadow_fading_type
        case 'none'
            % OK
        otherwise
            error('Right now, due to how the R3 frequency assignment is realized, just simulations WITHOUT shadow fading are supported for FFR simulations');
    end
    
    site_pos                        = reshape([sites.pos],2,[])';
    ROI_center                      = [mean(networkPathlossMap.roi_x) mean(networkPathlossMap.roi_y)];
    distance_to_center              = sqrt(sum([site_pos(:,1)-ROI_center(1) site_pos(:,2)-ROI_center(2)].^2,2));
    [null_var, min_distance_eNodeB] = min(distance_to_center);
    
    % Frequency assignment
    LTE_config.scheduler_params.frequency_assignment = utils.ffrUtils.assign_frequencies_to_hex_grid(sites(min_distance_eNodeB).sectors(1).eNodeB_id,eNodeBs,networkPathlossMap.sector_centers);
    
    if ~LTE_config.FFR_override
        tx_mode_string_long = utils.miscUtils.tx_mode_to_string_long(LTE_config.tx_mode,LTE_config.nTX,LTE_config.nRX);
        [optimum_BFR, optimum_FR_SINR_switching_dB, ~] = utils.ffrUtils.load_FFR_optimum_BFRs(tx_mode_string_long);
    else
        optimum_BFR                  = LTE_config.FFR_params.beta_FR;
        optimum_FR_SINR_switching_dB = LTE_config.FFR_params.SINR_threshold_value;
    end
    
    %% Calculate frequency assignment
    FFR_UE_mapping = utils.ffrUtils.assign_FFR_band_to_UEs(UEs,optimum_FR_SINR_switching_dB,networkPathlossMap);    
    LTE_config.scheduler_params.FFR_UE_mapping = FFR_UE_mapping;
    LTE_config.scheduler_params.beta_FR        = optimum_BFR;
end

%% Initialize schedulers (eNB)
LTE_init_add_schedulers(LTE_config,sites,UEs,CQI_mapper,BLER_curves);

%% Initialize schedulers (UE)

%% Other UE initialization, including adding a downlink and uplink channel object to each user
% The downlink will contain pathloss maps, so depending on the user's position, it will 'see' a certain pathloss.
% Add also the penetration loss and noise figure.
% The uplink is simply a delay between the UE and the eNodeB.

for u_=1:length(UEs)
    
    % Add receiver antenna gain
    UEs(u_).antenna_gain = LTE_config.UE.antenna_gain;
    
    % Add noise figure
    UEs(u_).receiver_noise_figure = LTE_config.UE.receiver_noise_figure;
    
    % Thermal noise (receiver) for the link quality model (in linear: watts)
    UEs(u_).thermal_noise_W_RB = 10^(0.1*LTE_config.UE.thermal_noise_density)/1000 * LTE_config.RB_bandwidth * 10^(UEs(u_).receiver_noise_figure/10);
    
    % Default tx mode for feedback (for the old trace format -v1- this
    % sets the only tx mode that can be used)
    UEs(u_).default_tx_mode = LTE_config.tx_mode;
       
    % Set signaling channel (eNodeB to UE)
    UEs(u_).eNodeB_signaling = network_elements.eNodebSignaling;
    
    % Set signaling channel (eNodeB to UE)
    UEs(u_).VeNodeB_signaling = network_elements.eNodebSignaling;
    
    % Number of RX antennas
    UEs(u_).nRX = LTE_config.nRX;
    
    % Set BLER curves for ACK/NACK calculation
    UEs(u_).BLER_curves = BLER_curves;
    
    % Clock
    UEs(u_).clock = networkClock;
    
    % CQI mapper
    UEs(u_).CQI_mapper = CQI_mapper;
        
    % Configure extra tracing
    UEs(u_).trace_SINR = LTE_config.trace_SINR;
    
    % Adaptive RI
    UEs(u_).adaptive_RI = LTE_config.adaptive_RI;
        
end
%% Safeguard against having no UEs
if length(UEs)==1 && isempty(UEs(1).id)
    no_UEs = true;
else
    no_UEs = false;
end

if ~no_UEs
    % Assign the UEs to their nearest (in terms of SINR) eNodeB and assign some extra parameters
    sector_UE = false(1,3);
    for u_ = 1:length(UEs)        
        if ~UEs(u_).trace_UE
            % Attach UE to eNodeB
            [ site_id, sector_num, eNodeB_id] = networkPathlossMap.cell_assignment(UEs(u_).pos); %#ok<ASGLU>
            eNodeBs(eNodeB_id).attachUser(UEs(u_));
        else
            % Do not attach the UEs yet (do that when we move and/or activate them)
        end
              
        % Check whether this UE should be deactivated to speed-up simulation
        if ~isempty(LTE_config.compute_only_UEs_from_this_eNodeBs)
            if isempty(find(UEs(u_).attached_eNodeB.eNodeB_id==LTE_config.compute_only_UEs_from_this_eNodeBs,1))
                % Deactivate UE
                UEs(u_).deactivate_UE = true;
            else
                % Activate UE (already activated by default, but just in case)
                UEs(u_).deactivate_UE = false;
            end
        end 
                
        % Append traffic model to users
            UEs(u_).traffic_model = LTE_trafficmodel(LTE_config.traffic_models,UEs(u_),max(LTE_config.feedback_channel_delay,0));            
    end 
      
            % Choose as many points per cell as users
        if ~exist('extra_info','var')
            extra_info = [];
        end
       
end

function LTE_init_add_schedulers(LTE_config,sites,UEs,CQI_mapper,BLER_curves)
% Adds the needed scheduler type and resource block grid to each eNodeb's sector
% input:   eNodeBs  ... array of eNodeBs
%          UEs      ... array of UEs

% Check whether this scheduler exists
schedulers.schedulerFactory.chek_whether_scheduler_is_defined(LTE_config.scheduler);

if LTE_config.debug_level>=1
    switch LTE_config.scheduler
        case 'FFR'
            fprintf('Creating %s scheduler (FR: %s, PR: %s) and resource block grids\n',LTE_config.scheduler,LTE_config.scheduler_params.FR_scheduler.scheduler,LTE_config.scheduler_params.PR_scheduler.scheduler);
        otherwise
            fprintf('Creating %s schedulers and resource block grids\n',LTE_config.scheduler);
    end
end

% No reason to use a different SINR averager instance for each scheduler, we can reuse the same one

switch LTE_config.SINR_averaging.algorithm
    case 'MIESM'
        the_SINR_averager = utils.miesmAveragerFast(LTE_config,LTE_config.SINR_averaging.BICM_capacity_tables,LTE_config.SINR_averaging.betas);
    case 'EESM'
        error('EESM SINR averaging is no longer supported supported');
    otherwise
        error('SINR averaging algorithm not supported');
end

% Add RB grid representation and scheduler to each sector.
% Set also homogeneous power load
for b_ = 1:length(sites)
    for s_=1:length(sites(b_).sectors)
        
        % Set whether the eNodeBs will always transmit, even if no UEs are attached.
        sites(b_).sectors(s_).always_on = LTE_config.always_on;
        
        max_data_power  = sites(b_).sectors(s_).max_power;
        signaling_power = sites(b_).sectors(s_).signaling_power;
        
        LTE_config.scheduler_params.max_power       = max_data_power; % For backwards compatibility
        LTE_config.scheduler_params.CQI_params      = LTE_config.CQI_params;
        LTE_config.scheduler_params.default_tx_mode = LTE_config.tx_mode;
        
        % RB grid creation and initialization
        sites(b_).sectors(s_).RB_grid = network_elements.resourceBlockGrid(LTE_config.N_RB,LTE_config.sym_per_RB_nosync,LTE_config.sym_per_RB_sync);
        sites(b_).sectors(s_).RB_grid.set_homogeneous_power_allocation(sites(b_).sectors(s_).max_power,sites(b_).sectors(s_).signaling_power);
        
        % Continue with Scheduler initialization
        sites(b_).sectors(s_).scheduler = schedulers.schedulerFactory.create_scheduler(LTE_config.scheduler,LTE_config.scheduler_params,sites(b_).sectors(s_));
        
        % Set scheduler SINR averaging algorithm
        sites(b_).sectors(s_).scheduler.set_SINR_averager(the_SINR_averager);

        % Other data required to perform SINR averaging at the transmitter side
        sites(b_).sectors(s_).scheduler.set_CQI_mapper(CQI_mapper);
        sites(b_).sectors(s_).scheduler.set_BLER_curves(BLER_curves);
        
        % Add genie information
        sites(b_).sectors(s_).scheduler.set_genie_UEs(UEs);
        sites(b_).sectors(s_).scheduler.set_genie_eNodeBs(sites);
        
        % Add TTI delay information
        sites(b_).sectors(s_).scheduler.set_feedback_delay_TTIs(LTE_config.feedback_channel_delay);
    end
end
