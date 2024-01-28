function [ sites,eNodeBs ] = LTE_add_schedulers_to_Vehicles(LTE_config,UEs,sites,eNodeBs,CQI_mapper,BLER_curves)
  
%% Initialize schedulers
LTE_init_add_schedulers(LTE_config,sites,UEs,CQI_mapper,BLER_curves);


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
