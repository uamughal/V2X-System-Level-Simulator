function pregenerated_fast_fading = LTE_init_get_microscale_fading_SL_trace(LTE_config)
% Generate the fading parameters that model the fast (microscale) fading at system level.

%% Config

% Possible Tx modes are:
%   1: Single Antenna
%   2: Transmit Diversity
%   3: Open Loop Spatial Multiplexing
%   4: Closed Loop SM
% Number of antenna ports can either be 2 or 4
% nLayers specifies how many layers (symbols) are going to be transmitted.
% Either 1, 2, 3 or 4

standalone = false;

if standalone
    LTE_config.debug_level            = 1;
    % Initial config (Winner II+)
    config.system_bandwidth           = 1.4e6;
    config.channel_type               = 'winner+';
    config.nTX                        = 2;
    config.nRX                        = 2;
    config.trace_length_s             = 1;
    config.UE_speed                   = 400/3.6; % converted to m/s. High speed ~ uncorrelated
    config.parallel_toolbox_installed = true; % change it to false for testing purposes, as if not you will not be able to debug properly
    config.feedback_channel_delay     = 0;
    config.correlated_fading          = true;
    config.f                          = 2.1400e+009;
    config.TTI_length                 = 1.0000e-003;
    config.tx_mode                    = 4; % OLSM
    
    winner_config                     = load('winner_trace_config');
    config.trace_params               = winner_config.config_trace_params;
    config.trace_params.speed         = config.UE_speed;
    clear winner_config
else
    % Initial config (simulator-linked)
    config.system_bandwidth           = LTE_config.bandwidth;
    config.channel_type               = LTE_config.channel_model.type;
    config.nTX                        = LTE_config.nTX;
    config.nRX                        = LTE_config.nRX;
    config.trace_length_s             = LTE_config.channel_model.trace_length;
    config.UE_speed                   = LTE_config.UE_speed; % converted to m/s
    config.parallel_toolbox_installed = LTE_config.parallel_toolbox_installed; % change it to false for testing purposes, as if not you will not be able to debug properly
    config.feedback_channel_delay     = LTE_config.feedback_channel_delay;
    config.correlated_fading          = LTE_config.channel_model.correlated_fading;
    config.f                          = LTE_config.frequency;
    config.trace_params               = LTE_config.trace_params;
    config.TTI_length                 = LTE_config.TTI_length;
    config.tx_mode                    = LTE_config.tx_mode;
    config.non_parallel_channel_trace = LTE_config.non_parallel_channel_trace;
    
    % Option for wideband precoding
    if LTE_config.tx_mode==4
        config.wideband_precoding = LTE_config.wideband_precoding;
    end
end

% sigma_n2 = 10^((LTE_config.UE.receiver_noise_figure + LTE_config.UE.thermal_noise_density)/10)/1000;    % Receiver noise variance in Watt

% We now have all of the possible precoding combinations stored
precoding_configs = phy_modeling.miscUtils.get_all_precoding_combinations;

% Channel trace for the target and interfering channels
switch config.channel_type
    case 'winner+'
        channel_factory_H0 = channel_gain_wrappers.winnerChannelFactory(config.system_bandwidth,config.trace_params,LTE_config.winner_antenna_params);
        channel_factory_H1 = channel_gain_wrappers.winnerChannelFactory(config.system_bandwidth,config.trace_params,LTE_config.winner_antenna_params);
end
if LTE_config.debug_level>=1
    fprintf('Generating %dx%d channel trace of length %3.2fs\n',config.nTX,config.nRX,ceil(config.trace_length_s));
end
H_trace0 = channel_factory_H0.generate_FF_trace(config.trace_length_s/config.TTI_length);

% Interfering channel trace
if LTE_config.debug_level>=1
    fprintf('Generating %dx%d interfering channel trace of length %3.2fs\n',config.nTX,config.nRX,ceil(config.trace_length_s));
end
H_trace1 = channel_factory_H1.generate_FF_trace(config.trace_length_s/config.TTI_length);


%% Channel normalization

% Note: each MIMO channel is normalized to a mean power of one
H_trace_normalized        = H_trace0.H_RB_samples;
H_trace_interf_normalized = H_trace1.H_RB_samples;

% Free up memory
clear H_trace0;
clear H_trace1;

switch LTE_config.trace_version
    case 'v1'
        % v1 trace
        pregenerated_fast_fading = phy_modeling.channelTraceFactory_v1.generate_channel_trace(config,precoding_configs,H_trace_normalized,H_trace_interf_normalized);
    otherwise
        % v2 trace
        pregenerated_fast_fading = phy_modeling.channelTraceFactory_v2.generate_channel_trace(config,precoding_configs,H_trace_normalized,H_trace_interf_normalized);
end
end

