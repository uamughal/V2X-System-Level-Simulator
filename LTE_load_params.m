function LTE_config = LTE_load_params(varargin)

%% Default simulation option
default_simulation = 'LTE_V';

%% If no simulation option is defined, use the default one
if ~isempty(varargin)
    simulation_type = varargin{1};
else
    simulation_type = default_simulation; %  Default value
end

fprintf('Using "%s" simulation configuration.\n',simulation_type)

%% Load the corresponding simulation parameters
switch simulation_type
    case 'LTE_V'
        LTE_config = simulation_config.LTE_V_.apply_parameters;
    otherwise
        warning('Simulation type not defined: using default one instead');
        simulation_type = default_simulation;
        LTE_config = simulation_config.hex_grid_tilted.apply_parameters;
end
