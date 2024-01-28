% function [UE_throughput] = LTE_show_aggregate_results(simulation_data)
% 
% 
% 
% if ~isempty(cells_to_plot)
%     [UEs_to_use cell_sum_throughput] = utils.resultsFileReader.get_UEs_in_given_cells(cells_to_plot,the_UE_traces);
%     
%     the_UE_traces_to_plot = the_UE_traces(UEs_to_use);
%     wideband_SINRs_all = reshape([the_UE_traces.wideband_SINR],simulation_data.LTE_config.simulation_time_tti,[]);
%     wideband_SINRs_all = wideband_SINRs_all(1,:);
%     
%     % To get the average RB occupancy
%     nRB                     = unique([the_eNodeB_traces(cells_to_plot).RB_grid_size]); % It should be a unique value!
%     RB_occupancy_ratio      = double([the_eNodeB_traces(cells_to_plot).scheduled_RBs])/nRB;
%     mean_RB_occupancy_ratio = mean(RB_occupancy_ratio);
%     
%     % To get the axes limits
%     throughput_Mbps_ecdf_all = utils.miscUtils.ecdf([the_UE_traces.average_throughput_Mbps]);
%     spectral_eff_ecdf_all    = utils.miscUtils.ecdf([the_UE_traces.average_spectral_efficiency_bit_per_cu]);
%     wideband_SINR_ecdf_all   = utils.miscUtils.ecdf(wideband_SINRs_all);
%       
%     % The actual values
%     throughput_Mbps_ecdf = utils.miscUtils.ecdf([the_UE_traces_to_plot.average_throughput_Mbps]);
%     fairness_index       = sum(throughput_Mbps_ecdf.input_data).^2 / sum(throughput_Mbps_ecdf.input_data.^2) / sum(isfinite(throughput_Mbps_ecdf.input_data));
%     spectral_eff_ecdf    = utils.miscUtils.ecdf([the_UE_traces_to_plot.average_spectral_efficiency_bit_per_cu]);
%     wideband_SINR_ecdf   = utils.miscUtils.ecdf(wideband_SINRs_all(UEs_to_use));
%      UE_throughput = throughput_Mbps_ecdf.x;
%     
%     
%     
%     % Wideband SINR-to-throughput plot
%     wideband_SINR_vector                       = wideband_SINR_ecdf.input_data;
%     spectral_eff_vector                        = spectral_eff_ecdf.input_data;
%     throughput_vector                          = throughput_Mbps_ecdf.input_data;
%     [wideband_SINR_binned  throughput_binned   numel_bins] = utils.miscUtils.fit_scatterplot_data(wideband_SINR_vector,throughput_vector,50);
%     [wideband_SINR_binned2 spectral_eff_binned numel_bins] = utils.miscUtils.fit_scatterplot_data(wideband_SINR_vector,spectral_eff_vector,50);
%     
%     % Throughput ECDF
%     hold(UE_throughput_axes,'all');
%     plot(UE_throughput_axes,throughput_Mbps_ecdf.x,throughput_Mbps_ecdf.f,'blue');
%     grid(UE_throughput_axes,'on');
%     title(UE_throughput_axes,'UE average throughput');
%     xlabel(UE_throughput_axes,'average UE throughput [Mb/s]');
%     ylabel(UE_throughput_axes,'F(x)');
%     if constant_axes_limits
%         xlim(UE_throughput_axes,[throughput_Mbps_ecdf_all.min throughput_Mbps_ecdf_all.max]);
%     else
%         xlim(UE_throughput_axes,'auto');
%     end
%     hold(UE_throughput_axes,'off');
%     
%     % Spectral efficiency ECDF
%     hold(UE_spectral_eff_axes,'all');
%     plot(UE_spectral_eff_axes,spectral_eff_ecdf.x,spectral_eff_ecdf.f,'blue');
%     grid(UE_spectral_eff_axes,'on');
%     title(UE_spectral_eff_axes,'UE average spectral efficiency');
%     xlabel(UE_spectral_eff_axes,'average UE spectral efficiency [bit/cu]');
%     ylabel(UE_spectral_eff_axes,'F(x)');
%     if constant_axes_limits
%         xlim(UE_spectral_eff_axes,[spectral_eff_ecdf_all.min spectral_eff_ecdf_all.max]);
%     else
%         xlim(UE_spectral_eff_axes,'auto');
%     end
%     hold(UE_spectral_eff_axes,'off');
%     
%     % Wideband SINR ECDF
%     hold(UE_wideband_SINR_axes,'all');
%     plot(UE_wideband_SINR_axes,wideband_SINR_ecdf.x,wideband_SINR_ecdf.f,'blue');
%     grid(UE_wideband_SINR_axes,'on');
%     title(UE_wideband_SINR_axes,'UE wideband SINR');
%     xlabel(UE_wideband_SINR_axes,'UE wideband SINR [dB]');
%     ylabel(UE_wideband_SINR_axes,'F(x)');
%     if constant_axes_limits
%         xlim(UE_wideband_SINR_axes,[wideband_SINR_ecdf_all.min wideband_SINR_ecdf_all.max]);
%     else
%         xlim(UE_wideband_SINR_axes,'auto');
%     end
%     hold(UE_wideband_SINR_axes,'off');
%     
%     % SINR-to-throughput
%     hold(UE_SINR_to_throughput_axes,'all');
%     if plot_points_scatter
%         if plot_in_GUI
%             scatter(UE_SINR_to_throughput_axes,wideband_SINR_vector,throughput_vector,'.b');
%         else
%             if plot_mean_scatter
%                 scatter(UE_SINR_to_throughput_axes,wideband_SINR_vector,throughput_vector,'.b','SizeData',50);
%             else
%                 scatter(UE_SINR_to_throughput_axes,wideband_SINR_vector,throughput_vector,'.b','SizeData',150);
%             end
%         end
%     end
%     if plot_mean_scatter
%         if plot_in_GUI
%             scatter(UE_SINR_to_throughput_axes,wideband_SINR_binned,throughput_binned,'.r');
%         else
%             scatter(UE_SINR_to_throughput_axes,wideband_SINR_binned,throughput_binned,'.r','SizeData',350);
%         end
%     end
%     grid(UE_SINR_to_throughput_axes,'on');
%     title(UE_SINR_to_throughput_axes,'UE wideband SINR-to-throughput mapping');
%     xlabel(UE_SINR_to_throughput_axes,'UE wideband SINR [dB]');
%     ylabel(UE_SINR_to_throughput_axes,'average UE throughput [Mb/s]');
%     if constant_axes_limits
%         xlim(UE_SINR_to_throughput_axes,[wideband_SINR_ecdf_all.min wideband_SINR_ecdf_all.max]);
%         ylim(UE_SINR_to_throughput_axes,[throughput_Mbps_ecdf_all.min throughput_Mbps_ecdf_all.max]);
%     else
%         xlim(UE_SINR_to_throughput_axes,'auto');
%         ylim(UE_SINR_to_throughput_axes,'auto');
%     end
%     hold(UE_SINR_to_throughput_axes,'off');
%     
%     % SINR-to-spectral efficiency
%     hold(UE_SINR_to_spectral_eff_axes,'all');
%     if plot_points_scatter
%         if plot_in_GUI
%             scatter(UE_SINR_to_spectral_eff_axes,wideband_SINR_vector,spectral_eff_vector,'.b');
%         else
%             if plot_mean_scatter
%                 scatter(UE_SINR_to_spectral_eff_axes,wideband_SINR_vector,spectral_eff_vector,'.b','SizeData',50);
%             else
%                 scatter(UE_SINR_to_spectral_eff_axes,wideband_SINR_vector,spectral_eff_vector,'.b','SizeData',150);
%             end
%         end
%     end
%     if plot_mean_scatter
%         if plot_in_GUI
%             scatter(UE_SINR_to_spectral_eff_axes,wideband_SINR_binned,spectral_eff_binned,'.r');
%         else
%             scatter(UE_SINR_to_spectral_eff_axes,wideband_SINR_binned,spectral_eff_binned,'.r','SizeData',350);
%         end
%     end
%     grid(UE_SINR_to_spectral_eff_axes,'on');
%     title(UE_SINR_to_spectral_eff_axes,'UE wideband SINR-to-spectral efficiency mapping');
%     xlabel(UE_SINR_to_spectral_eff_axes,'UE wideband SINR [dB]');
%     ylabel(UE_SINR_to_spectral_eff_axes,'average UE specctral efficiency [bit/cu]');
%     if constant_axes_limits
%         xlim(UE_SINR_to_spectral_eff_axes,[wideband_SINR_ecdf_all.min wideband_SINR_ecdf_all.max]);
%         ylim(UE_SINR_to_spectral_eff_axes,[spectral_eff_ecdf.min spectral_eff_ecdf.max]);
%     else
%         xlim(UE_SINR_to_spectral_eff_axes,'auto');
%         ylim(UE_SINR_to_spectral_eff_axes,'auto');
%     end
%     hold(UE_SINR_to_spectral_eff_axes,'off');
%     
%     cell_sum_throughput_non_NaN      = cell_sum_throughput(isfinite(cell_sum_throughput));
%     cell_sum_throughput_non_NaN_mean = mean(cell_sum_throughput_non_NaN);
%     ignored_cells = sum(isnan(cell_sum_throughput));
%     
%     sep = '---------------------------------------';
%     statistics_text = sprintf([
%         '%s\nSimulations statistics:\n\n',...
%         '%g cells, %g UEs\n',...
%         'Simulation length: %g TTIs\n',...
%         'Scheduler: %s\n',...
%         'Mode: %gx%g, %s\n',...
%         '%s\n',...
%         'Cell statistics:\n\n',...
%         'Fairness index: %g\n',...
%         'Peak/Avg/Edge UE throughput:\n',...
%         '%3.2f/%3.2f/%3.2f Mb/s\n',...
%         'Average cell throughput: %3.2fMb/s\n',...
%         'Ignored cells (disabled): %g\n',...
%         'mean RB occupancy: %3.2f%%\n%s',...
%         ],...
%         sep,...
%         length(cells_to_plot),...
%         sum(UEs_to_use),...
%         simulation_data.LTE_config.simulation_time_tti,...
%         simulation_data.LTE_config.scheduler,...
%         simulation_data.LTE_config.nTX,...
%         simulation_data.LTE_config.nRX,...
%         utils.miscUtils.tx_mode_to_string(simulation_data.LTE_config.tx_mode),...
%         sep,...
%         fairness_index,...
%         throughput_Mbps_ecdf.p95,...
%         throughput_Mbps_ecdf.mean_x,...
%         throughput_Mbps_ecdf.p05,...
%         cell_sum_throughput_non_NaN_mean,...
%         ignored_cells,...
%         mean_RB_occupancy_ratio*100,sep);
%     if plot_in_GUI
%         set(handles.cell_statistics_text,'String',statistics_text);
%     else
%         set(handles.cell_statistics_text,'String',statistics_text);
%         fprintf('%s\n',statistics_text);
%     end
% end

function [WideBand_SINR,UE_Throughput] = LTE_show_aggregate_results(simulation_data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

UE_list = cell(1,length(simulation_data.UEs));
for u_=1:length(simulation_data.UEs)
    UE_list{u_} = sprintf('%g',u_);
end

cell_list = cell(1,length(simulation_data.eNodeBs));
for c_=1:length(simulation_data.eNodeBs)
    cell_list{c_} = sprintf('%g',c_);
end

simulation_data.UE_list   = UE_list;
simulation_data.cell_list = cell_list;

if isfield(simulation_data.LTE_config,'default_shown_GUI_cells') && ~isempty(simulation_data.LTE_config.default_shown_GUI_cells)
    cells_to_plot = simulation_data.LTE_config.default_shown_GUI_cells;
    if simulation_data.LTE_config.compact_results_file
        cells_to_plot = cells_to_plot(cells_to_plot<=length(simulation_data.the_eNodeB_traces)); % Filtero out possible out of range values
    else
        cells_to_plot = cells_to_plot(cells_to_plot<=length(simulation_data.simulation_traces.eNodeB_tx_traces)); % Filtero out possible out of range values
    end
 end

% Plot main plot


% Load data
eNodeB_sites           = simulation_data.eNodeB_sites;
eNodeBs                = simulation_data.eNodeBs;
UEs                    = simulation_data.UEs;
networkPathlossMap     = simulation_data.networkPathlossMap;

if isfield(simulation_data,'simulation_traces')
    % non-compact results
    the_UE_traces     = [simulation_data.simulation_traces.UE_traces];
    the_eNodeB_traces = [simulation_data.simulation_traces.eNodeB_tx_traces];
else
    the_UE_traces     = simulation_data.the_UE_traces;
    the_eNodeB_traces = simulation_data.the_eNodeB_traces;
end
cells_to_plot         = cells_to_plot(cells_to_plot<=length(the_eNodeB_traces)); % Filtero out possible out of range values
N_UEs = length(UEs);

%     UE_throughput_axes           = handles.UE_throughput_axes;
%     UE_spectral_eff_axes         = handles.UE_spectral_eff_axes;
%     UE_wideband_SINR_axes        = handles.UE_wideband_SINR_axes;
%     UE_SINR_to_throughput_axes   = handles.UE_SINR_to_throughput_axes;
%     UE_SINR_to_spectral_eff_axes = handles.UE_SINR_to_spectral_eff_axes;


% Clear old plots
% cla(UE_throughput_axes);
% cla(UE_spectral_eff_axes);
% cla(UE_wideband_SINR_axes);
% cla(UE_SINR_to_throughput_axes);
% cla(UE_SINR_to_spectral_eff_axes);

if ~isempty(cells_to_plot)
    [UEs_to_use cell_sum_throughput] = utils.resultsFileReader.get_UEs_in_given_cells(cells_to_plot,the_UE_traces);
    
    the_UE_traces_to_plot = the_UE_traces(UEs_to_use);
    wideband_SINRs_all = reshape([the_UE_traces.wideband_SINR],simulation_data.LTE_config.simulation_time_tti,[]);
    wideband_SINRs_all = wideband_SINRs_all(1,:);
    
    % To get the average RB occupancy
    nRB                     = unique([the_eNodeB_traces(cells_to_plot).RB_grid_size]); % It should be a unique value!
    RB_occupancy_ratio      = double([the_eNodeB_traces(cells_to_plot).scheduled_RBs])/nRB;
    mean_RB_occupancy_ratio = mean(RB_occupancy_ratio);
    
    % To get the axes limits
    throughput_Mbps_ecdf_all = utils.miscUtils.ecdf([the_UE_traces.average_throughput_Mbps]);
    spectral_eff_ecdf_all    = utils.miscUtils.ecdf([the_UE_traces.average_spectral_efficiency_bit_per_cu]);
    wideband_SINR_ecdf_all   = utils.miscUtils.ecdf(wideband_SINRs_all);
        
    % The actual values
    throughput_Mbps_ecdf = utils.miscUtils.ecdf([the_UE_traces_to_plot.average_throughput_Mbps]);
    fairness_index       = sum(throughput_Mbps_ecdf.input_data).^2 / sum(throughput_Mbps_ecdf.input_data.^2) / sum(isfinite(throughput_Mbps_ecdf.input_data));
    spectral_eff_ecdf    = utils.miscUtils.ecdf([the_UE_traces_to_plot.average_spectral_efficiency_bit_per_cu]);
    wideband_SINR_ecdf   = utils.miscUtils.ecdf(wideband_SINRs_all(UEs_to_use));
    WideBand_SINR=wideband_SINRs_all(UEs_to_use);
    UE_Throughput=[the_UE_traces_to_plot.average_throughput_Mbps];
    UE_Efficiency=[the_UE_traces_to_plot.average_spectral_efficiency_bit_per_cu];
    % Wideband SINR Text file
    
%     fid = fopen('Rx SINR.txt','wt');
%     fprintf(fid,'%0.6f %0.6f\n',wideband_SINR_ecdf.x);
%     fclose(fid);
    
%     % Wideband SINR-to-throughput plot
%     wideband_SINR_vector                       = wideband_SINR_ecdf.input_data;
%     spectral_eff_vector                        = spectral_eff_ecdf.input_data;
%     throughput_vector                          = throughput_Mbps_ecdf.input_data;
%     [wideband_SINR_binned  throughput_binned   numel_bins] = utils.miscUtils.fit_scatterplot_data(wideband_SINR_vector,throughput_vector,50);
%     [wideband_SINR_binned2 spectral_eff_binned numel_bins] = utils.miscUtils.fit_scatterplot_data(wideband_SINR_vector,spectral_eff_vector,50);
%    
%           
%     hold(UE_SINR_to_throughput_axes,'all');
%     %plot(UE_SINR_to_throughput_axes,wideband_SINR_ecdf.x,throughput_Mbps_ecdf.x(1:length(wideband_SINR_ecdf.x),:),'k');
%     plot(UE_SINR_to_throughput_axes,wideband_SINR_ecdf.x(1:length(throughput_Mbps_ecdf.x),:),throughput_Mbps_ecdf.x,'k'); % For Macro Only case
%     grid(UE_SINR_to_throughput_axes,'on');
%     title(UE_SINR_to_throughput_axes,'UE Rx SINR-to-throughput mapping');
%     xlabel(UE_SINR_to_throughput_axes,'UE Rx SINR [dB]');
%     ylabel(UE_SINR_to_throughput_axes,'UE Average throughput [Mb/s]');
%     %plot(wideband_SINR_ecdf.x,throughput_Mbps_ecdf.x(1:length(wideband_SINR_ecdf.x),:),'k'); % To plot the 
%          
%     
%     % Throughput ECDF
%     hold(UE_throughput_axes,'all');
%     plot(UE_throughput_axes,throughput_Mbps_ecdf.x,throughput_Mbps_ecdf.f,'blue');
%     grid(UE_throughput_axes,'on');
%     title(UE_throughput_axes,'UE Average Throughput');
%     xlabel(UE_throughput_axes,' UE Average Throughput [Mb/s]');
%     ylabel(UE_throughput_axes,'CDF');
%    
%     xlim(UE_throughput_axes,'auto');
%    
%     hold(UE_throughput_axes,'off');
%     
%     % Spectral efficiency ECDF
%     hold(UE_spectral_eff_axes,'all');
%     plot(UE_spectral_eff_axes,spectral_eff_ecdf.x,spectral_eff_ecdf.f,'blue');
%     grid(UE_spectral_eff_axes,'on');
%     title(UE_spectral_eff_axes,'UE average spectral efficiency');
%     xlabel(UE_spectral_eff_axes,'UE average spectral efficiency [bit/cu]');
%     ylabel(UE_spectral_eff_axes,'CDF');
%     
%    
%         xlim(UE_spectral_eff_axes,'auto');
%     
%     hold(UE_spectral_eff_axes,'off');
%     
%     % Wideband SINR ECDF
%     hold(UE_wideband_SINR_axes,'all');
%     plot(UE_wideband_SINR_axes,wideband_SINR_ecdf.x,wideband_SINR_ecdf.f,'blue');
%     grid(UE_wideband_SINR_axes,'on');
%     title(UE_wideband_SINR_axes,'UE Rx SINR');
%     xlabel(UE_wideband_SINR_axes,'UE Rx SINR [dB]');
%     ylabel(UE_wideband_SINR_axes,'CDF');
%     
%         xlim(UE_wideband_SINR_axes,[wideband_SINR_ecdf_all.min wideband_SINR_ecdf_all.max]);
%     
%     
%     hold(UE_wideband_SINR_axes,'off');
%     
%     % SINR-to-throughput
%     hold(UE_SINR_to_throughput_axes_s,'all');
%     
   
end
