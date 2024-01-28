clc
clear all
simulation_type = 'LTE_V';

%% Load parameters
LTE_config = LTE_load_params(simulation_type); % Choose simulation type
LTE_config  = LTE_load_params_dependant(LTE_config); % Load many many parameters
DEBUG_LEVEL = LTE_config.debug_level;  % we can find from LTE_V_.m
varargin{4} = []; % To avoid having an empty varargin variable

%% Parameter Initiallized & Added
LTE_config.traffic_models.usetraffic_model = true;
LTE_config.TTI_resolution                = 10; 
LTE_config.inter_eNodeB_distance         = 500; 
LTE_config.move_LTE_V_UEs                = true;  % give necessary parameter for this fuction
LTE_config.UE_direction                  = 90;
LTE_config.vehicle_UEs                   = 82; 
LTE_config.lane_width                    = 3.5;
LTE_config.CoMP_configuration1            = 1;
LTE_config.SIR_threshold              = 5;
LTE_config.compute_only_UEs_from_this_eNodeBs =   [13:15]; 
LTE_config.default_shown_GUI_cells =           [13:15];
%% Load and plot BLER curves
[ BLER_curves, CQI_mapper ] = LTE_init_load_BLER_curves(LTE_config);

%% Get the eNodeBs, the Macroscopic Pathloss and the Shadow Fading
[eNodeB_sites, eNodeBs, networkPathlossMap, networkShadowFadingMap] = LTE_init_network_generation(LTE_config,varargin);
varargin{2} = []; % clear varargin{2}

%% Add a clock to each network element
networkClock = network_elements.clock(LTE_config.TTI_length);
for b_=1:length(eNodeB_sites)
    eNodeB_sites(b_).clock = networkClock;
end


%% Create users (UEs) and add schedulers
% Vehicle only as UEs
[UEs, extra_UE_cache_info] = LTE_init_generate_users_and_add_schedulers(LTE_config,eNodeB_sites,eNodeBs,networkPathlossMap,CQI_mapper,BLER_curves,networkClock);

%% CoMP site Initialization
if isfield(LTE_config,'CoMP_configuration')
CoMP.initialize_CoMP_sites(LTE_config,eNodeBs);
end


%% Plot network
%    utils.plotUtils.plot_sector_SINR_cdfs(LTE_config.plots.sector_SINR_cdf,networkPathlossMap,networkShadowFadingMap,LTE_config,eNodeBs);
 %   LTE_plot_loaded_network(LTE_config,eNodeB_sites,eNodeBs,UEs,networkPathlossMap,CQI_mapper,true,networkShadowFadingMap);
    
%% Generate the fast fading traces for V2N
 pregenerated_ff = LTE_init_get_microscale_fading_SL_trace(LTE_config);
 
theta=pregenerated_ff.traces{1,1}.trace.theta;
zeta=pregenerated_ff.traces{1,1}.trace.zeta;
psi=pregenerated_ff.traces{1,1}.trace.psi;


 avg_theta_ff = mean(mean(pregenerated_ff.traces{1,1}.trace.theta))
 avg_psi_ff = mean(mean(pregenerated_ff.traces{1,1}.trace.psi))
 avg_zeta_ff = mean(mean(pregenerated_ff.traces{1,1}.trace.zeta))
                
%% Generate the fast fading traces for V2V 
   LTE_config.trace_params = channel_gain_wrappers.winnerChannelFactory1.get_default_config(LTE_config.frequency1,LTE_config.nTX,LTE_config.nRX,LTE_config.UE_speed);
   pregenerated_ff1 = LTE_init_get_microscale_fading_SL_trace1(LTE_config);
                                
%% SINR average algorithm
        switch LTE_config.SINR_averaging.algorithm
        case 'MIESM'
            the_SINR_averager = utils.miesmAveragerFast(LTE_config,LTE_config.SINR_averaging.BICM_capacity_tables,LTE_config.SINR_averaging.betas);
        otherwise
            error('SINR averaging algorithm %s not supported',LTE_config.SINR_averaging.algorithm);
        end
         
 %% UE initialization for eNB uplink & downlink
for u_=1:length(UEs)        
        % Add downlink channel (includes macroscopic pathloss, shadow fading and fast fading models)
        UEs(u_).downlink_channel = channel_models.downlinkChannelModel(UEs(u_));
        
        % Add eNodebs to the downlink channel (may be needed)
        UEs(u_).downlink_channel.eNodeBs = eNodeBs;
        
        % Set fast fading from the eNodeB to an attached UE.
        UEs(u_).downlink_channel.set_fast_fading_model_model(channel_gain_wrappers.fastFadingWrapper(pregenerated_ff,'random',length(eNodeBs)));
        
        % Macroscopic pathloss
        UEs(u_).downlink_channel.set_macroscopic_pathloss_model(networkPathlossMap);
        
        % Shadow fading (data obtained from planning tools already have this information incorporated)
        if LTE_config.macroscopic_pathloss_is_model
            UEs(u_).downlink_channel.set_shadow_fading_model(networkShadowFadingMap);
        end
        
        % Cache the RB_grid object in the UE object, so as to avoid too many calls to the function. This will have to be taken into account when implementing handover
        UEs(u_).RB_grid = UEs(u_).downlink_channel.RB_grid;
        
        % Uplink channel
        UEs(u_).uplink_channel = channel_models.uplinkChannelModel(...
            UEs(u_),...
            LTE_config.N_RB,...
            LTE_config.maxStreams,...
            LTE_config.feedback_channel_delay);
    
        % Set UE SINR averaging algorithm
        UEs(u_).SINR_averager = the_SINR_averager;
end
    
%% UE initialization for Sidelink Channel
 for u_=1:length(UEs)        
        % Add downlink channel (includes macroscopic pathloss, shadow fading and fast fading models)
        UEs(u_).side_downlink_channel = channel_models.side_downlinkChannelModel(UEs(u_));
        
        % Add eNodebs to the downlink channel (may be needed)
        UEs(u_).side_downlink_channel.eNodeBs = eNodeBs;
        
        % Set fast fading from the eNodeB to an attached UE.
        UEs(u_).side_downlink_channel.set_fast_fading_model_model(channel_gain_wrappers.fastFadingWrapper(pregenerated_ff1,'random',length(eNodeBs)));
                
        % Shadow fading (data obtained from planning tools already have this information incorporated)
        if LTE_config.macroscopic_pathloss_is_model
            UEs(u_).side_downlink_channel.set_shadow_fading_model(networkShadowFadingMap);
        end

        % Uplink channel
        UEs(u_).side_uplink_channel = channel_models.side_uplinkChannelModel(...
            UEs(u_),...
            LTE_config.N_RB,...
            LTE_config.maxStreams,...
            LTE_config.feedback_channel_delay);        
 end
 
  
 %% Initialise the tracing for V2V
    
    % Global traces for V2V
    simulation_traces = tracing.simTraces;
    simulation_traces.TTIs_to_ignore_when_calculating_aggregates = LTE_config.feedback_channel_delay;       

   % Traces from received UE feedbacks (eNodeB side)
    simulation_traces.eNodeB_rx_feedback_traces = tracing.receivedFeedbackTrace(...
        LTE_config.simulation_time_tti,...
        length(UEs),...
        LTE_config.N_RB,...
        LTE_config.maxStreams,...
        LTE_config.unquantized_CQI_feedback);
    simulation_traces.eNodeB_rx_feedback_traces.parent_results_object = simulation_traces;    
    
    % eNodeB traces
    for c_=1:length(eNodeBs)
        if c_==1
            simulation_traces.eNodeB_tx_traces = tracing.enodebTrace(eNodeBs(c_).RB_grid,LTE_config.maxStreams,LTE_config.simulation_time_tti,networkClock);
        else
            simulation_traces.eNodeB_tx_traces(c_) = tracing.enodebTrace(eNodeBs(c_).RB_grid,LTE_config.maxStreams,LTE_config.simulation_time_tti,networkClock);
        end
        simulation_traces.eNodeB_tx_traces(c_).parent_results_object = simulation_traces;
        eNodeBs(c_).sector_trace   = simulation_traces.eNodeB_tx_traces(c_);
        eNodeBs(c_).feedback_trace = simulation_traces.eNodeB_rx_feedback_traces;
    end
   
    % UE traces
    for u_=1:length(UEs)
        UEs(u_).trace = tracing.ueTrace(...
            LTE_config.simulation_time_tti,...
            LTE_config.N_RB,...
            LTE_config.nTX,...
            LTE_config.nRX,...
            LTE_config.maxStreams,...
            LTE_config.unquantized_CQI_feedback,...
            LTE_config.trace_SINR,...
            LTE_config.scheduler_params.av_window,...
            LTE_config.TTI_length,...
            LTE_config.reduced_feedback_logs);
        
        UEs(u_).trace.parent_results_object = simulation_traces;
        
        if u_==1
            simulation_traces.UE_traces     = UEs(u_).trace;
        else
            simulation_traces.UE_traces(u_) = UEs(u_).trace;
        end
    end
            
    LTE_config.UE_count     = length(UEs);
    LTE_config.site_count   = length(eNodeB_sites);
    LTE_config.eNodeB_count = length(eNodeBs);
    
%% Add Center UE eNodeB Traces from received UE feedbacks (eNodeB side)
  for u_ = 1:length(UEs)   
           simulation_traces.UE_eNodeB_rx_feedback_traces = tracing.receivedFeedbackTrace(...
                LTE_config.simulation_time_tti,...
                length(UEs),...
                LTE_config.N_RB,...
                LTE_config.maxStreams,...
                LTE_config.unquantized_CQI_feedback);                        
                
            % Initialize the UE eNodeB trace
            % simulation_traces.eNodeB_tx_traces= tracing.enodebTrace(UEs(u_).attached_eNodeB.RB_grid,LTE_config.maxStreams,LTE_config.simulation_time_tti,networkClock);
            simulation_traces.UE_eNodeB_tx_traces = tracing.enodebTrace(UEs(u_).RB_grid,LTE_config.maxStreams,LTE_config.simulation_time_tti,networkClock);
            UEs(u_).sectors(1).sector_trace   = simulation_traces.UE_eNodeB_tx_traces;
            UEs(u_).sectors(1).feedback_trace = simulation_traces.UE_eNodeB_rx_feedback_traces;
            simulation_traces.UE_eNodeB_rx_feedback_traces.parent_results_object = simulation_traces;
           % simulation_traces.UE_eNodeB_tx_traces(u_).parent_results_object = simulation_traces;
            
           % Add the network clock to UE eNodeB site
            UEs(u_).sectors(1).parent_eNodeB.clock=UEs(u_).attached_eNodeB.parent_eNodeB.clock;         
 end                   
        
    %% Give the schedulers access to the UE traces
    % Then they can make decisions base on their received throughput. More
    % complex and realistic solutions may be possible, but then the eNodeB
    % should dynamically allocate resources to store UE-related data (and then
    % release once the UE is not attached to it anymore). It is easier like this :P
    for c_=1:length(eNodeBs)
        eNodeBs(c_).scheduler.set_UE_traces(simulation_traces.UE_traces);
    end

    
 %% Main simulation loop
    if DEBUG_LEVEL>=1
        fprintf('Entering main simulation loop, %5.0f TTIs\n',LTE_config.simulation_time_tti);
    end 
    % Inititialize timer
    ticID_start_sim_loop = tic;
    starting_time = toc(ticID_start_sim_loop);
    
    num_markers = 5;
    s_markings  = round(linspace(1,length(eNodeBs),num_markers));
    u_markings  = round(linspace(1,length(UEs),num_markers));
    
    % Network clock is initialised to 0
    while networkClock.current_TTI < LTE_config.simulation_time_tti
        % First of all, advance the network clock
        networkClock.advance_1_TTI;
        if DEBUG_LEVEL>=1
            fprintf(' TTI %5.0f/%d: ',networkClock.current_TTI,LTE_config.simulation_time_tti);
        end

  
   %% Move Vehicle UEs
      %  move_vehicles.move_LTE_east_UE(LTE_config,UEs,networkPathlossMap,eNodeBs);
      %  move_vehicles.move_LTE_north_UE(LTE_config,UEs,networkPathlossMap,eNodeBs);
      %  move_vehicles.move_LTE_south_UE(LTE_config,UEs,networkPathlossMap,eNodeBs);
      %  move_vehicles.move_LTE_west11_UE(LTE_config,UEs,networkPathlossMap,eNodeBs);
      %  move_vehicles.move_LTE_west12_UE(LTE_config,UEs,networkPathlossMap,eNodeBs);
      %  move_vehicles.move_LTE_west21_UE(LTE_config,UEs,networkPathlossMap,eNodeBs);
      %  move_vehicles.move_LTE_west22_UE(LTE_config,UEs,networkPathlossMap,eNodeBs);
      %  move_vehicles.move_LTE_west23_UE(LTE_config,UEs,networkPathlossMap,eNodeBs);  
        
        
   %% Set Broadcast Range
    % Get the connected UEs of current UE according to the lanes deployment
    UEs = Vehicle_connected_UEs(LTE_config,UEs);


    %%  Scheduling & Feedback for V2V broadcast   
        if LTE_config.trace_simulation_mode
            attach_UEs_to_eNodeBs_according_to_trace(LTE_config,UEs,networkPathlossMap,eNodeBs);
        end   
        
        % The eNodeBs receive the feedbacks from the UEs
        for s_ = 1:length(eNodeBs)
               % Receives and stores the received feedbacks from the UEs
               eNodeBs(s_).receive_UE_feedback;              
               % Schedule users
               eNodeBs(s_).schedule_users;            
              if ~isempty(find(s_==s_markings,1))
                  if DEBUG_LEVEL>=1
                     fprintf('*');
                  end
              end
        end
        
    for u_=1:length(UEs)
        UEs(u_).deactivate_UE  =0;
    end
         
      % Call link quality model(non-0 delay case)
      % Call link performance model (evaluates whether TBs are received
      % corretly according to the information conveyed by the link quality
      % model. Additionally, send feedback (channel quality indicator +
      % ACK/NACK)       
                          
         if LTE_config.feedback_channel_delay~=0
           % Broadcast_UEs_order = randperm(size(UEs,2)); % permute UEs order
            for u_ = 1:length(UEs)         
            % Assume broadcasted UEs attached to the center UE eNodeB
             if UEs(u_).traffic_model.bit_count > 0 && ~isempty(UEs(u_).connected_UEs(1).id)
               UEs(u_).sectors(1).connected_UEs=network_elements.UE;       
                for i=1:length(UEs(u_).connected_UEs)
                    UEs(u_).sectors(1).connected_UEs(i)=UEs(u_).connected_UEs(i);
                end  
                                             
            % Center UE Receives and stores the received feedbacks from the broadcast range UEs

                  UEs(u_).sectors(1).receive_UE_broadcast_feedback; 

                % Measure SINR and prepare CQI feedback
                if ~UEs(u_).deactivate_UE                 
                    for i=1:length(UEs(u_).connected_UEs)
                        if ~UEs(u_).connected_UEs(i).deactivate_UE
                        % The UE(u_) as the UE eNodeB
                        UEs(u_).connected_UEs(i).attached_UE=UEs(u_);
                        % Call link quality model for the 1 broadcast UE
                        UEs(u_).connected_UEs(i).link_quality_model(LTE_config);
                        % Call link performance model for the 1 broadcast UE
                        UEs(u_).connected_UEs(i).link_performance_model;
                        % Send Channel feedback (channel quality indicator + ACK/NACK)
                        UEs(u_).connected_UEs(i).send_broadcast_feedback;
                        % deactivate the receivers
                        UEs(u_).connected_UEs(i).deactivate_UE = 1;
                        end
                    end  

                else
                    UEs(u_).dummy_link_quality_model(LTE_config);
                    UEs(u_).dummy_link_performance_model;
                    UEs(u_).send_feedback;
                end
                
                if ~isempty(find(u_==u_markings,1))
                    if DEBUG_LEVEL>=1
                        fprintf('+');
                    end
                end 
             end
           end
         end
                          

%%
        if DEBUG_LEVEL>=1
            fprintf('\n');
        end       
        if mod(networkClock.current_TTI,10)==0
            elapsed_time = toc(ticID_start_sim_loop);
            time_per_iteration = elapsed_time / networkClock.current_TTI;
            estimated_time_to_finish = (LTE_config.simulation_time_tti - networkClock.current_TTI)*time_per_iteration;
            estimated_time_to_finish_h = floor(estimated_time_to_finish/3600);
            estimated_time_to_finish_m = estimated_time_to_finish/60 - estimated_time_to_finish_h*60;
            fprintf('Time to finish: %3.0f hours and %3.2f minutes\n',estimated_time_to_finish_h,estimated_time_to_finish_m);
        end
     
    end
    
%%   Performance Matrics 
    %[prr Av_PRR0 Av_PRR]= PRR_calculation(UEs); %  Av_PRR0 Av_prr
    % prr = PRR_calculation(UEs); % 
    
%% 
    
    if ~isempty(extra_UE_cache_info)
        simulation_traces.extra_UE_info = extra_UE_cache_info;
    end
    
    if DEBUG_LEVEL>=1
        fprintf('Saving results to %s\n',LTE_config.results_file);
    end
    
    % Some post-processing
    simulation_traces.calculate_UE_aggregates;

        % Some options to save space
        if LTE_config.delete_ff_trace_at_end
            pregenerated_ff = [];
        end
        
        if LTE_config.delete_pathloss_at_end
            networkPathlossMap.pathloss = [];
            if exist('networkShadowFadingMap','var')
                networkPathlossMap.pathloss = [];
            else
                % Do nothing
            end
        end
        
%% plot simulation results


%%
outage=load('Outage probability.txt');
outage1=load('Outage probability1.txt');
    figure;
    plot(10:20:250,outage);
    hold on;
    plot(10:20:250,outage1);
    hold on;



%%
    Rx_SNR=load('Rx SNR.txt'); 
    Rx_SNR1=load('Rx SNR1.txt'); 
    Rx_SNR2=load('Rx SNR2.txt');
    Rx_SNR3=load('Rx SNR3.txt');
    
    Rx_SNR1=Rx_SNR(find(Rx_SNR>-8.378)); % for 150m
    Rx_SNR2=Rx_SNR(find(Rx_SNR> 1.356)); %for 50m 

    
    Rx_SNR4=Rx_SNR(find(Rx_SNR> -7.287));
    Rx_SNR5= Rx_SNR(find(Rx_SNR> 2.735));
    
    
    figure;
    cdfplot(Rx_SNR);
    hold on;
    cdfplot(Rx_SNR3);
    hold on;
    cdfplot(Rx_SNR1);
    hold on;
    cdfplot(Rx_SNR4);
    hold on;
    cdfplot(Rx_SNR2);
    hold on;
    cdfplot(Rx_SNR5);
    hold on;
    

    prr=load('Packet_RR.txt');
    prr1=load('Packet_RR1.txt');
    figure;
    plot(10:20:250,prr);
    hold on;
    plot(10:20:250,prr1);
    hold on;
    
    open ('Interference.fig');
    hold on;
    open ('Outage Probability.fig');
    hold on;
    open ('UE Throughput.fig');
    hold on;
    open ('Spectral Efficiency.fig');
    hold on;
        

    
    
%%
    for i=1:length(UEs)
        if UEs(i).trace.average_throughput_Mbps > 0
          Average_throughput{i}=UEs(i).trace.average_throughput_Mbps;
        else
          Average_throughput{i}=[];  
        end   
        if ~isempty(UEs(i).wideband_SINR)
           UE_RxSINR{i}=UEs(i).wideband_SINR;
        else
           UE_RxSINR{i}=[];          
        end
    end

     Average_throughput = cell2mat(Average_throughput);        
     fid = fopen('UE Throughput1.txt','at');
     fprintf(fid,'%0.6f \n',Average_throughput);
     fclose(fid);
    
    UE_RxSINR = cell2mat(UE_RxSINR);
    ind1=find(UE_RxSINR(:)< 23 );
    UE_RxSINR = UE_RxSINR(ind1);
    ind2=find(UE_RxSINR(:)> -20);
    UE_RxSINR = UE_RxSINR(ind2);

               
    
%%                 
    out_p= load('outage probability.txt');
    figure;
    plot(10:20:250,out_p);
    hold on;
    


%%     
    fid = fopen('Rx SNR.txt','at');
    fprintf(fid,'%0.6f \n',UE_RxSINR);
    fclose(fid);
            
    UE_throughput=load('UE Throughput.txt');
    Rx_SNR=load('Rx SNR.txt'); % for 250m
    PRR= load('Packet_RR.txt');
    
figure;
plot(10:20:250,PRR);
hold on;

figure;
  cdfplot(UE_throughput);
  hold on;  
  title('UE Average Throughput')
  hold on;
  
% outage probability 
Rx_SNR1=Rx_SNR(find(Rx_SNR>-8.378)); % for 150m
Rx_SNR2=Rx_SNR(find(Rx_SNR> 1.356)); %for 50m 


figure;
  cdfplot(Rx_SNR);
  hold on;
  cdfplot(Rx_SNR1);
  hold on;
  cdfplot(Rx_SNR2);
  hold on;
  title('UE Rx SNR');
  hold on;
  
  
snr=load('SINR_without CoMP.txt');
figure;
plot(1:length(snr),snr);

snr1=load('SINR_with CoMP.txt');
figure;
plot(1:length(snr1),snr1);


%
fid = fopen('Rx SINR.txt','wt');
fprintf(fid,'%0.6f \n',UE_RxSINR{:});
fclose(fid);   
    
fid = fopen('UE Throughput.txt','wt');
fprintf(fid,'%0.6f \n',[Average_throughput{:}]);
fclose(fid);


RxSINR=load('Rx SINR.txt'); 

figure;
 UE_Rx_SINR=utils.miscUtils.ecdf(RxSINR);
 title('UE Rx SINR');
 cdfplot(RxSINR);
 title('UE Rx SINR');
 hold on;
 plot(UE_Rx_SINR.x,UE_Rx_SINR.f,'red'); 







