classdef PropFair_Traffic_CoMP < CoMP.schedulers.CoMP_scheduler
    % (c) Thomas Blazek, 2014,ITC,  thomas.blazek@nt.tuwien.ac.at
    % This is an example comp scheduler. It mirrors the behaviour of the
    % standard Round Robin scheduler (extraction from ringbuffer), but does
    % so simultaneously for all
    
    properties
        
        UE_matrix
        last_scheduled
    end
    
    methods
        function obj =  PropFair_Traffic_CoMP(site)
            obj = obj@CoMP.schedulers.CoMP_scheduler(site);
            obj.user_allocations = cell(size(site.attached_eNodeBs));
            
            obj.UE_matrix = [];
            obj.last_scheduled = ones(1, length(site.attached_eNodeBs));
        end
        
        
        
        function ues = extract_users(obj,nr_ues, eNodeB_idx)
            % Extract from ringbuffer
            mod_length = length(obj.UE_matrix{eNodeB_idx});
            idx_array = mod(((1:nr_ues)+obj.last_scheduled(eNodeB_idx)-1), mod_length)+1;
            ue_temp = obj.UE_matrix{eNodeB_idx}(idx_array);
            obj.last_scheduled(eNodeB_idx) = idx_array(end);
            ues = [ue_temp.id];
        end
        
        function schedule_users(obj, eNodeB,attached_UEs,last_received_feedbacks)
            eNodeB_array = obj.CoMP_site.attached_eNodeBs;
            for i_ = 1:length(eNodeB_array) % set the values
                sched(i_) = eNodeB_array(i_).scheduler;
                RB_grid(i_) = sched(i_).RB_grid;
                obj.user_allocations{i_}= zeros(size(RB_grid(i_).user_allocation));
                RB_grid(i_).user_allocation = zeros(size(RB_grid(i_).user_allocation));
                RB_grid(i_).size_bits = 0;
                tx_mode(i_)           = sched(i_).default_tx_mode;
                current_TTI(i_)       = sched(i_).clock.current_TTI;
                %% added by me
                obj.UE_matrix{i_}=eNodeB_array(i_).attached_UEs_vector;
            end
            
            site = obj.CoMP_site;
           
            for rb_ = 1:size(RB_grid(1).user_allocation,1)
                current_UEs = []; % array of UEs scheduled in same RB (in different eNodeBs)
                for i_ =randperm(length(eNodeB_array)) %randperm so no eNodeB is always first
                    forbidden_users = site.get_forbidden_UEs(eNodeB_array(i_));
                    % checks which users are not allowed alongside this
                    % eNodeB
                    %%I add to avoid that there is no UEs in some eNodeB
                 if ~isempty(eNodeB_array(i_).attached_UEs_vector)
                    if isempty(current_UEs) ||isempty(forbidden_users)||...
                            sum(sum(bsxfun(@eq, forbidden_users, current_UEs)))==0
                        
                        next_ue =  obj.extract_users(1, i_);
                        dont_schedule = 0;
                        if(~isempty(site.cooperation_rules))
                            % check which eNodeBs are not allowed along the
                            % chosen UE
                            forbidden_eNodeBs_idx = [[site.cooperation_rules.UE_id] == next_ue];
                            forbidden_eNodeBs = [site.cooperation_rules(forbidden_eNodeBs_idx).dont_schedule];
                            for j_ = 1:length(forbidden_eNodeBs)
                                if obj.user_allocations{ [eNodeB_array.eNodeB_id] == forbidden_eNodeBs(j_).eNodeB_id}~=0
                                    % if a forbidden UE already has a RB
                                    % assigned, dont schedule
                                    dont_schedule =1;
                                    break;
                                end
                                
                            end
                        end
                        if ~dont_schedule
                            % Schedule user
                            obj.user_allocations{i_}(rb_) = next_ue;
                            current_UEs = [current_UEs next_ue];
                        else
                            % Push the ringbuffer back, so the user is not
                            % left out on the next run
                            obj.last_scheduled(i_) = mod(obj.last_scheduled(i_)-2, length(obj.UE_matrix{i_}))+1;
                        end
                    end
                 end
                end
            end
                          
            
        end
        
     function RBs = Propfair_Traffic_scheduler(obj,N_UE,N_RB,c,user_ind,attached_UEs)
           % core scheduling function (same in LL and SL)
           
           if ~mod(obj.clock.current_TTI-1,5)
               overhead = obj.overhead_ref+obj.overhead_sync;
           else
               overhead = obj.overhead_ref;
           end
           alpha_temp   = 1;
           RBs          = zeros(N_RB*N_UE,1);
           RB_UEs       = zeros(N_RB,N_UE);
           true_RB_UEs  = RB_UEs;
           rb_vect      = zeros(N_UE,1)';
           bits_left    = zeros(N_UE,1);
           isbits       = false(N_UE,1);
           metric       = ones(N_RB,N_UE)*-Inf;
           RB_set     = true(N_RB,1);
           RB_UEs     = false(N_RB,N_UE);
           
           for ii = 1:N_UE         % Check if data is available
                if strcmp(attached_UEs(user_ind(ii)).traffic_model.type,'fullbuffer')
                    bits_left(user_ind(ii)) = 1;
                else
                    bits_left(user_ind(ii)) = attached_UEs(user_ind(ii)).traffic_model.bit_count;
                end
                isbits = logical(bits_left);
           end          
           
           % Precalculated values taken out from the loop (speeds up simulations)
           cleaned_c_log_matrix = log10(max(c,eps)*2*12*7);
           avgd_UE_throughputs  = (obj.av_const-1)*obj.av_throughput(user_ind);

           % Calculate metric for each RB and attached user
           
           %% Add for LTE-R UE
            ue_r_ind=find([attached_UEs.id]==841);
         
            prb_r=0; %assigned PRB for LTE-R
            if ~isempty(ue_r_ind)
                bits_left_r=attached_UEs(ue_r_ind).traffic_model.bit_count;
              for rr = 1:N_RB
               if sum(bits_left_r)
                   prb_r=prb_r+1;                
                   res                    = find(RB_set);
                   metric                 = -Inf(N_RB,N_UE);
                   UE_avgd_pre_metric     = -alpha_temp*log10(max(avgd_UE_throughputs+sum(RB_UEs.*c,1)*2*12*7,eps));
                   UE_avgd_pre_metric_mat = UE_avgd_pre_metric(ones(1,N_RB),:);
                   
                   metric(res(1:sum(RB_set)),:) = cleaned_c_log_matrix(res(1:sum(RB_set)),:)+UE_avgd_pre_metric_mat(res(1:sum(RB_set)),:);
                   % for u_ = 1:N_UE
                   % for r_ = 1:N_RB
                   % metric(res(r_),u_) = c(res(r_),u_)*12*7/((obj.av_const-1)*obj.av_throughput(user_ind(u_))+RB_UEs(:,u_).'*c(:,u_)*12*7);      % 12*7 equals the number of elements in a RB
                   % metric(res(r_),u_) = log10(max(c(res(r_),u_),eps)*12*7)-alpha_temp*log10(max((obj.av_const-1)*obj.av_throughput(user_ind(u_))+RB_UEs(:,u_).'*c(:,u_)*12*7,eps)); % Old implementation
                   % metric(res(r_),u_) = log10(max(c(res(r_),u_)*12*7,eps))-alpha_temp*log10(max(obj.av_throughput(user_ind(u_)),eps));
                   % end
                   % end                 
                                             
                   metric(:,~isbits(user_ind)) = -Inf;          % if there are no bits left, set metric to -Inf% only set the metirc -inf ii instead of users_ind(ii)
                   
                 %  maxi            = max(metric(:));
                   
                   UE_idx=find(user_ind==ue_r_ind);
                   metri=metric(:,UE_idx);
                   RB_idxs=find(metri==max(metri));
                   RB_indx=randi(length(RB_idxs));
                   RB_idx=RB_idxs(RB_indx);
                   
                   ind             = randi(length(RB_idx));
                   
                   tmp_UE          = UE_idx(ind);
                   tmp_RB          = RB_idx(ind);
                   
                   RB_set(tmp_RB)               = false;
                   RB_UEs(tmp_RB,tmp_UE)        = true;
                   
                   % coarse decrease for UE who got the current RB and check if there are still bits left
%                    if ~strcmp(attached_UEs(tmp_UE).traffic_model.type,'fullbuffer')
%                        
%                        if sum(RB_UEs(:,tmp_UE)) <= 1       %coarse decrease with crc-bits subtracted (only for non-fullbuffer)
%                            attached_UEs(tmp_UE).traffic_model.coarse_decrease(c(tmp_RB,tmp_UE)*(2*12*7-overhead-24));
%                        else                                     %crc is subtracted only once
%                            attached_UEs(tmp_UE).traffic_model.coarse_decrease(c(tmp_RB,tmp_UE)*(2*12*7-overhead));
%                        end
%                        bits_left(tmp_UE) = attached_UEs(tmp_UE).traffic_model.bit_count;
%                        isbits(tmp_UE)    = logical(bits_left(tmp_UE));
%                    end
                   if ~strcmp(attached_UEs(user_ind(tmp_UE)).traffic_model.type,'fullbuffer')
                       
                       if sum(RB_UEs(:,tmp_UE)) <= 1       %coarse decrease with crc-bits subtracted (only for non-fullbuffer)
                           attached_UEs(user_ind(tmp_UE)).traffic_model.coarse_decrease(c(tmp_RB,user_ind(tmp_UE))*(2*12*7-overhead-24));
                       else                                     %crc is subtracted only once
                           attached_UEs(user_ind(tmp_UE)).traffic_model.coarse_decrease(c(tmp_RB,user_ind(tmp_UE))*(2*12*7-overhead));
                       end
                       bits_left(user_ind(tmp_UE)) = attached_UEs(user_ind(tmp_UE)).traffic_model.bit_count;
                       isbits(user_ind(tmp_UE))    = logical(bits_left(user_ind(tmp_UE)));
                       bits_left_r= bits_left(user_ind(tmp_UE));
                   end
               end
              end
           
           end
           
            %%
                   
           
           
           
           for rr = 1:N_RB-prb_r
               if sum(bits_left)
                   
                   
                   res                    = find(RB_set);
                   metric                 = -Inf(N_RB,N_UE);
                   UE_avgd_pre_metric     = -alpha_temp*log10(max(avgd_UE_throughputs+sum(RB_UEs.*c,1)*2*12*7,eps));
                   UE_avgd_pre_metric_mat = UE_avgd_pre_metric(ones(1,N_RB),:);
                   
                   metric(res(1:sum(RB_set)),:) = cleaned_c_log_matrix(res(1:sum(RB_set)),:)+UE_avgd_pre_metric_mat(res(1:sum(RB_set)),:);
                   % for u_ = 1:N_UE
                   % for r_ = 1:N_RB
                   % metric(res(r_),u_) = c(res(r_),u_)*12*7/((obj.av_const-1)*obj.av_throughput(user_ind(u_))+RB_UEs(:,u_).'*c(:,u_)*12*7);      % 12*7 equals the number of elements in a RB
                   % metric(res(r_),u_) = log10(max(c(res(r_),u_),eps)*12*7)-alpha_temp*log10(max((obj.av_const-1)*obj.av_throughput(user_ind(u_))+RB_UEs(:,u_).'*c(:,u_)*12*7,eps)); % Old implementation
                   % metric(res(r_),u_) = log10(max(c(res(r_),u_)*12*7,eps))-alpha_temp*log10(max(obj.av_throughput(user_ind(u_)),eps));
                   % end
                   % end                 
                                             
                   metric(:,~isbits(user_ind)) = -Inf;          % if there are no bits left, set metric to -Inf% only set the metirc -inf ii instead of users_ind(ii)
                   
                   maxi            = max(metric(:));
                   [RB_idx UE_idx] = find(metric == maxi);
                   ind             = randi(length(RB_idx));
                   
                   tmp_UE          = UE_idx(ind);
                   tmp_RB          = RB_idx(ind);
                   
                   RB_set(tmp_RB)               = false;
                   RB_UEs(tmp_RB,tmp_UE)        = true;
                   
                   % coarse decrease for UE who got the current RB and check if there are still bits left
%                    if ~strcmp(attached_UEs(tmp_UE).traffic_model.type,'fullbuffer')
%                        
%                        if sum(RB_UEs(:,tmp_UE)) <= 1       %coarse decrease with crc-bits subtracted (only for non-fullbuffer)
%                            attached_UEs(tmp_UE).traffic_model.coarse_decrease(c(tmp_RB,tmp_UE)*(2*12*7-overhead-24));
%                        else                                     %crc is subtracted only once
%                            attached_UEs(tmp_UE).traffic_model.coarse_decrease(c(tmp_RB,tmp_UE)*(2*12*7-overhead));
%                        end
%                        bits_left(tmp_UE) = attached_UEs(tmp_UE).traffic_model.bit_count;
%                        isbits(tmp_UE)    = logical(bits_left(tmp_UE));
%                    end
                   if ~strcmp(attached_UEs(user_ind(tmp_UE)).traffic_model.type,'fullbuffer')
                       
                       if sum(RB_UEs(:,tmp_UE)) <= 1       %coarse decrease with crc-bits subtracted (only for non-fullbuffer)
                           attached_UEs(user_ind(tmp_UE)).traffic_model.coarse_decrease(c(tmp_RB,user_ind(tmp_UE))*(2*12*7-overhead-24));
                       else                                     %crc is subtracted only once
                           attached_UEs(user_ind(tmp_UE)).traffic_model.coarse_decrease(c(tmp_RB,user_ind(tmp_UE))*(2*12*7-overhead));
                       end
                       bits_left(user_ind(tmp_UE)) = attached_UEs(user_ind(tmp_UE)).traffic_model.bit_count;
                       isbits(user_ind(tmp_UE))    = logical(bits_left(user_ind(tmp_UE)));
                   end
               end
           end
            RB_UEs = RB_UEs';
            RBs = RB_UEs(:);           
     end
     end
  end 
        
        
        
        
        
 


