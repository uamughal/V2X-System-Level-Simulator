function [eNodeBs, eNodeBs_sectors] = LTE_add_signaling_power(LTE_config,eNodeBs,eNodeBs_sectors)
% Network generation. Either from a file (cache) or calling the necessary function.
% (c) Josep Colom Ikuno, INTHFT, 2009
% www.nt.tuwien.ac.at
    
%% Just Add signaling power
          
        % Add the power separation. X% to signaling/pilots (always on) and the rest for data
        set_signaling_power(LTE_config,eNodeBs);
        
        % Store the other eNodeBs as (potential) interferers
        for s_ = 1:length(eNodeBs_sectors)
            eNodeBs_sectors(s_).neighbors_eNodeB = eNodeBs_sectors([1:(s_-1) (s_+1):length(eNodeBs_sectors)]);
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

