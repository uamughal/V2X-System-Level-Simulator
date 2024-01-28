classdef downlinkChannelModel < handle
% Represents the downlink channel model that a specific user possesses. Each UE
% instance will have its own specific channel model.
% (c) Josep Colom Ikuno, INTHFT, 2008

   properties
       % These could actually be maps or an implementation that directly
       % calculates every time it is invoked.
       macroscopic_pathloss_model_is_set = false;
       macroscopic_pathloss_model
       shadow_fading_model_is_set = false;
       shadow_fading_model
       fast_fading_model_is_set = false;
       fast_fading_model
       
       % All eNodeBs
       eNodeBs
       % All UEs
       UEs
       % User to which this channel model is attached to
       attached_UE
       
       % This variables model the downlink signaling and the data that was
       % transmitted
       % RB_grid (retrieved via a function call)
       
   end

   methods
       % class constructor
       function obj = downlinkChannelModel(aUE)
           obj.attached_UE = aUE;
       end
       % Returns the macroscopic pathloss in dB between the given user's
       % position and his eNodeB sector. Returns 0 in case no model is specified.
       function [pathloss is_set] = macroscopic_pathloss(obj)
           if ~obj.macroscopic_pathloss_model_is_set
               pathloss = 0;
               is_set   = false;
               return
           else
               % Get eNodeB id and sector number
               eNodeB_id = obj.attached_UE.attached_eNodeB.eNodeB_id;
               pos = obj.attached_UE.pos;
               % Now get the pathloss
               pathloss = obj.macroscopic_pathloss_model.get_pathloss_eNodeB(pos,eNodeB_id);
               is_set   = true;
           end
       end
       
       % Returns the macroscopic pathloss in dB between the given user's
       % position and a given eNodeB.
       function [pathloss is_set] = interfering_macroscopic_pathloss(obj,interferingEnodeBids)
           if ~obj.macroscopic_pathloss_model_is_set
               pathloss = zeros([length(interferingEnodeBids) 1]); % Column vector
               is_set   = false;
               return
           else
               pos = obj.attached_UE.pos;
               % Now get the pathloss
               pathloss = reshape(obj.macroscopic_pathloss_model.get_pathloss_eNodeB(pos,interferingEnodeBids),[],1);
               is_set   = true;
           end
       end
       
       % Returns the shadow fading pathloss in dB between the given user's
       % position and his eNodeB sector. Returns 0 in case no model is specified (for
       % example, when Odyssey data would be used).
       function [pathloss is_set] = shadow_fading_pathloss(obj)
           if ~obj.shadow_fading_model_is_set
               pathloss = 0;
               is_set   = false;
               return
           else
               % Get eNodeB id and sector number
               if isa(obj.shadow_fading_model,'channel_gain_wrappers.shadowFadingMapClaussen')&&~obj.shadow_fading_model.oneMapPerSite
                   % Take the cell ID
                   eNodeB_id = obj.attached_UE.attached_eNodeB.eNodeB_id;
               else
                   % Default
                   eNodeB_id = obj.attached_UE.attached_site.id;
               end
               pos      = obj.attached_UE.pos;
               pathloss = obj.shadow_fading_model.get_pathloss(pos,eNodeB_id);
               is_set   = true;
           end
       end
       
       % Returns the shadow fading pathloss in dB between the given user's 
       % position and a given eNodeB. Returns 0 in case no model is specified (for
       % example, when Odyssey data would be used).
       function [pathloss is_set] = interfering_shadow_fading_pathloss(obj,interferingSiteIds,interfering_eNodeB_ids)
           if ~obj.shadow_fading_model_is_set
               pathloss = zeros([length(interferingSiteIds) 1]); % Column vector
               is_set   = false;
               return
           else
               % Get eNodeB id and sector number
               if isa(obj.shadow_fading_model,'channel_gain_wrappers.shadowFadingMapClaussen')&&~obj.shadow_fading_model.oneMapPerSite
                   % Take the cell ID
                   interferingIds = interfering_eNodeB_ids;
               else
                   % Default
                   interferingIds = interferingSiteIds;
               end
               pos      = obj.attached_UE.pos;
               pathloss = reshape(obj.shadow_fading_model.get_pathloss(pos,interferingIds),[],1);
               is_set   = true;
           end
       end
       
       % Returns the RB_grid so this UE can know what belongs to him
       function the_RB_grid = RB_grid(obj)
           if ~isempty(obj.attached_UE.attached_eNodeB)
               the_RB_grid = obj.attached_UE.attached_eNodeB.RB_grid;
           else
               the_RB_grid = [];
           end
       end
       % Set a macroscopic pathloss model
       function set_macroscopic_pathloss_model(obj,macroscopic_pathloss_model)
           obj.macroscopic_pathloss_model        = macroscopic_pathloss_model;
           obj.macroscopic_pathloss_model_is_set = true;
       end
       % Set a shadow fading model
       function set_shadow_fading_model(obj,shadow_fading_model)
           obj.shadow_fading_model        = shadow_fading_model;
           obj.shadow_fading_model_is_set = true;
       end
       % Set a fast fading model
       function set_fast_fading_model_model(obj,fast_fading_model)
           obj.fast_fading_model        = fast_fading_model;
           obj.fast_fading_model_is_set = true;
       end
   end
end 
