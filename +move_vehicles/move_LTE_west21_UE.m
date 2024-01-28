function move_LTE_west21_UE(LTE_config,UEs,networkPathlossMap,eNodeBs)

some_UE_out_of_ROI_this_TTI = false;
[ x_range, y_range ] = networkPathlossMap.valid_range;
ISD= LTE_config.inter_eNodeB_distance;
lane_width= LTE_config.lane_width;
 % Road Line range for Vehicle UEs, There are different kinds of road line.
y_range=[5*lane_width,ISD/sqrt(3)-5*lane_width];  % [17.5m,270.5m] (road line distance=253m )
 
 % If a UE went outside of ROI, relocate him from the begaining of the train line. 
% distance_resolution_y=(y_range(2)-y_range(1))/(LTE_config.simulation_time_tti/LTE_config.TTI_resolution);


%% 5+5 vehicles are moving to south direction 
for u_ =63:67
    UEs(u_).move;        
    ROI_teleport       = ~UEs(u_).is_in_roi(x_range,y_range);
    handover_requested = UEs(u_).cell_change.requested;
    
    if ROI_teleport || handover_requested 
        old_eNodeB_id = UEs(u_).attached_eNodeB.eNodeB_id;
        if ROI_teleport
            
            x_init=-ISD/3-lane_width/2; 
            new_UE_position = [x_init, y_range(2)-1];%when vehicles go out, come in from opposite direction
            
            % Actually it should not be done like this. Measure all the neighboring cells' SINR and then decide which one is better
            [new_site_id, new_sector_id, new_eNodeB_id] = networkPathlossMap.cell_assignment(new_UE_position,'LTE_V');
            
            % Teleport UE
            UEs(u_).pos = new_UE_position;       
            
        elseif handover_requested
            new_eNodeB_id                     = UEs(u_).cell_change.target_eNodeB;
            UEs(u_).cell_change.requested     = false; % reset the handover request field
            UEs(u_).cell_change.target_eNodeB = [];
        end
        
        % Deattach UE from old eNodeB and reattach to new one
        UEs(u_).start_handover(eNodeBs(new_eNodeB_id));
        
        % Accordingly activate or deactivate UE
        % Check whether this UE should be deativated to speed-up
        % simulation. Allows for run-time changes.
        if ~isempty(LTE_config.compute_only_UEs_from_this_eNodeBs)
            if isempty(find(UEs(u_).attached_eNodeB.eNodeB_id==LTE_config.compute_only_UEs_from_this_eNodeBs,1))
                % Deactivate UE
                UEs(u_).deactivate_UE = true;
            else
                % Activate UE (already activated by default, but just in case)
                UEs(u_).deactivate_UE = false;
            end
        end
        
        % Print some debug
        if ~some_UE_out_of_ROI_this_TTI
            if LTE_config.debug_level>=1
                fprintf(1,'\n');
            end
            some_UE_out_of_ROI_this_TTI = true;
        end
        if LTE_config.debug_level>=1
            if ROI_teleport
                fprintf(' UE %g going out of ROI, teleporting to %g %g. eNodeB %g -> eNodeB %g\n',UEs(u_).id,new_UE_position(1),new_UE_position(2),old_eNodeB_id,new_eNodeB_id);
            elseif handover_requested
                fprintf(' UE %g handover request. eNodeB %g -> eNodeB %g\n',UEs(u_).id,old_eNodeB_id,new_eNodeB_id);
            end
        end
    end
end

for u_ =70:74
    UEs(u_).move;        
    ROI_teleport       = ~UEs(u_).is_in_roi(x_range,y_range);
    handover_requested = UEs(u_).cell_change.requested;
    
    if ROI_teleport || handover_requested 
        old_eNodeB_id = UEs(u_).attached_eNodeB.eNodeB_id;
        if ROI_teleport
            
            x_init=-ISD/3-lane_width*3/2; 
            new_UE_position = [x_init, y_range(2)-1];%when vehicles go out, come in from opposite direction
            
            % Actually it should not be done like this. Measure all the neighboring cells' SINR and then decide which one is better
            [new_site_id, new_sector_id, new_eNodeB_id] = networkPathlossMap.cell_assignment(new_UE_position,'LTE_V');
            
            % Teleport UE
            UEs(u_).pos = new_UE_position;
            
        elseif handover_requested
            new_eNodeB_id                     = UEs(u_).cell_change.target_eNodeB;
            UEs(u_).cell_change.requested     = false; % reset the handover request field
            UEs(u_).cell_change.target_eNodeB = [];
        end
        
        % Deattach UE from old eNodeB and reattach to new one
        UEs(u_).start_handover(eNodeBs(new_eNodeB_id));
        
        % Accordingly activate or deactivate UE
        % Check whether this UE should be deativated to speed-up
        % simulation. Allows for run-time changes.
        if ~isempty(LTE_config.compute_only_UEs_from_this_eNodeBs)
            if isempty(find(UEs(u_).attached_eNodeB.eNodeB_id==LTE_config.compute_only_UEs_from_this_eNodeBs,1))
                % Deactivate UE
                UEs(u_).deactivate_UE = true;
            else
                % Activate UE (already activated by default, but just in case)
                UEs(u_).deactivate_UE = false;
            end
        end
        
        % Print some debug
        if ~some_UE_out_of_ROI_this_TTI
            if LTE_config.debug_level>=1
                fprintf(1,'\n');
            end
            some_UE_out_of_ROI_this_TTI = true;
        end
        if LTE_config.debug_level>=1
            if ROI_teleport
                fprintf(' UE %g going out of ROI, teleporting to %g %g. eNodeB %g -> eNodeB %g\n',UEs(u_).id,new_UE_position(1),new_UE_position(2),old_eNodeB_id,new_eNodeB_id);
            elseif handover_requested
                fprintf(' UE %g handover request. eNodeB %g -> eNodeB %g\n',UEs(u_).id,old_eNodeB_id,new_eNodeB_id);
            end
        end
    end
end


if some_UE_out_of_ROI_this_TTI
    if LTE_config.debug_level>=1
        fprintf('              ');
    end
end
end

%{
function new_ue_pos=new_UE_pos_on_road(LTE_config,UE)

road_resolution=LTE_config.inter_eNodeB_distance/10;
ang=LTE_config.UE_direction;
new_ue_pos=UE.pos+road_resolution*[cosd(ang),sind(ang)];
end
%}





