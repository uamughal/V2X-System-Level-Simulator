function UEs = Vehicle_connected_UEs(LTE_config, UEs)
% This function is to setting broadcast coommunication range
   for u=1:length(UEs)
     UEs(u).connected_UEs=network_elements.UE;
   end
   
   %%
       for u_=1:20
            UEs(u_).neighbors_UEs = UEs([1:(u_-1) (u_+1):20]); 
            j=1;
                 % Calculate the distance b/t Broadcast UE & Attached UE 
                 for i=1:19                 
                     UE_distance(i)=sqrt((UEs(u_).pos(1)-UEs(u_).neighbors_UEs(i).pos(1)).*(UEs(u_).pos(1)-UEs(u_).neighbors_UEs(i).pos(1))...
                             +(UEs(u_).pos(2)-UEs(u_).neighbors_UEs(i).pos(2)).*(UEs(u_).pos(2)-UEs(u_).neighbors_UEs(i).pos(2)));
                                
                 % calculate the relative Speed b/t Broadcast UE & Attached UE
                  if UEs(u_).walking_model.direction == UEs(u_).neighbors_UEs(i).walking_model.direction
                        UEs(u_).neighbors_UEs(i).relative_speed=0;
                  end
                  if UEs(u_).walking_model.direction == - UEs(u_).neighbors_UEs(i).walking_model.direction
                       UEs(u_).neighbors_UEs(i).relative_speed=2*LTE_config.UE_speed;
                  end

                  if UE_distance(i) < 150 
                      UEs(u_).connected_UEs(j)= UEs(u_).neighbors_UEs(i); 
                      %UEs(u_).connected_UEs(j).attached_UE=UEs(u_);
                      %UEs(u_).connected_UEs(j).UE_distance=UE_distance(i);
                      j=j+1;
                  end                
                end               
       end
        
   %%
        for u_=21:44
            UEs(u_).neighbors_UEs = UEs([21:(u_-1) (u_+1):44]); 
               j=1;
                 for i=1:23       
                     UE_distance1(i)=sqrt((UEs(u_).pos(1)-UEs(u_).neighbors_UEs(i).pos(1)).*(UEs(u_).pos(1)-UEs(u_).neighbors_UEs(i).pos(1))...
                             +(UEs(u_).pos(2)-UEs(u_).neighbors_UEs(i).pos(2)).*(UEs(u_).pos(2)-UEs(u_).neighbors_UEs(i).pos(2)));
                            
                         
                 % calculate the relative Speed b/t Broadcast UE & Attached UE
                  if UEs(u_).walking_model.direction == UEs(u_).neighbors_UEs(i).walking_model.direction
                        UEs(u_).neighbors_UEs(i).relative_speed=0;
                  end
                  if UEs(u_).walking_model.direction == - UEs(u_).neighbors_UEs(i).walking_model.direction
                       UEs(u_).neighbors_UEs(i).relative_speed=2*LTE_config.UE_speed;
                  end

                  if UE_distance1(i) < 150
                      UEs(u_).connected_UEs(j)= UEs(u_).neighbors_UEs(i); 
                      %UEs(u_).connected_UEs(j).attached_UE=UEs(u_);
                      %UEs(u_).connected_UEs(j).UE_distance=UE_distance(i);
                      j=j+1;
                  end                 
                end 
        end
        
        
   %% 
        for u_=45:76
            UEs(u_).neighbors_UEs = UEs([45:(u_-1) (u_+1):76]);  
            j=1;   
                 for i=1:31  
                     UE_distance2(i)=sqrt((UEs(u_).pos(1)-UEs(u_).neighbors_UEs(i).pos(1)).*(UEs(u_).pos(1)-UEs(u_).neighbors_UEs(i).pos(1))...
                             +(UEs(u_).pos(2)-UEs(u_).neighbors_UEs(i).pos(2)).*(UEs(u_).pos(2)-UEs(u_).neighbors_UEs(i).pos(2)));
                                
                 % calculate the relative Speed b/t Broadcast UE & Attached UE
                  if UEs(u_).walking_model.direction == UEs(u_).neighbors_UEs(i).walking_model.direction
                        UEs(u_).neighbors_UEs(i).relative_speed=0;
                  end
                  if UEs(u_).walking_model.direction == - UEs(u_).neighbors_UEs(i).walking_model.direction
                       UEs(u_).neighbors_UEs(i).relative_speed=2*LTE_config.UE_speed;
                  end
                  
                  % connected UEs for get broadcast message
                  if UE_distance2(i) < 150 
                      UEs(u_).connected_UEs(j)= UEs(u_).neighbors_UEs(i); 
                      %UEs(u_).connected_UEs(j).attached_UE=UEs(u_);
                      %UEs(u_).connected_UEs(j).UE_distance=UE_distance(i);
                       j=j+1;
                  end                  
                end 
        end
        
  %% 
        for u_=77:82
            UEs(u_).neighbors_UEs = UEs([77:(u_-1) (u_+1):82]);
             j=1;     
                 % Calculate the distance b/t Broadcast UE & Attached UE 
                 for i=1:5  
                     UE_distance3(i)=sqrt((UEs(u_).pos(1)-UEs(u_).neighbors_UEs(i).pos(1)).*(UEs(u_).pos(1)-UEs(u_).neighbors_UEs(i).pos(1))...
                             +(UEs(u_).pos(2)-UEs(u_).neighbors_UEs(i).pos(2)).*(UEs(u_).pos(2)-UEs(u_).neighbors_UEs(i).pos(2)));
                                
                 % calculate the relative Speed b/t Broadcast UE & Attached UE
                  if UEs(u_).walking_model.direction == UEs(u_).neighbors_UEs(i).walking_model.direction
                        UEs(u_).neighbors_UEs(i).relative_speed=0;
                  end
                  if UEs(u_).walking_model.direction == - UEs(u_).neighbors_UEs(i).walking_model.direction
                       UEs(u_).neighbors_UEs(i).relative_speed=2*LTE_config.UE_speed;
                  end
                  
                  % connected UEs for get broadcast message
                  if UE_distance3(i) < 150
                      UEs(u_).connected_UEs(j)= UEs(u_).neighbors_UEs(i); 
                      %UEs(u_).connected_UEs(j).attached_UE=UEs(u_);
                      %UEs(u_).connected_UEs(j).UE_distance=UE_distance(i);
                       j=j+1;
                  end                  
                end 
        end    
end