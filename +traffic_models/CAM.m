classdef CAM < traffic_models.generic_tm
% This class is used for video streaming simulations
% Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2010 by INTHFT 
properties
    type = 'CAM';
    packet_max = 300*8;    % max slice size in bits. 1 byte ->8 bits
    packet_min = 190*8;    % min slice size
    inter_time = 100;     % mean slice interarrival time (encoder delay),message generation period 100 ms
    c = 0.6;  % data according to RAN R1-070674
    d = 0.4;
    state ;
    iit_offset = 1; %randi(10)-1; 
    packet_counter;
    delay_constraint = 500; % in ms, mean arrival time for all slices in a frame is 500ms, interarrival between frames is 500ms
end

methods
    function obj = CAM(UE,HARQ_delay)
        obj = obj@traffic_models.generic_tm(UE,HARQ_delay,'CAM');
        % obj.iit_offset = randi(100)-1; 
        obj.is_fullbuffer = false;
    end
    
    % Generate CAM Traffic Model  300 bytes-190-190-190-190bytes
    % Interval is 100ms
    function check_TTI(obj)
        if ~mod(obj.UE.clock.current_TTI-1-obj.iit_offset,500) && obj.UE.clock.current_TTI > obj.iit_offset
                obj.generate_packet(obj.packet_max,obj.type);
                 obj.packet_counter = 1;
        end         
        if obj.UE.clock.current_TTI > obj.iit_offset +1
            if ~mod(obj.UE.clock.current_TTI-1-obj.iit_offset,100) && obj.packet_counter <=5
                obj.generate_packet(obj.packet_min,obj.type);
                obj.packet_counter =  obj.packet_counter+1;
            end     
        end
        obj.bit_count = obj.get_buffer_length;
    end
   
    
        
   function rate = get_arrival_rate(obj)
        rate = obj.packet_mean*8*5/0.5; % first 8: byte->bit second 5: 5 slices per frame
   end
    
   function packet_parts = decrease_packets(obj,N_data_bits)
        packet_parts = [];
        for bb = 1:length(obj.packet_buffer)
            read_temp = obj.read_start+bb-1;
            if read_temp > length(obj.packet_buffer)
                read_temp = mod(read_temp,length(obj.packet_buffer));
            end
            if obj.packet_buffer(read_temp).id 
                [N_data_bits,part] = obj.packet_buffer(read_temp).send_data(N_data_bits);
                packet_parts = [packet_parts,part];
                if N_data_bits <= 0
                    break;
                end
            end
        end
        obj.bit_count = obj.get_buffer_length;
   end
         
end

end