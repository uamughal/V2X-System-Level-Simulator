classdef WinnerC2PathlossModel < macroscopic_pathloss_models.generalPathlossModel
    % Propagation conditions as proposed by TS 36.942 V8.0.0
    % (c) Josep Colom Ikuno, INTHFT, 2009
    % www.nt.tuwien.ac.at
    properties
        frequency    % Frequency in HERTZs (for consistency, although this 
                     % model makes calculations with frequencies set in MHz
        environment  % Environment that this instance represents

        Hbs          % The base station antenna height in metres
        Hms          % The mobile station antenna height in metres
        
        pathloss_function_handle

    end

    methods
        % Class constructor
        function obj = WinnerC2PathlossModel(frequency,environment)
            obj.frequency   = frequency;
            obj.environment = environment;
            obj.name        = 'Urban macro-cell';
            
            % Dhb is the base station antenna height in metres, measured
            % from the average rooftop level (TS 36.942)
            switch environment
                case 'LOS'
                    obj.name = [obj.name ' LOS '];
                    obj.pathloss_function_handle = 1;
                    % Using suggested values in TS 36.942
                    obj.Hbs = 25;
                    obj.Hms = 1.5;
                case 'NLOS'
                    obj.name = [obj.name ' NLOS '];
                    obj.pathloss_function_handle = 2;
                    % Using suggested values in TS 36.942
                    obj.Hbs = 25;
                    obj.Hms = 1.5;
                otherwise
                    error(['"' environment '"" environment not valid']);
            end
        end
        
        % Returns the macroscopic pathloss in dB. Note: distance in METERS
        function pathloss_in_db = pathloss(obj,distance)
            % Restrict that pathloss must be bigger than MCL
            switch obj.pathloss_function_handle
                case 1
                    pathloss_in_db = obj.pathloss_los(distance);
                case 2
                    pathloss_in_db = obj.pathloss_nlos(distance);
            end
        end
        
        % LOS pathloss
        function pl = pathloss_los(obj,distance)
            % Macro cell propagation model for urban area (LOS case) 
            %           distance ... actual distance in m
            % output:   pl       ... NLOS pathloss in dB            
            % Calculations are done in m       
            hbs=25-1;
            hms=1.5-1;
            PL_bp=4*hbs*hms*obj.frequency/2.998e8;
            pl=40*log10(distance) + 13.47 +6*log10(obj.frequency/5) - 14*log10(hms) - 14*log10(hbs);
            %{
            distance0=reshape(distance,1,size(distance,1)*size(distance,2));
            for i=1:length(distance0)
                
            if distance0(i)<=PL_bp & distance0(i)>=10
                pL(i)= 26*log10(distance(i)) + 39 + 20*log10(obj.frequency/5);
            end
            if distance0(i) >= PL_bp & distance0(i) <= 5000
               pL(i)= 40*log10(distance(i)) + 13.47 +6*log10(obj.frequency/5) - 14*log10(hms) - 14*log10(hbs);
            end
            end
               pl=reshape(pL,size(distance,1),size(distance,2));
         %}
        end
        
        % NLOS pathloss
        function pl = pathloss_nlos(obj,distance)
            % Macro cell propagation model for urban area (NLOS case)
            %           distance ... actual distance in m
            % output:   pl       ... NLOS pathloss in dB   
            hbs=obj.Hbs-1;
            pl=(44.9-6.55*log10(hbs)).*log10(distance)+34.46+5.83*log10(hbs)+23*log10(obj.frequency/5.0);                  
        end
    end
end