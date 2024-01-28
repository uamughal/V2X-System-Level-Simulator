classdef miscUtils
    % Implements miscellaneous functions related to the trace generation.
    % Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at.
    % (c) 2011 by INTHFT
    % www.nt.tuwien.ac.at
    
    properties
    end
    
    methods(Static)
        function LTE_params = LTE_params_function
            % Re-create needed load_parameters data from Link level for the generation of the precoding matrices.
            % (c) Josep Colom Ikuno, INTHFT, 2008/2011
            % www.nt.tuwien.ac.at
            
            
            %% Create the Codebook for Precoding
            
            % Transmit diversity
            LTE_params.Z{1} =  [1, 0, 1i,  0;
                0,-1,  0, 1i;
                0, 1,  0, 1i;
                1, 0,-1i,  0];
            LTE_params.Z{2} =  [1, 0, 0, 0, 1i,  0,  0, 0;
                0, 0, 0, 0,  0,  0,  0, 0;
                0,-1, 0, 0,  0, 1i,  0, 0;
                0, 0, 0, 0,  0,  0,  0, 0;
                0, 1, 0, 0,  0, 1i,  0, 0;
                0, 0, 0, 0,  0,  0,  0, 0;
                1, 0, 0, 0,-1i,  0,  0, 0;
                0, 0, 0, 0,  0,  0,  0, 0;
                0, 0, 0, 0,  0,  0,  0, 0;
                0, 0, 1, 0,  0,  0, 1i, 0;
                0, 0, 0, 0,  0,  0,  0, 0;
                0, 0, 0,-1,  0,  0,  0,1i;
                0, 0, 0, 0,  0,  0,  0, 0;
                0, 0, 0, 1,  0,  0,  0,1i;
                0, 0, 0, 0,  0,  0,  0, 0;
                0, 0, 1, 0,  0,  0,-1i, 0];
            
            % Spatial multiplexing
            U_temp = [  1,-1,-1,-1;     % Matrix corresponding to vectors u0 ... u15 in Table 6.3.4.2.3-2
                1,-1i,1,1i;
                1,1,-1,1;
                1,1i,1,-1i;
                1,(-1-1i)/sqrt(2), -1i,(1-1i)/sqrt(2);
                1,(1-1i)/sqrt(2), 1i,(-1-1i)/sqrt(2);
                1,(1+1i)/sqrt(2), -1i,(-1+1i)/sqrt(2);
                1,(-1+1i)/sqrt(2), 1i,(1+1i)/sqrt(2);
                1,-1,1,1;
                1,-1i,-1,-1i;
                1,1,1,-1;
                1,1i,-1,1i;
                1,-1,-1,1;
                1,-1,1,-1;
                1,1,-1,-1;
                1,1,1,1;].';
            Wn = zeros(4,4,16);
            for ii = 1:16
                LTE_params.Wn(:,:,ii)=diag(ones(1,4))-2*U_temp(:,ii)*U_temp(:,ii)'/(U_temp(:,ii)'*U_temp(:,ii));
            end
            
            % W Matrix according to Table 6.3.4.2.3-1
            %  LTE_params.W{1} = cat(3,[1;0],[0;1],[1/sqrt(2);1/sqrt(2)],[1/sqrt(2);-1/sqrt(2)],...
            %         [1/sqrt(2);1i/sqrt(2)],[1/sqrt(2);-1i/sqrt(2)]);
            LTE_params.W{1} = cat(3,[1/sqrt(2);1/sqrt(2)],[1/sqrt(2);-1/sqrt(2)],...
                [1/sqrt(2);1i/sqrt(2)],[1/sqrt(2);-1i/sqrt(2)]);
            LTE_params.W{2} = cat(3,1/sqrt(2)*[1,0;0,1],1/(2)*[1,1;1,-1],1/(2)*[1,1;1i,-1i]);
            
            % Large delay CDD
            LTE_params.U_l{1} = 1;
            LTE_params.U_l{2} = 1/sqrt(2)*[1,1;1,exp(-1i*pi)];
            LTE_params.U_l{3} = 1/sqrt(3)*[1,1,1;1,exp(-1i*2*pi/3),exp(-1i*4*pi/3);1,exp(-1i*4*pi/3),exp(-1i*8*pi/3)];
            LTE_params.U_l{4} = 1/2*[1,1,1,1;1,exp(-1i*2*pi/4),exp(-1i*4*pi/4),exp(-1i*6*pi/4);...
                1,exp(-1i*4*pi/4),exp(-1i*8*pi/4),exp(-1i*12*pi/4);...
                1,exp(-1i*6*pi/4),exp(-1i*12*pi/4),exp(-1i*18*pi/4)];
            LTE_params.D_l{1} = 1;
            LTE_params.D_l{2} = [1,0;0,exp(-1i*pi)];
            LTE_params.D_l{3} = [1,0,0;0,exp(-1i*2*pi/3),0;0,0,exp(-1i*4*pi/3)];
            LTE_params.D_l{4} = [1,0,0,0;0,exp(-1i*2*pi/4),0,0;0,0,exp(-1i*4*pi/4),0;0,0,0,exp(-1i*6*pi/4)];
            
            % Note that as of v.8.3.0, small delay CDD is removed from the standard
            % (28/05/08	RAN_40	RP-080432	0043	-	Removal of small-delay CDD
            
            % Precoding matrix W columns to take for each layer mapping
            LTE_params.mapping{1} = ones(16,1);
            LTE_params.mapping{2}=[1 4;1 2;1 2;1 2;1 4;1 4;1 3;1 3;1 2;1 4;1 3;1 3;1 2;1 3;1 3;1 2];
            LTE_params.mapping{3}=[1 2 4;1 2 3;1 2 3;1 2 3;1 2 4;1 2 4;1 3 4;1 3 4;1 2 4;1 3 4;1 2 3;1 3 4;1 2 3;1 2 3;1 2 3;1 2 3];
            LTE_params.mapping{4}=[1 2 3 4;1 2 3 4;3 2 1 4;3 2 1 4;1 2 3 4;1 2 3 4;1 3 2 4;
                1 3 2 4;1 2 3 4;1 2 3 4;1 3 2 4;1 3 2 4;1 2 3 4;1 3 2 4;3 2 1 4;1 2 3 4];
        end
        
        function precoding_config = get_all_precoding_combinations
            % This small helper function returns all possible precoding options for LTE.
            % (c) Josep Colom Ikuno, INTHFT, 2008
            % www.nt.tuwien.ac.at
            
            precoding_config = cell(4,4,4); % Up to 4 layers, up to 4 TX antennas, 4 TX modes
            LTE_params       = phy_modeling.miscUtils.LTE_params_function;
            
            for tx_mode = 1:4
                switch tx_mode
                    case 1
                        % SISO
                        precoding_config{1,1,tx_mode}.tx_mode = tx_mode;
                        precoding_config{1,1,tx_mode}.nAtPort = 1;
                        precoding_config{1,1,tx_mode}.nLayers = 1;
                    case 2
                        % TxD
                        for nAtPort = [2 4]
                            precoding_config{nAtPort,nAtPort,tx_mode}.tx_mode = tx_mode;
                            precoding_config{nAtPort,nAtPort,tx_mode}.nAtPort = nAtPort;
                            precoding_config{nAtPort,nAtPort,tx_mode}.nLayers = nAtPort;
                            
                            %% Codebook setting
                            % We call the precoding matrix of TxD Z: Matrix corresponding to 36.211 section 6.3.4.3
                            precoding_config{nAtPort,nAtPort,tx_mode}.Z       = LTE_params.Z{nAtPort/2};
                            precoding_config{nAtPort,nAtPort,tx_mode}.name    = 'TxD';
                        end
                    case {3 4}
                        % OLSM, CLSM
                        for nAtPort = [2 4]
                            for nLayers = 1:nAtPort
                                precoding_config{nLayers,nAtPort,tx_mode}.tx_mode = tx_mode;
                                precoding_config{nLayers,nAtPort,tx_mode}.nAtPort = nAtPort;
                                precoding_config{nLayers,nAtPort,tx_mode}.nLayers = nLayers;
                                
                                %% Codebook setting
                                switch tx_mode
                                    case 4
                                        % CLSM
                                        switch nAtPort
                                            case 2
                                                switch nLayers
                                                    case 1
                                                        codebook_indexs = 0:3;
                                                    case 2
                                                        codebook_indexs = 1:2;
                                                end
                                            case 4
                                                codebook_indexs = 0:15;
                                        end
                                        
                                        % Closed loop spatial multiplexing, section 6.3.4.2.1
                                        W = zeros(nAtPort,nLayers,length(codebook_indexs));
                                        if (nAtPort == 2)
                                            if (min(codebook_indexs)<0 || max(codebook_indexs)>3) && nLayers ==1
                                                error('Only codebooks 0-3 are defined for %d layers (see TS.36.211, Table 6.3.4.2.3-1)',nLayers);
                                            elseif (min(codebook_indexs)<0 || max(codebook_indexs)>2) && nLayers ==2
                                                error('Only codebooks 0-2 are defined for %d layers (see TS.36.211, Table 6.3.4.2.3-1)',nLayers);
                                            end
                                            for cb_ = 1:length(codebook_indexs)
                                                codebook_index = codebook_indexs(cb_);
                                                W(:,:,cb_) = LTE_params.W{nLayers}(:,:,codebook_index+1);
                                            end
                                        else
                                            if min(codebook_indexs)<0 || max(codebook_indexs)>15
                                                error('Only codebooks 0-15 are defined (see TS.36.211, Table 6.3.4.2.3-2)');
                                            end
                                            for cb_ = 1:length(codebook_indexs)
                                                codebook_index = codebook_indexs(cb_);
                                                W_temp = 1/sqrt(nLayers)*LTE_params.Wn(:,:,codebook_index+1);
                                                W(:,:,cb_) = W_temp(:,LTE_params.mapping{nLayers}(codebook_index+1,:),1);
                                            end
                                        end
                                        precoding_config{nLayers,nAtPort,tx_mode}.W              = W;
                                        precoding_config{nLayers,nAtPort,tx_mode}.name           = 'CLSM';
                                        precoding_config{nLayers,nAtPort,tx_mode}.codebook_index = codebook_indexs;
                                        
                                    case 3
                                        % OLSM uses codebooks 12-15 in a cyclic way. Thus we set codebook_index to [12 13 14 15]
                                        switch nAtPort
                                            case 4
                                                codebook_index = [12 13 14 15]-1;
                                            case 2
                                                codebook_index = 0;
                                        end
                                        
                                        if (nAtPort == 2)
                                            W = LTE_params.W{nLayers}(:,:,codebook_index+1);
                                        else
                                            W_temp = 1/sqrt(nLayers)*LTE_params.Wn(:,:,codebook_index+1);
                                            % nLayers long Cyclic precoding matrix
                                            W = zeros(4,nLayers,4);
                                            for ii = 13:16
                                                W(:,:,ii-12) = W_temp(:,LTE_params.mapping{nLayers}(ii,:),ii-12);
                                            end
                                        end
                                        
                                        % Open loop spatial multiplexing, section 6.3.4.2.2 (Large CDD)
                                        precoding_config{nLayers,nAtPort,tx_mode}.W              = W;
                                        precoding_config{nLayers,nAtPort,tx_mode}.U              = LTE_params.U_l{nLayers};
                                        precoding_config{nLayers,nAtPort,tx_mode}.D              = LTE_params.D_l{nLayers};
                                        precoding_config{nLayers,nAtPort,tx_mode}.name           = 'OLSM';
                                        precoding_config{nLayers,nAtPort,tx_mode}.codebook_index = codebook_index;
                                end
                            end
                        end
                end
            end
        end
    end
    
end