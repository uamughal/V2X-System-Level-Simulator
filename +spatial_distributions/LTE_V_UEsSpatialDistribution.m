classdef LTE_V_UEsSpatialDistribution < spatial_distributions.networkElementSpatialDistribution
    % Homogeneous spatial distribution of network elements (HeNBs, Small Cells, Microcells, etc.). Assumes an already existing macrocell network and acts as overlay.
    % The ROI is the size of the networkPathlossMap in this case.
    % (c) Martin Taranetz, Josep Colom Ikuno, ITC, 2012
    
    properties
         nUEs   %[a] a LTE-V UEs /sector
        % values used to determine LTE-V UEs positions
    
         PS_isd   % ISD between among eNBs
         Vehicle_isd  % ISD between among Vehicles;
         lane_width % 
    end
    
    methods
         % Class constructor.
       function obj = LTE_V_UEsSpatialDistribution(networkPathlossMap, N_UEs, ps_ISD, vehicles_isd,tr_line_width)
         obj                    = obj@spatial_distributions.networkElementSpatialDistribution(networkPathlossMap);
         obj.nUEs   = N_UEs;                    
        
         obj.PS_isd=ps_ISD ;
         obj.Vehicle_isd=vehicles_isd;
         obj.lane_width=tr_line_width;                      
       end
       
       % Generate position [m]
       function elements_positions = generate_positions(obj)    
 %% LTE-V UEs (changed)
           n_V_UE=obj.nUEs;
           ps_ISD=obj.PS_isd;
           Lane_width=obj.lane_width;
           elements_positions=zeros(n_V_UE,3);
           reference_pos=[0,0];
           roi_min=obj.networkPathlossMap.coordinate_origin;
           data_res=obj.networkPathlossMap.data_res; 
           
%%   Spatial Poisson Distribution  
   %% For all vertical lanes (4)
    Lambda1 = 31;  % Lambda:Vehicle density
    u = unifrnd(0,1);
    M = 0;
   while u >= exp(-Lambda1)
        u = u*unifrnd(0,1);
        M=M+1;
   end 
   a = 0; b = 1299;
   Nall = M;
   while M > 0         %scatter in the [0,1299]
        M = M-1;
        u1 = unifrnd(0,1);
        A11(Nall-M) = ps_ISD/6-Lane_width*3/2;
        B11(Nall-M) = -433-(ps_ISD/sqrt(3)/2)+(b-a)*u1;
        A12(Nall-M) = ps_ISD/6-Lane_width/2;
        B12(Nall-M) = -433-(ps_ISD/sqrt(3)/2)+(b-a)*u1;
        A13(Nall-M) = ps_ISD/6+Lane_width/2;
        B13(Nall-M) = -433-(ps_ISD/sqrt(3)/2)+(b-a)*u1;
        A14(Nall-M) = ps_ISD/6+Lane_width*3/2;
        B14(Nall-M) = -433-(ps_ISD/sqrt(3)/2)+(b-a)*u1;
        
        A21(Nall-M) = 250+ps_ISD/6-Lane_width*3/2; 
        B21(Nall-M) = -433-(ps_ISD/sqrt(3)/2)+(b-a)*u1;
        A22(Nall-M) = 250+ps_ISD/6-Lane_width/2;
        B22(Nall-M) = -433-(ps_ISD/sqrt(3)/2)+(b-a)*u1;
        A23(Nall-M) = 250+ps_ISD/6+Lane_width/2;
        B23(Nall-M) = -433-(ps_ISD/sqrt(3)/2)+(b-a)*u1;
        A24(Nall-M) = 250+ps_ISD/6+Lane_width*3/2;
        B24(Nall-M) = -433-(ps_ISD/sqrt(3)/2)+(b-a)*u1;
        
        A31(Nall-M) = -250+ps_ISD/6-Lane_width*3/2;
        B31(Nall-M) = -433-(ps_ISD/sqrt(3)/2)+(b-a)*u1;
        A32(Nall-M) = -250+ps_ISD/6-Lane_width/2;
        B32(Nall-M) = -433-(ps_ISD/sqrt(3)/2)+(b-a)*u1;
        A33(Nall-M) = -250+ps_ISD/6+Lane_width/2;
        B33(Nall-M) = -433-(ps_ISD/sqrt(3)/2)+(b-a)*u1;
        A34(Nall-M) = -250+ps_ISD/6+Lane_width*3/2;
        B34(Nall-M) = -433-(ps_ISD/sqrt(3)/2)+(b-a)*u1;
        
        A41(Nall-M) = -2*250+ps_ISD/6-Lane_width*3/2;
        B41(Nall-M) = -433-(ps_ISD/sqrt(3)/2)+(b-a)*u1;
        A42(Nall-M) = -2*250+ps_ISD/6-Lane_width/2;
        B42(Nall-M) = -433-(ps_ISD/sqrt(3)/2)+(b-a)*u1;
        A43(Nall-M) = -2*250+ps_ISD/6+Lane_width/2;
        B43(Nall-M) = -433-(ps_ISD/sqrt(3)/2)+(b-a)*u1;
        A44(Nall-M) = -2*250+ps_ISD/6+Lane_width*3/2;
        B44(Nall-M) = -433-(ps_ISD/sqrt(3)/2)+(b-a)*u1;        
   end
   
   
 %% For all horizontal lanes(4)
    Lambda2 = 18;  % Lambda:Vehicle density
    u0 = unifrnd(0,1);
    M1 = 0;
   while u0 >= exp(-Lambda2)
        u0 = u0*unifrnd(0,1);
        M1=M1+1;
   end 
   c = 0; d = 750;
   Nall0 = M1;
   while M1 > 0         %scatter in the [0,433]
        M1 = M1-1;
        u2 = unifrnd(0,1);
        C11(Nall0-M1) = -250-ps_ISD/3+(d-c)*u2;
        D11(Nall0-M1) = -433-(ps_ISD/sqrt(3)/2-3*Lane_width/2);
        C12(Nall0-M1) = -250-ps_ISD/3+(d-c)*u2;
        D12(Nall0-M1) = -433-(ps_ISD/sqrt(3)/2-Lane_width/2);
        C13(Nall0-M1) = -250-ps_ISD/3+(d-c)*u2;
        D13(Nall0-M1) = -433-(ps_ISD/sqrt(3)/2+Lane_width/2);
        C14(Nall0-M1) = -250-ps_ISD/3+(d-c)*u2;
        D14(Nall0-M1) = -433-(ps_ISD/sqrt(3)/2+3*Lane_width/2);
        
        C21(Nall0-M1) = -250-ps_ISD/3+(d-c)*u2;
        D21(Nall0-M1) = -(ps_ISD/sqrt(3)/2-3*Lane_width/2);
        C22(Nall0-M1) = -250-ps_ISD/3+(d-c)*u2;
        D22(Nall0-M1) = -(ps_ISD/sqrt(3)/2-Lane_width/2);
        C23(Nall0-M1) = -250-ps_ISD/3+(d-c)*u2;
        D23(Nall0-M1) = -(ps_ISD/sqrt(3)/2+Lane_width/2);
        C24(Nall0-M1) = -250-ps_ISD/3+(d-c)*u2;
        D24(Nall0-M1) = -(ps_ISD/sqrt(3)/2+3*Lane_width/2);
        
        C31(Nall0-M1) = -250-ps_ISD/3+(d-c)*u2;
        D31(Nall0-M1) = 433-(ps_ISD/sqrt(3)/2-3*Lane_width/2);
        C32(Nall0-M1) = -250-ps_ISD/3+(d-c)*u2;
        D32(Nall0-M1) = 433-(ps_ISD/sqrt(3)/2-Lane_width/2);
        C33(Nall0-M1) = -250-ps_ISD/3+(d-c)*u2;
        D33(Nall0-M1) = 433-(ps_ISD/sqrt(3)/2+Lane_width/2);
        C34(Nall0-M1) = -250-ps_ISD/3+(d-c)*u2;
        D34(Nall0-M1) = 433-(ps_ISD/sqrt(3)/2+3*Lane_width/2);
        
        C41(Nall0-M1) = -250-ps_ISD/3+(d-c)*u2;
        D41(Nall0-M1) = 2*433-(ps_ISD/sqrt(3)/2-3*Lane_width/2);
        C42(Nall0-M1) = -250-ps_ISD/3+(d-c)*u2;
        D42(Nall0-M1) = 2*433-(ps_ISD/sqrt(3)/2-Lane_width/2);
        C43(Nall0-M1) = -250-ps_ISD/3+(d-c)*u2;
        D43(Nall0-M1) = 2*433-(ps_ISD/sqrt(3)/2+Lane_width/2);
        C44(Nall0-M1) = -250-ps_ISD/3+(d-c)*u2;
        D44(Nall0-M1) = 2*433-(ps_ISD/sqrt(3)/2+3*Lane_width/2);        
   end             
      A=[A11,A12,A13,A14,A21,A22,A23,A24,A31,A32,A33,A34,A41,A42,A43,A44,...
          C11,C12,C13,C14,C21,C22,C23,C24,C31,C32,C33,C34,C41,C42,C43,C44];
      B=[B11,B12,B13,B14,B21,B22,B23,B24,B31,B32,B33,B34,B41,B42,B43,B44,...
          D11,D12,D13,D14,D21,D22,D23,D24,D31,D32,D33,D34,D41,D42,D43,D44];
      for i=1:length(A)
          elements_positions(i,:)=[A(i),B(i),0];
      end
      
       %% Set moving direction for each vehicle (vertical lanes)
       for j=1:length(A11)+length(A12)
           elements_positions(j,3)= -90;
       end
       for j=j+1:j+1+length(A13)+length(A14)
           elements_positions(j,3)= 90;
       end
       for j=j+1:j+1+length(A21)+length(A22)
           elements_positions(j,3)= -90;
       end
       for j=j+1:j+1+length(A23)+length(A24)
           elements_positions(j,3)= 90;
       end
       for j=j+1:j+1+length(A31)+length(A32)
           elements_positions(j,3)= -90;
       end
       for j=j+1:j+1+length(A33)+length(A34)
           elements_positions(j,3)= 90;
       end
       for j=j+1:j+1+length(A41)+length(A42)
           elements_positions(j,3)= -90;
       end
       for j=j+1:j+1+length(A43)+length(A44)
           elements_positions(j,3)= 90;
       end
       
        %% Set moving direction for each vehicle (horizontal lanes)
       for j=j+1:j+1+length(C11)+length(C12)
           elements_positions(j,3)= 180;
       end
       for j=j+1:j+1+length(C13)+length(C14)
           elements_positions(j,3)= 0;
       end
       for j=j+1:j+1+length(C21)+length(C22)
           elements_positions(j,3)= 180;
       end
       for j=j+1:j+1+length(C23)+length(C24)
           elements_positions(j,3)= 0;
       end
       for j=j+1:j+1+length(C31)+length(C32)
           elements_positions(j,3)= 180;
       end
       for j=j+1:j+1+length(C33)+length(C34)
           elements_positions(j,3)= 0;
       end
       for j=j+1:j+1+length(C41)+length(C42)
           elements_positions(j,3)= 180;
       end
       for j=j+1:j+1+length(C43)+length(C44)
           elements_positions(j,3)= 0;
       end
       
           %% Give UE position on East side (4 lanes, 80 UEs)
           %{ 
           for i=1:5     
              r_ue_x=reference_pos(1)+ps_ISD/6-Lane_width*3/2;
              r_ue_y=reference_pos(2)+(ps_ISD/sqrt(3)/2)-(i-1)*50;
              elements_positions(i,:)=[r_ue_x,r_ue_y];
           end
           
           for i=6:10     
              r_ue_x=reference_pos(1)+ps_ISD/6-Lane_width/2;
              r_ue_y=reference_pos(2)+(ps_ISD/sqrt(3)/2)-(i-6)*50;
              elements_positions(i,:)=[r_ue_x,r_ue_y];
           end
           
           for i=11:15     
              r_ue_x=reference_pos(1)+ps_ISD/6+Lane_width/2;
              r_ue_y=reference_pos(2)-(ps_ISD/sqrt(3)/2)+(i-11)*50;
              elements_positions(i,:)=[r_ue_x,r_ue_y];
           end
           
           for i=16:20     
              r_ue_x=reference_pos(1)+ps_ISD/6+3*Lane_width/2;
              r_ue_y=reference_pos(2)-(ps_ISD/sqrt(3)/2)+(i-16)*50;
              elements_positions(i,:)=[r_ue_x,r_ue_y];
           end
            
           %% Give UE position on south side (4 lanes, 24 UEs)
           for i=21:26
              r_ue_x=reference_pos(1)+(ps_ISD/6-Lane_width*2)-(i-21)*50;
              r_ue_y=reference_pos(2)-(ps_ISD/sqrt(3)/2-3*Lane_width/2);
              elements_positions(i,:)=[r_ue_x,r_ue_y];
           end   
           for i=27:32
              r_ue_x=reference_pos(1)+(ps_ISD/6-Lane_width*2)-(i-27)*50;
              r_ue_y=reference_pos(2)-(ps_ISD/sqrt(3)/2-Lane_width/2);
              elements_positions(i,:)=[r_ue_x,r_ue_y];
           end
           for i=33:38
              r_ue_x=reference_pos(1)-(ps_ISD/2-Lane_width*2)+(i-33)*50;
              r_ue_y=reference_pos(2)-(ps_ISD/sqrt(3)/2+Lane_width/2);
              elements_positions(i,:)=[r_ue_x,r_ue_y];
           end
           for i=39:44
              r_ue_x=reference_pos(1)-(ps_ISD/2-Lane_width*2)+(i-39)*50;
              r_ue_y=reference_pos(2)-(ps_ISD/sqrt(3)/2+3*Lane_width/2);
              elements_positions(i,:)=[r_ue_x,r_ue_y];
           end
                      
           %% Give UE position on west side (4 lanes,2lanes 16UEs + 2UEs/2 lanes 10 + 4UEs)
           for i=45:52
              r_ue_x=reference_pos(1)-ps_ISD/3+3*Lane_width/2;
              r_ue_y=reference_pos(2)-(ps_ISD/sqrt(3)/2-2*Lane_width)+(i-45)*50;
              elements_positions(i,:)=[r_ue_x,r_ue_y];
           end
           for i=53
              r_ue_x=reference_pos(1)-ps_ISD/3+3*Lane_width/2;
              r_ue_y=reference_pos(2)-ps_ISD/sqrt(3)/2-2*Lane_width-50;
              elements_positions(i,:)=[r_ue_x,r_ue_y];
           end 
           
           %**************************************************************%
           for i=54:61
              r_ue_x=reference_pos(1)-ps_ISD/3+Lane_width/2;
              r_ue_y=reference_pos(2)-(ps_ISD/sqrt(3)/2-2*Lane_width)+(i-54)*50;
              elements_positions(i,:)=[r_ue_x,r_ue_y];
           end
           for i=62
              r_ue_x=reference_pos(1)-ps_ISD/3+Lane_width/2;
              r_ue_y=reference_pos(2)-ps_ISD/sqrt(3)/2-2*Lane_width-50;
              elements_positions(i,:)=[r_ue_x,r_ue_y];
           end
           
           %**************************************************************%
           for i=63:67
              r_ue_x=reference_pos(1)-ps_ISD/3-Lane_width/2;
              r_ue_y=reference_pos(2)+(ps_ISD/sqrt(3)-5*Lane_width)-(i-63)*50;
              elements_positions(i,:)=[r_ue_x,r_ue_y];
           end
           for i=68
              r_ue_x=reference_pos(1)-ps_ISD/3-Lane_width/2;
              r_ue_y=reference_pos(2)-4*Lane_width-50;
              elements_positions(i,:)=[r_ue_x,r_ue_y];
           end
           for i=69
              r_ue_x=reference_pos(1)-ps_ISD/3-Lane_width/2;
              r_ue_y=reference_pos(2)-ps_ISD/sqrt(3)/2-2*Lane_width-50;
              elements_positions(i,:)=[r_ue_x,r_ue_y];
           end
           %**************************************************************%
           for i=70:74
              r_ue_x=reference_pos(1)-ps_ISD/3-3*Lane_width/2;
              r_ue_y=reference_pos(2)+(ps_ISD/sqrt(3)-5*Lane_width)-(i-70)*50;
              elements_positions(i,:)=[r_ue_x,r_ue_y];
           end
           for i=75
              r_ue_x=reference_pos(1)-ps_ISD/3-3*Lane_width/2;
              r_ue_y=reference_pos(2)-4*Lane_width-50;
              elements_positions(i,:)=[r_ue_x,r_ue_y];
           end 
           for i=76
              r_ue_x=reference_pos(1)-ps_ISD/3-3*Lane_width/2;
              r_ue_y=reference_pos(2)-ps_ISD/sqrt(3)/2-2*Lane_width-50;
              elements_positions(i,:)=[r_ue_x,r_ue_y];
           end  
  
           %% Give UE position on north (2 lanes, 6 UEs)
           for i=77:79
              r_ue_x=reference_pos(1)-ps_ISD/3+(i-77)*50;
              r_ue_y=reference_pos(2)+ps_ISD/sqrt(3)-Lane_width/2;
              elements_positions(i,:)=[r_ue_x,r_ue_y];
           end
           for i=80:82
              r_ue_x=reference_pos(1)-ps_ISD/3+(i-80)*50;
              r_ue_y=reference_pos(2)+ps_ISD/sqrt(3)-3*Lane_width/2;
              elements_positions(i,:)=[r_ue_x,r_ue_y];
           end
        %}                                     
       end    
   end
end