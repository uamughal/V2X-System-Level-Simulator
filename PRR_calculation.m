function PRR = PRR_calculation(UEs) %  Av_PRR0 Av_PRR
%%

for u_=1:length(UEs)
 if ~isempty(UEs(u_).connected_UEs(1).id)
  x=0;
  y=length(UEs(u_).connected_UEs);
  for i=1:length(UEs(u_).connected_UEs)
       UEs_TB_size = UEs(u_).connected_UEs(i).trace.TB_size;
       ACK = UEs(u_).connected_UEs(i).trace.ACK;
       ind1 = find((UEs_TB_size(1,1:10)>0));
       ind2 = find (ACK(1,1:10)>0);
       if isempty(ind1) && isempty(ind2)
           y=y-1;
       end
       if length (ind1)==length(ind2)
         if ind1 ==ind2 
          x=x+1;
         end
       end
  end
  PRR(u_)=x./y;
  Packet_RR=fopen('Packet_RR.txt','at');
  fprintf(Packet_RR,'%f\n',PRR);
  fclose(Packet_RR);
  end
end

%%
%{
    for u_=1:length(UEs)
        if ~isempty(UEs(u_).connected_UEs(1).id)
        UEs_TB_size = UEs(u_).connected_UEs(1).trace.TB_size;
        expected_ACKs = UEs(u_).sectors(1).sector_trace.expected_ACKs;
        UE_acknowledged_datas = UEs(u_).sectors(1).sector_trace.acknowledged_data;
        ind = find((UEs_TB_size(1,1:10)>0));
        ind=ind(1:length(ind)-1);
        
        UE_acknowledged_datas0=double(UE_acknowledged_datas(1,ind));
        UEs_TB_size0=double(UEs_TB_size(1,ind));
        expected_ACKs0=double(expected_ACKs(1,ind));
        if ~isempty(ind)
            for u=1:length(ind)
                p_list(u_,u)=UE_acknowledged_datas0(u)./UEs_TB_size0(u)./ expected_ACKs0(u);            
            end
        end
       end
    end                
    PRR=p_list;
    [Inf_id1, Inf_id2]= find(PRR(:,:)==Inf);
    PRR(:,Inf_id2)=0;
    for i=1:size(PRR,1)
    non_zeros_id0 = find(PRR(i,:)>0);
    Av_PRR0(i,1) = sum(PRR(i,:),2)./length(PRR(i,:));    
    end
    
    Av_PRR0 = Av_PRR0(~isnan(Av_PRR0));
    non_zeros_id = find(Av_PRR0(:,1)>0);
    Av_PRR0 = Av_PRR0(non_zeros_id);
    non_zeros_id1 = find(Av_PRR0(:,1)<=2);
    Av_PRR = sum(Av_PRR0(non_zeros_id1))./length(non_zeros_id1);
  
    Packet_RR=fopen('Packet_RR.txt','at');
    fprintf(Packet_RR,'%f\n',Av_PRR);
    fclose(Packet_RR);
%}
end   
   
