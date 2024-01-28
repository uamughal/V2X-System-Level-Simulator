clc
clear

%%
   
    open('SINR to Spectral Efficiency.fig');
    h_line=get(gca,'Children');%get linehandles
    xdata=get(h_line,'Xdata');
    ydata=get(h_line,'Ydata');
    x1=xdata{1};
    y1=ydata{1};
    x2=xdata{6};
    y2=ydata{6};
    x3=xdata{7};
    y3=ydata{7};
    figure;
    plot(x1,y1);
    hold on;
    plot(x2,y2);
    hold on;
    plot(x3,y3); 
    hold on;
    
    x4=xdata{2};
    y4=ydata{2};
    scatter(x4,y4);
    hold on;
    
    x5=xdata{3};
    y5=ydata{3};
    scatter(x5,y5);
    hold on;
    
    x6=xdata{5};
    y6=ydata{5};
    scatter(x6,y6);
    hold on;
    
  


%     open('Outage Probability (2).fig');
%     h_line=get(gca,'Children');%get linehandles
%     xdata=get(h_line,'Xdata');
%     ydata=get(h_line,'Ydata');
%     x1=xdata{1};
%     y1=ydata{1};
%     x2=xdata{2};
%     y2=ydata{2};
%     x3=xdata{3};
%     y3=ydata{3};
%     
%     figure;
%     plot(x1,y1);
%     hold on;
%     plot(x2,y2);
%     hold on;
%     plot(x3,y3); 
%     hold on;
%     

        
    
%%

% outage=load('Outage probability.txt');
% outage1=load('Outage probability1.txt');
%     figure;
%     plot(10:20:250,outage);
%     hold on;
%     plot(10:20:250,outage1);
%     hold on;



%%
%     Rx_SNR=load('Rx SNR.txt'); 
%     Rx_SNR1=load('Rx SNR1.txt'); 
%     Rx_SNR2=load('Rx SNR2.txt');
%     Rx_SNR3=load('Rx SNR3.txt');
%     
%     Rx_SNR1=Rx_SNR(find(Rx_SNR>-8.378)); % for 150m
%     Rx_SNR2=Rx_SNR(find(Rx_SNR> 1.356)); %for 50m 
% 
%     
%     Rx_SNR4=Rx_SNR(find(Rx_SNR> -7.287));
%     Rx_SNR5= Rx_SNR(find(Rx_SNR> 2.735));
%     
%     
%     figure;
%     cdfplot(Rx_SNR);
%     hold on;
%     cdfplot(Rx_SNR3);
%     hold on;
%     cdfplot(Rx_SNR1);
%     hold on;
%     cdfplot(Rx_SNR4);
%     hold on;
%     cdfplot(Rx_SNR2);
%     hold on;
%     cdfplot(Rx_SNR5);
%     hold on;
    


%%
%     prr=load('Packet_RR.txt');
%     prr1=load('Packet_RR1.txt');
%     figure;
%     plot(10:20:250,prr);
%     hold on;
%     plot(10:20:250,prr1);
%     hold on;


% %%
% interference1= load('interference_without_CoMP.txt');
% interference1=interference1(find(interference1>-80));
% interference1=unique(interference1);
% interference2= load('interference_with_CoMP.txt');
% interference2=unique(interference2);
% interference3= load('interference_with_both.txt');
% interference3=unique(interference3);
% % interference4= load('interference_r_center.txt');
% 
% figure;
%   cdfplot(interference1);
%   hold on;  
%   cdfplot(interference2);
%   hold on;
%   cdfplot(interference3);
%   hold on;
  
%%
% snr1= load('snr_with_CoMP.txt');
% snr1=snr1(find(snr1>-20));
% snr1=unique(snr1);
% snr2= load('snr_without_CoMP.txt');
% snr2=unique(snr2);
% snr3= load('snr_with_both.txt');
% snr3=unique(snr3);
% 
% figure;
%   cdfplot(snr1);
%   hold on;  
%   cdfplot(snr2);
%   hold on;
%   cdfplot(snr3);
%   hold on;





%%
% UE_throughput=load('UE Throughput.txt');
% UE_throughput=unique(UE_throughput);
% UE_throughout=UE_throughput(find(UE_throughput<0.6));

%%

% UE_throughput1= load('UE Throughput1.txt');
% UE_throughput1=unique(UE_throughput1);
% UE_throughout1=UE_throughput1(find(UE_throughput1<0.6));

%%

% UE_throughput2= load('UE Throughput2.txt');
% UE_throughput2=unique(UE_throughput2);
% UE_throughout2=UE_throughput2(find(UE_throughput2<0.6));

%%
% figure;
%   cdfplot(UE_throughput2);
%   hold on;  
%   cdfplot(UE_throughput1);
%   hold on;
%   cdfplot(UE_throughput);
%   hold on;
%   title('UE Average Throughput')
%   hold on;
  
%%
  
% snr=load('Rx SNR.txt');
% 
%     ind1=find(snr(:)< 13 );
%     snr = snr(ind1);
%     ind2=find(snr(:)> -15);
%     snr = snr(ind2);
%     
% snr=snr(randperm(numel(snr)));
% 
%      fid = fopen('snr CoMP.txt','at');
%      fprintf(fid,'%0.6f \n',snr);
%      fclose(fid);
% 
% plot(1:length(snr),snr);

%%
  
 snr=load('SINR_with CoMP1.txt');
 figure;
 plot(6:length(snr)+5,snr);

 snr1=load('SINR_without CoMP.txt');
 figure;
 plot(6:length(snr1)+5,snr1);


     