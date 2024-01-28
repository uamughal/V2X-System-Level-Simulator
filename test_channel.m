clc
clear
%% Generate Pathloss for link Vehicle to Vehicle

BsHeight = 1.5;
MsHeight = 1.5;
CenterFrequency = 6e9;     
H_bs = BsHeight-1; %effective environment height
H_ms = MsHeight-1; %effective environment height
PL_bp = 4*H_bs.*H_ms*CenterFrequency/2.998e8; 

j=1;
k=1;
for i=10:1:PL_bp+1           
               loss_LoS1(i-9) = 22.7*log10(i) + 41.0 + 20*log10(CenterFrequency/5);
               distance1(j)=i;
               j=j+1;

end
figure 
plot(distance1,loss_LoS1);hold on;
grid on;

for i=22:1:5000
             loss_LoS2(i-21) = 40*log10(i) + 9.45 - 17.3*log10(H_bs) - 17.3*log10(H_ms) + 2.7*log10(CenterFrequency/5);  
             distance2(k)=i;
             k=k+1;
end

figure 
plot(distance2,loss_LoS2);hold on;
grid on;


xlabel('Distance (m)');
ylabel('Pathloss');
title('Urban micro-cell pathloss');



