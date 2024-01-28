% Arrays=example_syntetic_arrays
% returns some example Arrays
%
% Author: Martin K�ske (TUI)

function Arrays=example_syntetic_arrays

NAz=120; %3 degree sampling interval
Az=linspace(-180,180-1/NAz,NAz);

pattern(1,:,1,:)=dipole(Az,12); % slanted by 12 degree

Arrays(1)=AntennaArray('ULA',2,0.01,'FP-ECS',pattern,'Azimuth',Az); %ULA-2 1cm spacing
Arrays(2)=AntennaArray('ULA',4,0.01,'FP-ECS',pattern,'Azimuth',Az); %ULA-4 1cm spacing
Arrays(3)=AntennaArray('ULA',8,0.01,'FP-ECS',pattern,'Azimuth',Az); %ULA-8 1cm spacing
Arrays(4)=AntennaArray('UCA',4,0.01,'FP-ECS',pattern,'Azimuth',Az); %UCA-4 1cm radius
Arrays(5)=AntennaArray('UCA',8,0.01,'FP-ECS',pattern,'Azimuth',Az); %UCA-8 1cm radius

NAz=3*120; %3 degree sampling interval
Az=linspace(-180,180-1/NAz,NAz);
pattern=ones(2,2,1,NAz);
dist = 3e8/5.25e9*0.5;
Arrays(6)=AntennaArray('ULA',2,dist,'FP-ECS',pattern); % isotropic antenna
