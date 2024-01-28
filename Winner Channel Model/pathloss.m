function [loss, linkpar, fixpar, iterpar] = pathloss(wimpar,linkpar,fixpar,iterpar)
%PATHLOSS WIM pathloss models
%   PATH_LOSSES=PATHLOSS(WIMPAR,LINKPAR,FIXPAR,ITERPAR) returns path losses in dB scale
%   for all links defined in WIM input struct LINKPAR for the center
%   frequency and scenario given in WIMPAR. The output is a column vector
%   whose length is equal to the number of links defined in LINKPAR, e.g.
%   LENGTH(LINKPAR.MsBsDistance). The center frequencies and distances
%   must be specified in Herzes and meters, respectively.
%
%   Refs.   [1]: 3GPP TR 25.996 v6.1.0 (2003-09)
%           [2]: D1.1.2 V1.2, "WINNER II channel models"
%
%   See also WIMPARSET, LAYOUT2LINK, LAYOUTPARSET and SCENPARTABLES.

%   Authors: Lassi Hentil� (EBIT), Daniela Laselva (EBIT), Jari Salo (HUT),
%   Pekka Ky�sti (EBIT), Marko Milojevic (TUI), Mikko Alatossava (CWC/UOULU)

%   Revision history after WIM release:
%   Frequency and antenna height dependence added.    4.5.2006 HentLas
%   B5 scenarios added.                               4.5.2006 HentLas
%   B1 LOS&NLOS updated.                              8.9.2006 PekKy
%   D1.1.1 parameters, new scenarios added            12.1.2007 Marko
%   LoS/NloS condition renewed, D1.1.1 parameters     8.2.2007 MikkoA
%   Updated according to the D1.1.2 [2]               8.10.2007 HentLas



% extract required parameters from the input structs
NumLinks = length(iterpar.UserIndeces);
MsBsDistance=linkpar.MsBsDistance(iterpar.UserIndeces);
Scenario=iterpar.Scenario;
NumPenetratedFloors = linkpar.NumPenetratedFloors(iterpar.UserIndeces);
NumFloors = linkpar.NumFloors(iterpar.UserIndeces);
CenterFrequency=wimpar.CenterFrequency*1e-9; % Center frequency -> GHz

PropagCondition = iterpar.PropagCondition;
LoSConnectionLinks = find(PropagCondition);
NumLoSConnectionLinks = length(LoSConnectionLinks);
NLoSConnectionLinks = find(PropagCondition==0);
NumNLoSConnectionLinks = length(NLoSConnectionLinks);

SF_sigma = [];
if ~isempty(LoSConnectionLinks)
    if MsBsDistance(LoSConnectionLinks) > iterpar.LoS.PL_range(2)
        error('MsBsDistance exceeds the maximum allowed cell radius, see [2] table 4-4!')
    elseif MsBsDistance(LoSConnectionLinks) < iterpar.LoS.PL_range(1)
        error('MsBsDistance is below the minimum allowed cell radius, see [2] table 4-4!')
    end
end

if ~isempty(NLoSConnectionLinks)
    if MsBsDistance(NLoSConnectionLinks) > iterpar.NLoS.PL_range(2)
        error('MsBsDistance exceeds the maximum allowed cell radius, see [2] table 4-4!')
    elseif MsBsDistance(NLoSConnectionLinks) < iterpar.NLoS.PL_range(1)
        error('MsBsDistance is below the minimum allowed cell radius, see [2] table 4-4!')
    end
end


if isnan(linkpar.BsHeight(iterpar.UserIndeces(1)))
    warning('A1 Bs/MsHeights not defined by the user --> defauls taken from [2] table 4-4')
end

switch Scenario

    case {'A1'}
        %Check antenna heights
        if isnan(linkpar.BsHeight(iterpar.UserIndeces(1)))
            BsHeight = 1.5*ones(1,NumLinks);
            MsHeight = 1.5*ones(1,NumLinks);
        else
            BsHeight = linkpar.BsHeight(iterpar.UserIndeces);
            MsHeight = linkpar.MsHeight(iterpar.UserIndeces);
        end
        if NumPenetratedFloors > 0
            FL = 17 + 4*(NumPenetratedFloors - 1); % Floor loss
        else FL = 0;
        end

        if NumLoSConnectionLinks
            loss(LoSConnectionLinks) = iterpar.LoS.PL_B + iterpar.LoS.PL_A*log10(MsBsDistance(LoSConnectionLinks)) + iterpar.LoS.PL_C*log10(CenterFrequency/5) + FL;
            SF_sigma(LoSConnectionLinks) = 3;
        end

        if NumNLoSConnectionLinks
            if  strcmpi(wimpar.PathLossOption,'CR_light') % Corridor-to-Room light walls
                loss(NLoSConnectionLinks) = iterpar.NLoS.PL_A(1)*log10(MsBsDistance(NLoSConnectionLinks)) + iterpar.NLoS.PL_B(1)+ iterpar.LoS.PL_C*log10(CenterFrequency/5) + iterpar.NLoS.PL_X(1)+ FL;
                SF_sigma(NLoSConnectionLinks) = 4;
            end
            if  strcmpi(wimpar.PathLossOption,'CR_heavy') % Corridor-to-Room heavy walls
                loss(NLoSConnectionLinks) = iterpar.NLoS.PL_A(1)*log10(MsBsDistance(NLoSConnectionLinks)) + iterpar.NLoS.PL_B(1)+ iterpar.LoS.PL_C*log10(CenterFrequency/5) + iterpar.NLoS.PL_X(2)+ FL;
                SF_sigma(NLoSConnectionLinks) = 4;
            end
            if  strcmpi(wimpar.PathLossOption,'RR_light') % Room-to-Room through light walls
                loss(NLoSConnectionLinks) = 20*log10(MsBsDistance(NLoSConnectionLinks)) + iterpar.NLoS.PL_B(2) + floor(MsBsDistance(NLoSConnectionLinks)./10).*linkpar.WallLoss(1) + iterpar.LoS.PL_C*log10(CenterFrequency/5)+ iterpar.NLoS.PL_X(1)+ FL;
                SF_sigma(NLoSConnectionLinks) = 6;
            end
            if  strcmpi(wimpar.PathLossOption,'RR_heavy') % Room-to-Room through heavy walls
                loss(NLoSConnectionLinks) = 20*log10(MsBsDistance(NLoSConnectionLinks)) + iterpar.NLoS.PL_B(3) + floor(MsBsDistance(NLoSConnectionLinks)./10).*linkpar.WallLoss(2) + iterpar.LoS.PL_C*log10(CenterFrequency/5)+ iterpar.NLoS.PL_X(2)+ FL;
                SF_sigma(NLoSConnectionLinks) = 8;
            end
        end


    case {'A2', 'B4'}
        if strcmp(upper(Scenario),'A2')
            BsHeight = (2+3*(NumFloors-1)).*ones(1,NumLinks);
            MsHeight = 1.5*ones(1,NumLinks);
        else %B4
            BsHeight = 10*ones(1,NumLinks);
            MsHeight = (1.5+3*(NumFloors-1)).*ones(1,NumLinks);
        end

        StreetWidth = linkpar.StreetWidth(iterpar.UserIndeces);
        if isnan(linkpar.Dist1) % NaN default -> will be drawn randomly
            Dist_out = 1;
            while any(Dist_out > MsBsDistance | Dist_out < StreetWidth/2)
                Dist_out = MsBsDistance - (MsBsDistance-StreetWidth/2).*rand(1,length(MsBsDistance));
            end

        else %Dist1 is defined by the user
            Dist_out = linkpar.Dist1(iterpar.UserIndeces);
            Dist_out = Dist_out(NLoSConnectionLinks);
            for linkNum = 1:length(MsBsDistance) % check applicability and change Dist1 if needed
                if Dist_out(linkNum) > MsBsDistance(linkNum) || Dist_out(linkNum) < StreetWidth(linkNum)/2
                    Dist_out(linkNum) = MsBsDistance(linkNum) - (MsBsDistance(linkNum)-StreetWidth(linkNum)/2).*rand(1);
                end
            end
        end
        Dist_in = MsBsDistance - Dist_out; %indoor distance
        Theta = acos(StreetWidth./2./Dist_out)*180/pi; %angle from the BS to the normal of the wall

        %Free space loss
        FSL = 32.4+20*log10((Dist_out+Dist_in)/1000) + 20*log10(CenterFrequency*1000);
        %Outdoor loss
        LossOut = max((41 + 20*log10(CenterFrequency/5) + 22.7*log10(Dist_out+Dist_in)), FSL);
        %indoor loss
        LossIn = 0.5*Dist_in; %alpha is 0.5 dB/meter
        %Through wall loss
        LossWall = 14+15*(1-cos(Theta)).^2;
        %Total loss
        loss = LossOut+LossIn+LossWall;

        SF_sigma(NLoSConnectionLinks) = 7;


    case {'B1', 'B2'}         % scenario B1 with d1 & d2 PL model

        if isnan(linkpar.BsHeight(iterpar.UserIndeces(1)))
            BsHeight = 10*ones(1,NumLinks);
            MsHeight = 1.5*ones(1,NumLinks);

        else
            BsHeight = linkpar.BsHeight(iterpar.UserIndeces);
            MsHeight = linkpar.MsHeight(iterpar.UserIndeces);
        end

        if NumLoSConnectionLinks
            MsBsDistance_LoS = MsBsDistance(LoSConnectionLinks);
            Dist1 = MsBsDistance_LoS;
            H_bs = BsHeight(LoSConnectionLinks)-1; %effective environment height
            H_ms = MsHeight(LoSConnectionLinks)-1; %effective environment height
            PL_bp = 4*H_bs.*H_ms*wimpar.CenterFrequency/2.998e8;
            SF_sigma_tmp = [];
            loss_LoS = [];

            ind1 = (Dist1 <= PL_bp & Dist1 >= iterpar.LoS.PL_range(1));
            if sum(ind1)>0
                loss_LoS(ind1) = iterpar.LoS.PL_A(1)*log10(Dist1(ind1)) + iterpar.LoS.PL_B(1) + iterpar.LoS.PL_C*log10(CenterFrequency/5);
            end
            SF_sigma_tmp(ind1) = 3;

            ind2 = (Dist1 >= PL_bp & Dist1 <= iterpar.LoS.PL_range(2));
            if sum(ind2)>0
                loss_LoS(ind2) = iterpar.LoS.PL_A(2)*log10(Dist1(ind2)) + iterpar.LoS.PL_B(2) - 17.3*log10(H_bs(ind2)) - 17.3*log10(H_ms(ind2)) + 2.7*log10(CenterFrequency/5);
            end
            SF_sigma_tmp(ind2) = 3;

            SF_sigma(LoSConnectionLinks) = SF_sigma_tmp;
            loss(LoSConnectionLinks) = loss_LoS;

        end

        if NumNLoSConnectionLinks

            MsBsDistance_NLoS = MsBsDistance(NLoSConnectionLinks);
            H_bs = BsHeight(NLoSConnectionLinks)-1; %effective environment height
            H_ms = MsHeight(NLoSConnectionLinks)-1; %effective environment height
            PL_bp = 4*H_bs.*H_ms*wimpar.CenterFrequency/2.998e8;
            StreetWidth = linkpar.StreetWidth(iterpar.UserIndeces);
            StreetWidth = StreetWidth(NLoSConnectionLinks);
            loss_LoS = [];

            if isnan(linkpar.Dist1) % NaN default -> will be drawn randomly
                Dist1 = 1;
                while any(Dist1>MsBsDistance_NLoS-StreetWidth/2 | Dist1<StreetWidth/2)
                    Dist1 = MsBsDistance_NLoS-StreetWidth/2 - (MsBsDistance_NLoS-StreetWidth).*rand(1,length(MsBsDistance_NLoS));
                end
            else %Dist1 is defined by the user
                Dist1 = linkpar.Dist1(iterpar.UserIndeces);
                Dist1 = Dist1(NLoSConnectionLinks);
                for linkNum = 1:length(MsBsDistance_NLoS) % check applicability and change Dist1 if needed
                    if Dist1(linkNum) > (MsBsDistance_NLoS(linkNum)-StreetWidth(linkNum)/2) || Dist1(linkNum) < StreetWidth(linkNum)/2
                        Dist1(linkNum) = MsBsDistance_NLoS-StreetWidth/2 - (MsBsDistance_NLoS-StreetWidth).*rand(1);
                    end
                end
            end
            Dist2 = MsBsDistance_NLoS - Dist1;

            loss_a = B1_NLOS_PL(Dist1,Dist2,loss_LoS,NLoSConnectionLinks,PL_bp,iterpar,CenterFrequency,wimpar,StreetWidth,H_bs,H_ms);
            loss_b = B1_NLOS_PL(Dist2,Dist1,loss_LoS,NLoSConnectionLinks,PL_bp,iterpar,CenterFrequency,wimpar,StreetWidth,H_bs,H_ms);
            loss(NLoSConnectionLinks) = min(loss_a(NLoSConnectionLinks),loss_b(NLoSConnectionLinks));

            SF_sigma(NLoSConnectionLinks) = 4;
        end


    case {'B3'}

        BsHeight = 6*ones(1,NumLinks);
        MsHeight = 1.5*ones(1,NumLinks);

        if NumLoSConnectionLinks % LOS
            loss(LoSConnectionLinks) = iterpar.LoS.PL_B+iterpar.LoS.PL_A*log10(MsBsDistance(LoSConnectionLinks)) + iterpar.LoS.PL_C*log10(CenterFrequency/5);
            SF_sigma(LoSConnectionLinks) = 3;
        end

        if NumNLoSConnectionLinks % NLOS
            loss(NLoSConnectionLinks) = iterpar.NLoS.PL_B+iterpar.NLoS.PL_A*log10(MsBsDistance(NLoSConnectionLinks)) + iterpar.NLoS.PL_C*log10(CenterFrequency/5);
            SF_sigma(NLoSConnectionLinks) = 4;
        end


    case {'B5a'}
        BsHeight = 25*ones(1,NumLinks);
        MsHeight = 25*ones(1,NumLinks);
        loss = iterpar.LoS.PL_A*log10(MsBsDistance(LoSConnectionLinks))+iterpar.LoS.PL_B + iterpar.LoS.PL_C*log10(CenterFrequency/5);
        SF_sigma(LoSConnectionLinks) = 4;


    case {'B5c'}
        % same as B1 LoS, except for antenna heights
        % Check antenna heights
        if isnan(linkpar.BsHeight(iterpar.UserIndeces(1)))
            BsHeight = 10*ones(1,NumLinks);
            MsHeight = 5*ones(1,NumLinks);
        else
            BsHeight = linkpar.BsHeight(iterpar.UserIndeces);
            MsHeight = linkpar.MsHeight(iterpar.UserIndeces);
        end

        MsBsDistance_LoS = MsBsDistance(LoSConnectionLinks);
        Dist1 = MsBsDistance_LoS;
        H_bs = BsHeight(LoSConnectionLinks)-1; %effective environment height
        H_ms = MsHeight(LoSConnectionLinks)-1; %effective environment height
        PL_bp = 4*H_bs.*H_ms*wimpar.CenterFrequency/2.998e8;
        SF_sigma_tmp = [];
        loss_LoS = [];

        ind1 = (Dist1 <= PL_bp & Dist1 >= iterpar.LoS.PL_range(1));
        if sum(ind1)>0
            loss_LoS(ind1) = iterpar.LoS.PL_A(1)*log10(Dist1(ind1)) + iterpar.LoS.PL_B(1) + iterpar.LoS.PL_C*log10(CenterFrequency/5);
        end
        SF_sigma_tmp(ind1) = 3;

        ind2 = (Dist1 >= PL_bp & Dist1 <= iterpar.LoS.PL_range(2));
        if sum(ind2)>0
            loss_LoS(ind2) = iterpar.LoS.PL_A(2)*log10(Dist1(ind2)) + iterpar.LoS.PL_B(2) - 17.3*log10(H_bs(ind2)) - 17.3*log10(H_ms(ind2)) + 2.7*log10(CenterFrequency/5);
        end
        SF_sigma_tmp(ind2) = 3;

        SF_sigma(LoSConnectionLinks) = SF_sigma_tmp;
        loss(LoSConnectionLinks) = loss_LoS;


    case{'B5f'}     % NLOS
        BsHeight = 25*ones(1,NumLinks);
        MsHeight = 15*ones(1,NumLinks);
        loss = iterpar.NLoS.PL_A*log10(MsBsDistance) + iterpar.NLoS.PL_B + iterpar.NLoS.PL_C*log10(CenterFrequency/5);
        SF_sigma(NLoSConnectionLinks) = 8;


    case {'C1'}
        %Check antenna heights
        if isnan(linkpar.BsHeight(iterpar.UserIndeces(1)))
            BsHeight = 25*ones(1,NumLinks);
            MsHeight = 1.5*ones(1,NumLinks);
        else
            BsHeight = linkpar.BsHeight(iterpar.UserIndeces);
            MsHeight = linkpar.MsHeight(iterpar.UserIndeces);
        end

        PL_bp = 4*BsHeight.*MsHeight*wimpar.CenterFrequency/2.998e8;

        if NumLoSConnectionLinks %LOS
            MsBsDistance_LoS = MsBsDistance(LoSConnectionLinks);
            BsHeight_LoS = BsHeight(LoSConnectionLinks);
            MsHeight_LoS = MsHeight(LoSConnectionLinks);
            SF_sigma_tmp = [];
            loss_LoS = [];

            ind1 = (MsBsDistance_LoS <= PL_bp(LoSConnectionLinks) & MsBsDistance_LoS >= iterpar.LoS.PL_range(1));
            if sum(ind1)>0
                loss_LoS(ind1) = iterpar.LoS.PL_A(1)*log10(MsBsDistance_LoS(ind1)) + iterpar.LoS.PL_B(1) + iterpar.LoS.PL_C(1)*log10(CenterFrequency/5);
            end
            SF_sigma_tmp(ind1) = 4;
            ind2 = MsBsDistance_LoS >= PL_bp(LoSConnectionLinks) & MsBsDistance_LoS <= iterpar.LoS.PL_range(2);
            if sum(ind2)>0
                loss_LoS(ind2) = iterpar.LoS.PL_A(2)*log10(MsBsDistance_LoS(ind2)) + iterpar.LoS.PL_B(2) + iterpar.LoS.PL_C(2)*log10(CenterFrequency/5) - 16.2*log10(MsHeight_LoS(ind2)) - 16.2*log10(BsHeight_LoS(ind2));
            end
            SF_sigma_tmp(ind1) = 6;

            loss(LoSConnectionLinks) = loss_LoS;
            SF_sigma(LoSConnectionLinks)=SF_sigma_tmp;
        end

        if NumNLoSConnectionLinks %NLOS
            loss(NLoSConnectionLinks) = (44.9-6.55*log10(BsHeight(NLoSConnectionLinks))).*log10(MsBsDistance(NLoSConnectionLinks)) + 31.46 + 5.83*log10(BsHeight(NLoSConnectionLinks)) + iterpar.NLoS.PL_C*log10(CenterFrequency/5);
            SF_sigma(NLoSConnectionLinks) = 8;
        end


    case {'C2', 'C3'}
        %Check antenna heights
        if isnan(linkpar.BsHeight(iterpar.UserIndeces(1)))
            BsHeight = 25*ones(1,NumLinks);
            MsHeight = 1.5*ones(1,NumLinks);
        else
            BsHeight = linkpar.BsHeight(iterpar.UserIndeces);
            MsHeight = linkpar.MsHeight(iterpar.UserIndeces);
        end
        
        PL_bp = 4*BsHeight.*MsHeight*wimpar.CenterFrequency/2.998e8;

        if NumLoSConnectionLinks %LOS
            MsBsDistance_LoS = MsBsDistance(LoSConnectionLinks);
            BsHeight_LoS = BsHeight(LoSConnectionLinks);
            MsHeight_LoS = MsHeight(LoSConnectionLinks);
            SF_sigma_tmp = [];
            loss_LoS = [];

            ind1 = (MsBsDistance_LoS <= PL_bp(LoSConnectionLinks) & MsBsDistance_LoS >= iterpar.LoS.PL_range(1));
            if sum(ind1)>0
                loss_LoS(ind1) = iterpar.LoS.PL_A(1)*log10(MsBsDistance_LoS(ind1)) + iterpar.LoS.PL_B(1) + iterpar.LoS.PL_C(1)*log10(CenterFrequency/5);
            end
            SF_sigma_tmp(ind1) = 4;
            ind2 = MsBsDistance_LoS >= PL_bp(LoSConnectionLinks) & MsBsDistance_LoS <= iterpar.LoS.PL_range(2);
            if sum(ind2)>0
                loss_LoS(ind2) = iterpar.LoS.PL_A(2)*log10(MsBsDistance_LoS(ind2)) + iterpar.LoS.PL_B(2) + iterpar.LoS.PL_C(2)*log10(CenterFrequency/5) - 14*log10(MsHeight_LoS(ind2)) - 14*log10(BsHeight_LoS(ind2));
            end
            SF_sigma_tmp(ind1) = 6;

            loss(LoSConnectionLinks) = loss_LoS;
            SF_sigma(LoSConnectionLinks) = SF_sigma_tmp;
        end

        if NumNLoSConnectionLinks %NLOS
            loss(NLoSConnectionLinks) = (44.9-6.55*log10(BsHeight(NLoSConnectionLinks))).*log10(MsBsDistance(NLoSConnectionLinks)) + 34.46 + 5.83*log10(BsHeight(NLoSConnectionLinks)) + iterpar.NLoS.PL_C*log10(CenterFrequency/5);
            SF_sigma(NLoSConnectionLinks) = 8;
        end

        
        
    case {'C4'}
        % Check antenna heights
        if isnan(linkpar.BsHeight(iterpar.UserIndeces(1)))
            BsHeight = 10*ones(1,NumLinks);
            MsHeight = (1.5 + 3*(NumFloors-1)).*ones(1,NumLinks);
        else
            BsHeight = linkpar.BsHeight(iterpar.UserIndeces);
            MsHeight = linkpar.MsHeight(iterpar.UserIndeces);
        end

        StreetWidth = linkpar.StreetWidth(iterpar.UserIndeces);
        if isnan(linkpar.Dist1) % NaN default -> will be drawn randomly
            Dist_out = 1;
            while any(Dist_out > MsBsDistance | Dist_out < StreetWidth/2)
                Dist_out = MsBsDistance - (MsBsDistance-StreetWidth/2).*rand(1,length(MsBsDistance));
            end

        else %Dist1 is defined by the user
            Dist_out = linkpar.Dist1(iterpar.UserIndeces);
            Dist_out = Dist_out(NLoSConnectionLinks);
            for linkNum = 1:length(MsBsDistance) % check applicability and change Dist1 if needed
                if Dist_out(linkNum) > MsBsDistance(linkNum) || Dist_out(linkNum) < StreetWidth(linkNum)/2
                    Dist_out(linkNum) = MsBsDistance(linkNum) - (MsBsDistance(linkNum)-StreetWidth(linkNum)/2).*rand(1);
                end
            end
        end
        Dist_in = MsBsDistance - Dist_out; %indoor distance
        
        PL_C2_NLOS = (44.9-6.55*log10(BsHeight)).*log10(Dist_out+Dist_in) + 31.46 + 5.83*log10(BsHeight) + iterpar.NLoS.PL_C*log10(CenterFrequency/5);
        loss = PL_C2_NLOS + 17.4 + 0.5*Dist_in;
        SF_sigma(NLoSConnectionLinks) = 8;

        

    case {'D1'}
        %Check antenna heights
        if isnan(linkpar.BsHeight(iterpar.UserIndeces(1)))
            BsHeight = 32*ones(1,NumLinks);
            MsHeight = 1.5*ones(1,NumLinks);
        else
            BsHeight = linkpar.BsHeight(iterpar.UserIndeces);
            MsHeight = linkpar.MsHeight(iterpar.UserIndeces);
        end

        PL_bp = 4*BsHeight.*MsHeight*wimpar.CenterFrequency/2.998e8;

        if NumLoSConnectionLinks %LOS
            MsBsDistance_LoS = MsBsDistance(LoSConnectionLinks);
            BsHeight_LoS = BsHeight(LoSConnectionLinks);
            MsHeight_LoS = MsHeight(LoSConnectionLinks);
            SF_sigma_tmp = [];
            loss_LoS = [];

            ind1 = (MsBsDistance_LoS <= PL_bp(LoSConnectionLinks) & MsBsDistance_LoS >= iterpar.LoS.PL_range(1));
            if sum(ind1)>0
                loss_LoS(ind1) = iterpar.LoS.PL_A(1)*log10(MsBsDistance_LoS(ind1)) + iterpar.LoS.PL_B(1) + iterpar.LoS.PL_C(1)*log10(CenterFrequency/5);
            end
            SF_sigma_tmp(ind1) = 4;
            ind2 = MsBsDistance_LoS >= PL_bp(LoSConnectionLinks) & MsBsDistance_LoS <= iterpar.LoS.PL_range(2);
            if sum(ind2)>0
                loss_LoS(ind2) = iterpar.LoS.PL_A(2).*log10(MsBsDistance_LoS(ind2)) + iterpar.LoS.PL_B(2) - 18.5*log10(BsHeight_LoS(ind2)) - 18.5*log10(MsHeight_LoS(ind2))+ iterpar.LoS.PL_C(2)*log10(CenterFrequency/5);
            end
            SF_sigma_tmp(ind1) = 6;

            loss(LoSConnectionLinks) = loss_LoS;
            SF_sigma(LoSConnectionLinks) = SF_sigma_tmp;
        end

        if NumNLoSConnectionLinks % NLOS
            loss(NLoSConnectionLinks) =  iterpar.NLoS.PL_A*log10(MsBsDistance(NLoSConnectionLinks)) + iterpar.NLoS.PL_B - 0.13*(BsHeight(NLoSConnectionLinks)-25).*log10(MsBsDistance(NLoSConnectionLinks)/100) - 0.9*(MsHeight(NLoSConnectionLinks)-1.5) + iterpar.NLoS.PL_C*log10(CenterFrequency/5);

            SF_sigma(NLoSConnectionLinks) = 8;
        end


    case {'D2a'}
        %Check antenna heights
        if isnan(linkpar.BsHeight(iterpar.UserIndeces(1)))
            BsHeight = 32*ones(1,NumLinks);
            MsHeight = 1.5*ones(1,NumLinks);
        else
            BsHeight = linkpar.BsHeight(iterpar.UserIndeces);
            MsHeight = linkpar.MsHeight(iterpar.UserIndeces);
        end

        PL_bp = 4*BsHeight.*MsHeight*wimpar.CenterFrequency/2.998e8;


        MsBsDistance_LoS = MsBsDistance(LoSConnectionLinks);
        BsHeight_LoS = BsHeight(LoSConnectionLinks);
        MsHeight_LoS = MsHeight(LoSConnectionLinks);
        SF_sigma_tmp = [];
        loss_LoS = [];

        ind1 = (MsBsDistance_LoS <= PL_bp(LoSConnectionLinks) & MsBsDistance_LoS >= iterpar.LoS.PL_range(1));
        if sum(ind1)>0
            loss_LoS(ind1) = iterpar.LoS.PL_A(1)*log10(MsBsDistance_LoS(ind1)) + iterpar.LoS.PL_B(1) + iterpar.LoS.PL_C(1)*log10(CenterFrequency/5);
        end
        SF_sigma_tmp(ind1) = 4;
        ind2 = MsBsDistance_LoS >= PL_bp(LoSConnectionLinks) & MsBsDistance_LoS <= iterpar.LoS.PL_range(2);
        if sum(ind2)>0
            loss_LoS(ind2) = iterpar.LoS.PL_A(2).*log10(MsBsDistance_LoS(ind2)) + iterpar.LoS.PL_B(2) - 18.5*log10(BsHeight_LoS(ind2)) - 18.5*log10(MsHeight_LoS(ind2))+ iterpar.LoS.PL_C(2)*log10(CenterFrequency/5);
        end
        SF_sigma_tmp(ind1) = 6;

        loss(LoSConnectionLinks) = loss_LoS;
        SF_sigma(LoSConnectionLinks) = SF_sigma_tmp;

        

    otherwise       % all other scenarios with one d PL model
        
        BsHeight = 10*ones(1,NumLinks);
        MsHeight = 1.5*ones(1,NumLinks);
        SF_sigma = 3*ones(1,NumLinks);
        loss(LoSConnectionLinks) = iterpar.LoS.PL_B + iterpar.LoS.PL_A*log10(MsBsDistance(LoSConnectionLinks));
        loss(NLoSConnectionLinks) = iterpar.NLoS.PL_B + iterpar.NLoS.PL_A*log10(MsBsDistance(NLoSConnectionLinks));

end     % end switch Scenario


% output
linkpar.MsHeight(iterpar.UserIndeces) = MsHeight;
linkpar.BsHeight(iterpar.UserIndeces) = BsHeight;
iterpar.NLoS.SF_sigma=SF_sigma(1);
loss=loss(:);

% function d=distrnd(num,rmax)
% % DISTRND Distance from BS in a circular cell
% %   D=DISTRND(K,RMAX) generates K random variables from the pdf
% %   p(r)=2*r/RMAX^2. This is the pdf of distance from base station when
% %   users are uniformly (in area) distributed in a cell with radius RMAX
% %   1 x num vector.
%
% %   Authors: Jari Salo (HUT), Marko Milojevic (TUI)
%
% % create random variables from triangular pdf whose width is 2*rmax
% a=sum(repmat(rmax,2,1).*rand(2,num));
%
% % fold the random variables about the rmax
% inds=find(a>rmax);
% a(inds)=-a(inds)+2*rmax(inds);
%
% d=a(:).';


function loss = B1_NLOS_PL(Dist1,Dist2,loss_LoS,NLoSConnectionLinks,PL_bp,iterpar,CenterFrequency,wimpar,StreetWidth,H_bs,H_ms);


% Same as B1 LOS
ind1 = (Dist1 <= PL_bp & Dist1 >= (iterpar.NLoS.PL_range(1)-StreetWidth/2));
if sum(ind1)>0
    loss_LoS(ind1) = iterpar.NLoS.PL_A(1)*log10(Dist1(ind1)) + iterpar.NLoS.PL_B(1) + 20*log10(CenterFrequency/5);
end

ind2 = (Dist1 >= PL_bp & Dist1 <= iterpar.NLoS.PL_range(2)-StreetWidth/2);
if sum(ind2)>0
    loss_LoS(ind2) = iterpar.NLoS.PL_A(2)*log10(Dist1(ind2)) + iterpar.NLoS.PL_B(2) - 17.3*log10(H_bs(ind2)) - 17.3*log10(H_ms(ind2)) + 2.7*log10(CenterFrequency/5);
end

% plus component for B1 NLOS
nj = max(2.8-0.0024*Dist1 , 1.84);
loss(NLoSConnectionLinks) = loss_LoS + 20 - 12.5*nj + 10*nj.*log10(Dist2) + 3*log10(CenterFrequency/5);
