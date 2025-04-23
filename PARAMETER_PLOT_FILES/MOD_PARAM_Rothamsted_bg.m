%%%%%%%%%%%%%%%%%%% PARAMETERS AND INITIAL CONDITION %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% - only for biogeochemistry part - %%%%%%%%%%%%%%%%%%%%%
% - data from initial spinup must already be loaded in struct `prev_data` - 
%%%%%%%%%%%%% - data from bg spinup in struct `prev_bg` - %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: High vegetation ignored here!

% TODO: Anything that gets overwritten in MAIN_FRAME before calling this file should be defined here properly
%       Consider calling `MOD_PARAM_Rothamsted_Spinup_common.m` here, and loading only relevant data from previous spinups
Ci_sunL(1,:) = [Ca(1)]; % %% [umolCO2/mol]
Ci_sunH(1,:) = [Ca(1)]; %% [umolCO2/mol]
Ci_shdL(1,:) = [Ca(1)]; % %% [umolCO2/mol]
Ci_shdH(1,:) = [Ca(1)]; %% [umolCO2/mol]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
LAI_L(1,:)=prev_data.LAI_L(end,:);
B_L(1,:,:)= prev_data.B_L(end,:,:);
Rrootl_L(1,:)=prev_data.Rrootl_L(end,:);
PHE_S_L(1,:)=prev_data.PHE_S_L(end,:);
dflo_L(1,:)=prev_data.dflo_L(1,:);
AgeL_L(1,:)=prev_data.AgeL_L(end,:);
e_rel_L(1,:)=prev_data.e_rel_L(end,:);
hc_L(1,:) =prev_data.hc_L(end,:);
SAI_L(1,:)=prev_data.SAI_L(end,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BLit(1)=prev_data.BLit(end);  %% [kg DM /m2 PFT] Litter Biomass



if OPT_SoilBiogeochemistry == 1  % Values from spinup
    Nreserve_L(1,:) = 0.547872612893958; %
    Preserve_L(1,:) = 0.057805595467867 ; %
    Kreserve_L(1,:) = 0.036037541469013 ; %
    
    FNC_L(1,:) = 1;
    NupI_L(1,:,:) = repmat([0.012738244408113 0.001046744728782  0.001927354959625], 1, 1);
    RexmyI(1,:,:) = [0.000072416067436 0.005891840450278  0.000231813245917];
    NavlI(1,:,:) =[  3.929546067507789   2.521608870474120  0.172079615252432];% guess the two above are not crop specific ?

    %% From spinup
    P(1,:) = load(['..' filesep 'results' filesep 'Rothamsted_biogeo_spinup-B_end-curr.mat']).B_end;
    P(P<0) = 0;

    %% Previous
    % P(1,:)=   1000*[ 0.124083987036102   0.325573628166407   0.054069975887170                   0                   0   0.101512860630231   0.126564601766440,...
    %     0.014062733529605   0.755595303823738   1.096103576915454   6.402023138646436   0.007894696102488   0.013135775003470   0.000162321193831,...
    %     0.000104278919687   0.000081160596916   0.000156418379530   0.038738455769914   0.133541727544331   0.065957977030849                   0,...
    %     0.002528243856056   0.009999140372626                   0   0.004492687719394   0.913446366042410   0.005555237467077   0.011170245415232,...
    %     0.003668130264076                   0   0.000185349798914   0.000017774714467   0.000036321130404   0.000252824385605   0.001098054538178,...
    %     0   0.000476491894740   0.137001078180096   0.001107194360735   0.001675533186330   0.000550219539611                   0,...
    %     0.000007588306019   0.150000000000000   0.000172964791821   0.015000467306710   0.000000368010905   0.005528056261498                   0,...
    %     0.000603437980692   0.016621912220522   0.000015673158772   0.000009163307352   0.000808451498098   0.500387085728242];

    Upl=0.01; %% Soil production [mm/yr]
    HIST=0;
    %[B_IO]=Biogeochemistry_IO(Zs(ms+1)*1000,Lat,Lon,Upl,HIST,FertN,FertP,FertK,DepN,DepP,DepK);
    [B_IO]=Biogeochemistry_IO(Zs(ms+1),Lat,Lon,Upl);
    B_IO.SC_par=[1 1 1 1];
    PHs=7;
end
%%%%%%%%%%%%%%%%%
Vx_H=0;  %% [mm/ m2 PFT];
Vl_H=0;  %% [mm/ m2 PFT];
Vx_L=0;   %% [mm/ m2 PFT];
Vl_L=100;   %% [mm/ m2 PFT];  % From Buckley et al., 2024 
%%%
TBio_L=1;  %%[ton DM / ha ]
TBio_H=0;  %[ton DM / ha ]
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%5%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Initial Conditions
SWE(1)=0; %% [mm]
SND(1)=0;
Ts(1)=Ta(1)+2;
Tdamp(1)=3;   % From Buckley et al., 2024 
Tdp(1,:)= 3*ones(1,ms);  % From Buckley et al., 2024 
Sdp(1,:)= (0.001*dz).*cv_s.*(Tdp(1,:)); %% [J �C/m^2 K]
TdpI_H(1,:)=3;  % From Buckley et al., 2024 
TdpI_L(1,:)=3;  % From Buckley et al., 2024 
%%% Snow_alb = soil_alb initial 
snow_alb.dir_vis = 0.2;
snow_alb.dif_vis = 0.2;
snow_alb.dir_nir = 0.2;
snow_alb.dif_nir = 0.2;
In_L(1,:)=0; In_H(1,:)=0;
In_urb(1)=0; In_rock(1)= 0;
SP_wc(1)=0 ; %%[mm]
In_Litter(1)=0;
In_SWE(1)= 0;
ros(1)= 0;
t_sls(1)= 0;
e_sno(1) = 0.97;
tau_sno(1) = 0;
EK(1)=0; 
WAT(1) = 0;
ICE(1) = 0; 
IP_wc(1)=0; 
ICE_D(1)= 0;
FROCK(1)=0;
Ws_under(1)=1; 
%%%%%%%%%%%%%% Volume [mm]
%O(1,:)=  [0.3269    0.3280    0.3290    0.3314    0.3341    0.3367    0.3402    0.3442    0.3476    0.3505    0.3541  0.3564]; 
O(1,:)=  Ofc;  
%dz(ms)=ZWT(1)-Zs(ms);
%%%%%%%%%%%%%%%%%%%%%%%%
%Vmin= 0*dz;
%Vmax= (Osat-Ohy)*dz;
%Vmax(ms) = Vmax(ms) + (Zbed-ZWT(1))*(Osat-Ohy);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V(1,:) = (O(1,:) -Ohy).*dz;
%V(1,ms)= 0 + (Osat-Ohy)*(Zbed-ZWT(1));
%%%%%%
cd(cur_dir)