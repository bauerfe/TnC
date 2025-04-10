%%%%%%%%%%%%%%%%%%% PARAMETERS AND INITIAL CONDITION %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SOIL AND HYDROLOGICAL PARAMETER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
%%% Rainfall Disaggregation

Vmax_L = [96]; %35



a_dis = NaN ;
pow_dis = NaN;
%%%%%%%%%%%%%%%%%%%%
fpr=1;
SvF=1; %% Sky View Factor
SN=0; %% Stream Identifier
Slo_top=0;  %% [fraction dy/dx]
Slo_pot=zeros(1,ms); %% [fraction dy/dx]
Asur = 1./cos(atan(Slo_top)); %% Real Area/Projected Area [m^2/m^2]
Ared = 1; %% 1-Frock
aR =1; %%% anisotropy ratio
%Kh=Ks*aR;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
cellsize=1; %%[m^2];
aTop = 1000*cellsize^2./cellsize; %% [mm] Ratio betweeen Area/ContourLenght
Kbot = NaN; %% [mm/h] Conductivity at the bedrock layer % at NaN - 
Krock = NaN; %% [mm/h] Conductivity of Fractured Rock
zatm = 2.41; %% Reference Height
%%%%%%%%%%%%%%%%%%
%%%%%%% VEG. SPECIES  --- Low Grasses -- High Decidous
%%%% LAND COVER PARTITION
Cwat = 0.0; Curb = 0.0 ; Crock = 0.0;
Cbare = 0.0; Ccrown = [1.0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SOIL INPUT %%% Sandy Loam
Pcla= 0.29;
Psan= 0.35; % 0.254 old value before sept  %% 35/40/25 sand silt clay is what they have on the paper 
Porg= 0.085;
Color_Class = 0;
%%%%%%%%%%%%%%%%%%%
[Osat,L,Pe,Ks,O33,rsd,lan_dry,lan_s,cv_s,K_usle]=Soil_parameters(Psan,Pcla,Porg);
%%%%%%%%%%%%
rsd=rsd*ones(1,ms);
lan_dry=lan_dry*ones(1,ms);
lan_s =lan_s*ones(1,ms);
cv_s = cv_s*ones(1,ms);
%%%%%%%%%%%%%%%
SPAR=2; %%% SOIL PARAMETER TYPE 1-VanGenuchten 2-Saxton-Rawls
%nVG=L+1;
%alpVG = 1/(-101.9368*Pe); %%[1/mm]%;
p=3+2/L;
m=2/(p-1); nVG= 1/(1-m);
alpVG=(((-101.9368*Pe)*(2*p*(p-1))/(p+3))*((55.6+7.4*p+p^2)/(147.8+8.1*p+0.092*p^2)))^-1; %%[1/mm]%;
%%%
Osat=Osat*ones(1,ms);
%Ohy = Ohy*ones(1,ms) ; %% [-]
L=L*ones(1,ms);
Pe = Pe*ones(1,ms);
O33 = O33*ones(1,ms);
alpVG= alpVG*ones(1,ms); %% [1/mm]
nVG= nVG*ones(1,ms); %% [-]
Ks_Zs= Ks*ones(1,ms); %%[mm/h]
%%%%%%%%%%% Matric Potential
Kfc = 0.2; %% [mm/h]
Phy = 10000; %% [kPa]
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
[Ofc,Oss_Lp,Owp_Lp,Ohy]=Soil_parametersII(ms,Osat,L,Pe,Ks_Zs,O33,nVG,alpVG,Kfc,1,1,Phy);
%%%%%%%%%%%%%%
clear Oss_Lp Owp_Lp
%%%%%%%%%%%%%%
Oice = 0;
%%%%%%%%%%%%%
Zs = [ 0 10 20 50 100 150 200 300 400 500 600 800 1000 1300]; %%% [ms+1]
Zdes = 10;
Zinf=  10;
Zbio = 250;
if  not(length(Zs)==ms+1)
    disp('SOIL LAYER MESH INCONSISTENT')
    return
end
[EvL_Zs]=Evaporation_layers(Zs,Zdes); %%% Evaporation Layer fraction
[Inf_Zs]=Evaporation_layers(Zs,Zinf); %%% Infiltration Depth Layer fraction
[Bio_Zs]=Evaporation_layers(Zs,Zbio); %%% Infiltration Depth Layer fraction
dz= diff(Zs); %%%% [mm]  Thickness of the Layers
Dz=zeros(1,ms);
for i = 1:ms
    if i>1
        Dz(i)= (dz(i)+ dz(i-1))/2; %%% Delta Depth Between Middle Layer  [mm]
    else
        Dz(i)=dz(1)/2; %%% Delta Depth Between First Middle Layer and soil surface [mm]
    end
end
%%%%%%%%%%%%%%%%% OTHER PARAMETER
In_max_urb=5;
In_max_rock=0.1; %% [mm]
%%%%%%%%%%%%% SNOW PARAMETER
TminS=-0.8;%% Threshold temperature snow
TmaxS= 2.8;%% Threshold temperature snow
ros_max1=580; %600; %%% [kg/m^3]
ros_max2=300; %450; %%% [kg/m^3]
Th_Pr_sno = 8.0; %%% [mm/day] Threshold Intensity of snow to consider a New SnowFall
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ICE Parameter
Ice_wc_sp =0.01; %% [-] Specific Maximum water content ice
ros_Ice_thr = 500 ; %% [kg/m^3] Density Thrshold to transform snow into ice
Aice = 0.28; %% [-] Ice albedo
WatFreez_Th = -8; %% [°C] Threshold for freezing lake water
dz_ice = 0.45; %% [mm / h] Water Freezing Layer progression without snow-layer
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%
ExEM = 0.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% cc -- number of crown area
%%% Root Depth
CASE_ROOT=1;  %%% Type of Root Profile
%%%%%
ZR95_H = [0]; %% [mm]
ZR95_L = [250]; %% [mm] 
ZR50_H = NaN;
ZR50_L = NaN;
ZRmax_H = NaN;
ZRmax_L = NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kct=0.75; %%% Factor Vegetation Cover --- for throughfall
%5 Interception Parameter
gcI=3.7; %%% [1/mm]
KcI=0.06; %%%% [mm] -- Mahfouf and Jacquemin 1989
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Interception Parameter
Sp_SN_In= 5.9; %% [mm/LAI]
Sp_LAI_L_In= [0.2]; %%[mm/LAI]
Sp_LAI_H_In= [0.2]; %%[mm/LAI]
%%%%%%%%%%% Leaf Dimension
d_leaf_H= [3.5]; %%[cm]
d_leaf_L= [0.8];  %% [cm]
%%%%%%%% Biochemical parameter
KnitH=0.2; %%% Canopy Nitrogen Decay
KnitL=0.15; %%% Canopy Nitrogen Decay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mSl_H = 0.0;%% [m2 PFT /gC]  Linear increase in Sla with LAI
mSl_L = 0.0; % 0.001; %% [m2 LAI/gC * m2 PFT / m2 LAI]  0.0 - 0.004  Brod. Dec. Tree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  Photosynthesis Parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------
FI_H=[0.081];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_H=[1000]; %%[Pa]
a1_H=[7];
go_H=[0.01];% [mol / s m^2] minimum Stomatal Conductance
CT_H=[3]; %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_H =[0.649];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_H =[72]; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_H=[Inf]; %% [mol CO2 / s m^2 ];  mesophyll conductance
rjv_H=[1.97]; %%% Scaling Jmax - Vmax  [umol electrons / umolCO2 ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------
FI_L=[0.081];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_L=[1000]; %%[Pa]
a1_L=[6];
go_L=[0.01];% % [mol / s m^2] minimum Stomatal Conductance
CT_L=[3];  %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_L =[0.656];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_L =[55]; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_L=[Inf];
rjv_L= [2.4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pss_H = [800]; Pwp_H = [3000]; %%% [kPa]
%Pss_L = [500]; Pwp_L = [3500]; %%% [kPa]
Psi_sto_50_H =  -2.0 ;%% [MPa]  Water Potential at 50% loss conductivity
Psi_sto_00_H =  -0.5; %% [MPa]  Water Potential at 2% loss conductivity
%%% Leaf
PsiL00_H =  -2.7 ;%%[MPa]  Water Potential at 50% loss conductivity
PsiL50_H =  -5.6; %% [MPa]  Water Potential at 2% loss conductivity
Kleaf_max_H = 5 ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_H  = 1200;  %%%  [500 - 3000]%  Leaf capacitance [mmolH20 / m^2 leaf MPa]
%%% Xylem
Axyl_H = 15.0 ; %% [cm^2 stem /m^2 PFT]
Kx_max_H = 80000;  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_H = -3.5; %%[MPa]  Water Potential at 50% loss conductivity
Cx_H= 150; %%% [kg / m^3 sapwood MPa]
%%------------------------
%
Psi_sto_00_L = -0.5;%  %% [MPa]  Water Potential at PLCs loss conductivity
Psi_sto_50_L = -2.0;%  %% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_L =  -0.9; %% [MPa]  Water Potential at PLCs% loss conductivity
PsiL50_L =  -3.0 ;%%[MPa]  Water Potential at 50% loss conductivity % reduce this to fix rapid drying
Kleaf_max_L = 5 ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_L  = 1200;  %%%  [500 - 3000]%  [mmolH20 / m^2 leaf MPa]
%%% Xylem
Axyl_L = 0.0 ; %% [cm^2 stem /m^2 PFT]
Kx_max_L = 80000;  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_L = -4.5; %%[MPa]  Water Potential at 50% loss conductivity
Cx_L= 150; %%% [kg / m^3 sapwood MPa]
%%%%%%%%%%%%%%%% Root Parameters
[RfH_Zs,RfL_Zs]=Root_Fraction_General(Zs,CASE_ROOT,ZR95_H,ZR50_H,ZR95_L,ZR50_L,ZRmax_H,ZRmax_L);
%%%% Growth Parameters
PsiG50_H= -0.45;  %%[MPa]
PsiG99_H= -1.2;  %%[MPa]
gcoef_H = 3.5; % [gC/m2 day]
%%------
PsiG50_L= -3.0;
PsiG99_L= -4.0;
gcoef_L = 3.5; % [gC/m2 day]

%%%%%%%% Vegetation Optical Parameter
[PFT_opt_H(1)]=Veg_Optical_Parameter(0);
[PFT_opt_L(1)]=Veg_Optical_Parameter(13);
OM_H=1;
OM_L=1;
%%%%%%
Sllit = 2 ; %%% Litter Specific Leaf area [m2 Litter / kg DM]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% VEGETATION PART %%%%%%
%%% HIGH VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%
Sl_H = [0.016]; % 0.018 0.05 -0.005 [m^2 gC] specific leaf area of  biomass [m^2 /gC]
Nl_H= [30]; %[gC/gN ] Leaf Carbon-Nitrogen ratio
[Stoich_H(1)]=Veg_Stoichiometric_Parameter(Nl_H(1));
%PLNR_H = 0.033; %% Percentage of Leaf N in Rubisco  [kgNRubisco / kgN]
r_H = [0.030];  %% [0.066 -0.011]respiration rate at 10° [gC/gN d ]
gR_H= [0.25]; % [0.22 - 0.28] growth respiration  [] -- [Rg/(GPP-Rm)]
%LAI_max_H= [2.0 10.0];
aSE_H= [1]; %%% Plant Type -- 1 Seasonal Plant --  0 Evergreen  -- 2 Grass species -- 3 Crops
dd_max_H= [1/365]; %%%0.005  [1/d]  0.0250 -- 0.005-0.025 death maximum for drought
dc_C_H =  2/365; %% [1/ d°C] -- [Factor of increasing mortality]
Tcold_H = [7]; %% [°C] Cold Leaf Shed
drn_H=  [1/1095]; %% turnover root  [1/d]
dsn_H= [1/365]; % normal transfer rate sapwood [1/d]
age_cr_H= [150]; %% [day] Critical Leaf Age
Bfac_lo_H= [0.95]; %% Leaf Onset Water Stress
Bfac_ls_H= [NaN]; %% Leaf Shed Water Stress [0-1]
Tlo_H = [12.9]; %% Mean Temperature for Leaf onset
Tls_H = [NaN NaN]; %% Mean Temperature for Leaf Shed
PAR_th_H= [NaN]; 
dmg_H= [35]; %%% Tree 30 Grasses Day of Max Growth
LAI_min_H = [0.01];
Trr_H = [3.5]; %% Translocation rate [gC /m^2 d]
mjDay_H = [180]; %% Maximum Julian day for leaf onset
LDay_min_H =[11.0]; %% Minimum Day duration for leaf onset
LtR_H = [1.0]; %%% Leaf to Root ratio maximum
Mf_H= [1/50]; %% fruit maturation turnover [1/d]
Wm_H= [1/16425] ; % wood turnover coefficient [1/d]
eps_ac_H = [1]; %% Allocation to reserve parameter [0-1]
LDay_cr_H = [12.30]; %%%  Threshold for senescence day light [h]
Klf_H =[1/15]; %% Dead Leaves fall turnover [1/d]
fab_H = [0.74]; %% fraction above-ground sapwood and reserve
fbe_H = [0.26]; %% fraction below-ground sapwood and reserve
ff_r_H= [0.1]; %% Reference allocation to Fruit and reproduction
[ParEx_H(1)]=Exudation_Parameter(0); 
[Mpar_H(1)]=Vegetation_Management_Parameter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOW VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sl_L = [0.035]; % 0.05 -0.005 [m^2 gC] specific leaf area of  biomass [m^2 /gC]
Nl_L= [23]; %[gC/gN ] Leaf Carbon-Nitrogen ratio
[Stoich_L(1)]=Veg_Stoichiometric_Parameter(Nl_L(1));
%PLNR_L = 0.033; %% Percentage of Leaf N in Rubisco  [kgNRubisco / kgN]
r_L = [0.060];  %% [0.066 -0.011]respiration rate at 10° [gC/gN d ]
gR_L= [0.22]; % [0.22 - 0.28] growth respiration  [] -- [Rg/(GPP-Rm)]
%LAI_max_L= [2.0 10.0];2
aSE_L= [2]; %%% Plant Type -- 1 Seasonal Plant --  0 Evergreen  -- 2 Grass species -- 3 Crops
dd_max_L= [1/45];%%%0.005  [1/d]  0.0250 -- 0.005-0.025 death maximum for drought
dc_C_L = [7/365]; %% [1/ d°C] -- [Factor of increasing mortality]
Tcold_L = [-2.0]; %% [°C] Cold Leaf SLed
drn_L=  [1/450]; %% turnover root  [1/d]
dsn_L= [ 1/365]; % normal transfer rate sapwood [1/d]
age_cr_L= [180]; %% [day] Critical Leaf Age
Bfac_lo_L= [0.99]; %% Leaf Onset Water Stress
Bfac_ls_L= [NaN]; %% Leaf Shed Water Stress [0-1]
Tlo_L = [0.0]; %% Mean Temperature for Leaf onset
Tls_L = [NaN]; %% Mean Temperature for Leaf Shed
PAR_th_L= [NaN]; 
dmg_L= [20]; %%% Tree 30 Grasses Day of Max Growth
LAI_min_L = [0.1]; % was 0.1
Trr_L = [2.0]; %% Translocation rate [gC /m^2 d]
mjDay_L = [250]; %% Maximum Julian day for leaf onset
LDay_min_L =[9.7]; %% Minimum Day duration for leaf onset
LtR_L = [0.35]; %%% Leaf to Root ratio maximum
Mf_L= [1/50]; %% fruit maturation turnover [1/d]
Wm_L= [0] ; % wood turnover coefficient [1/d]
eps_ac_L = [0.2]; %% Allocation to reserve parameter [0-1]
LDay_cr_L = [10.7]; %%%  Threshold for senescence day light [h]
Klf_L =  1/50 ; % [1/83]; %% Dead Leaves fall turnover [1/d]
fab_L = 0.0; %% fraction above-ground sapwood and reserve
fbe_L = 1.0; %% fraction below-ground sapwood and reserve
ff_r_L= [0.1];
[ParEx_L(1)]=Exudation_Parameter(0); 
[Mpar_L(1)]=Vegetation_Management_Parameter;
% Mpar_L(1).jDay_cut =[ datenum([2006 6 1]) datenum([2006 8 1])];
Mpar_L(1).LAI_cut =[1.3]; %% LAI of grass after cut %% 2.65 Height of 15cm more or less - was 1.68


Mpar_L(1).jDay_cut = [ ...
    datenum([1990 5 4]) datenum([1991 5 4]) datenum([1992 5 4]) datenum([1993 5 4])...
    datenum([2005 5 1]) datenum([2005 6 28]) datenum([2005 7 27]) datenum([2005 7 23]) datenum([2005 8 31]) datenum([2005 9 1]) datenum([2005 9 21])...
    datenum([2006 5 4]) datenum([2006 6 14]) datenum([2006 7 10]) datenum([2006 8 8]) datenum([2006 9 4]) datenum([2006 10 10])...
    datenum([2007 4 25]) datenum([2007 5 25]) datenum([2007 7 10]) datenum([2007 8 20]) datenum([2007 9 27])...
    datenum([2008 5 8]) datenum([2008 6 19]) datenum([2008 7 24]) datenum([2008 8 28]) datenum([2008 9 29]) datenum([2008 10 26])...
    datenum([2009 5 20]) datenum([2009 7 02]) datenum([2009 8 06]) datenum([2009 8 31]) datenum([2009 10 21])...
    datenum([2010 5 24]) datenum([2010 6 29]) datenum([2010 8 22]) datenum([2010 10 12])...
    datenum([2011 4 20]) datenum([2011 6 6]) datenum([2011 7 12]) datenum([2011 8 24]) datenum([2011 9 28])...
    datenum([2012 6 18]) datenum([2012 7 9]) datenum([2012 8 7]) datenum([2012 8 27]) datenum([2012 10 4])...
    datenum([2013 6 6]) datenum([2013 7 11]) datenum([2013 8 21]) datenum([2013 10 19])...
    datenum([2014 5 17]) datenum([2014 6 19]) datenum([2014 9 3]) datenum([2014 10 19])...
    datenum([2015 4 21]) datenum([2015 6 2]) datenum([2015 8 21]) datenum([2015 9 28])...
    datenum([2016 5 26]) datenum([2016 7 4]) datenum([2016 8 13]) datenum([2016 9 22])...
    datenum([2017 5 18]) datenum([2017 6 20]) datenum([2017 8 02]) datenum([2017 9 21])...
    datenum([2018 4 27]) datenum([2018 6 9]) datenum([2018 9 10]) datenum([2018 10 22])...
    datenum([2019 5 6]) datenum([2019 6 17]) datenum([2019 8 8]) datenum([2019 9 17])...
    datenum([2020 4 22]) datenum([2020 5 25]) datenum([2020 10 8])];

%Mpar_L(1).jDay_cut = [50,100,150];
%2010 - 0622 (grazing)


%%%%%%
%%%%%%%%%%%% PRODUCTIVITY
%%% Maximum Rubisco Capacity
%[Vmax_H]=Maximum_Rubisco_Capacity(Sl_H,PLNR_H,Nl_H);
%[Vmax_L]=Maximum_Rubisco_Capacity(Sl_L,PLNR_L,Nl_L);
Vmax_H = [0]; %55
%%%%%%%%%%%%%%%%%%%%%%
L_day=zeros(NNd,1);
for j=2:24:NN
    [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day(j)]= SetSunVariables(Datam(j,:),DeltaGMT,Lon,Lat,t_bef,t_aft);
end
Lmax_day = max(L_day);
clear('h_S','delta_S','zeta_S','T_sunrise','T_sunset','L_day')
%%%%%%
%%%%%%%%% DORMANT 1 - MAX GROWTH 2 - NORMAL GROWTH 3 - SENESCENCE 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ci_sunL(1,:) = [Ca(1)]; % %% [umolCO2/mol]
Ci_sunH(1,:) = [Ca(1)]; %% [umolCO2/mol]
Ci_shdL(1,:) = [Ca(1)]; % %% [umolCO2/mol]
Ci_shdH(1,:) = [Ca(1)]; %% [umolCO2/mol]
%%%%%%%
LAI_H(1,:)=[0.0]; %
B_H(1,:,:)= [0 0 0 0 0 0 0 0]; %%
Rrootl_H(1,:) = [0] ;
PHE_S_H(1,:)=[0];
dflo_H(1,:)=[0];
AgeL_H(1,:)=[0];
e_rel_H(1,:)=[0];
hc_H(1,:) =[0]; %%
SAI_H(1,:) = [0]; %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
LAI_L(1,:)=[0.96]; %
B_L(1,:,:)= [29 0 445 320 0 0 29 0]; %%
Rrootl_L(1,:) = [4250] ;
PHE_S_L(1,:)=[3];
dflo_L(1,:)=[0];
AgeL_L(1,:)=[115];
e_rel_L(1,:)=[1];
hc_L(1,:) =[0.05]; %%
SAI_L(1,:) = [0.001]; %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BLit(1)=0.0;  %% [kg DM /m2 PFT] Litter Biomass
%%%%%%%%%%
Nreserve_H(1,:)= [0];
Preserve_H(1,:)= [0];
Kreserve_H(1,:)= [0];
FNC_H(1,:)=1;
NupI_H(1,:)= [0 0 0];
Nreserve_L(1,:)= [1000];
Preserve_L(1,:)= [100];
Kreserve_L(1,:)= [100];
FNC_L(1,:)=1;
NupI_L(1,:)= [0 0 0];
RexmyI(1,:)= [0 0 0]; 
%%%%%%%%%

if OPT_SoilBiogeochemistry == 1
    Nreserve_L(1,:)= 5.0; %
    Preserve_L(1,:)= 1.0 ; %
    Kreserve_L(1,:)= 6.4 ; %
    FNC_L(1,:)=1;
    NupI_L(1,:)= [ 0.100826 0.01101   0.043436947];
    RexmyI(1,:)= [ 0.033995   0.1972202             0];
     NavlI(1,:)=[  0.244   0.0115  0.0618]; 
    %%%
    %ORIGINAL
%     P(1,:)=   1000*[ 0.124083987036102   0.325573628166407   0.054069975887170                   0                   0   0.101512860630231   0.126564601766440,...
%         0.014062733529605   0.755595303823738   1.096103576915454   6.402023138646436   0.007894696102488   0.013135775003470   0.000162321193831,...
%         0.000104278919687   0.000081160596916   0.000156418379530   0.038738455769914   0.133541727544331   0.065957977030849                   0,...
%         0.002528243856056   0.009999140372626                   0   0.004492687719394   0.913446366042410   0.005555237467077   0.011170245415232,...
%         0.003668130264076                   0   0.000185349798914   0.000017774714467   0.000036321130404   0.000252824385605   0.001098054538178,...
%         0   0.000476491894740   0.137001078180096   0.001107194360735   0.001675533186330   0.000550219539611                   0,...
%         0.000007588306019   0.150000000000000   0.000172964791821   0.015000467306710   0.000000368010905   0.005528056261498                   0,...
%         0.000603437980692   0.016621912220522   0.000015673158772   0.000009163307352   0.000808451498098   0.500387085728242];
%     %%%
    %My test to run from the spin up onwards:
    %could restate here P values from 
Pini = load("C:\Users\jordi\OneDrive - Imperial College London\Desktop\Paper 1\TC_4\RUN_RESULTS\RUN_RESULTSCH_OE2_a.mat",'P');
P(1,:) = Pini.P(15000,:);
%     DepN= 0.00136; %
%     DepP=1.578e-005; %  gP/m2 day
%     DepK=1.27e-004;
% 
% 
%         FertN=0*ones(1,366);
%         FertN([50 100 200]) = 200;
%         FertP=0*ones(1,366);
%         FertP([50 100 200]) = 100;
%         FertK=0*ones(1,366);
%         FerK([50 100 200]) = 100;

%DepN= 0.00136; %
    %DepP=1.578e-005; %  gP/m2 day
    %DepK=1.27e-004;
    
    % OLD SHOULD STILL WORK
        FertN=0*ones(1,366);
        FertN([50 100 200]) = 30;
        FertP=0*ones(1,366);
        FertP([50 100 200]) = 20;
        FertK=0*ones(1,366);
        FerKN([50 100 200]) = 20;
    
%     %% NEW
%     
%     Time = datetime( (datenum([1959 1 1]):datenum([2020 1 1]))', "convertfrom", "datenum");
%     FertN = timetable(Time,...
%         zeros(size(Time)), 'VariableNames', {'Nrate'});
%     
%     FertP = timetable(Time,...
%         zeros(size(Time)), 'VariableNames', {'Prate'});
%     
%     FertK = timetable(Time,...
%         zeros(size(Time)), 'VariableNames', {'Krate'});
%     
%     %add fertilizer quantities and timing here...
% 
% %%%%CAREFUL !!! THE BELOW ONLY WORKS IF WE START POST SPIN UP!!!!!

% FertApplication = readtimetable("C:\Users\jordi\OneDrive - Imperial College London\Desktop\CurrentWork\Paper 1 - GSM\Data\Sites\CH-Oe2(Oensingen crop)\Observations\CH-OE2.Management.allyears.xlsx");
% FertDates = datenum(FertApplication.DateYear);
% Nrate = FertApplication.N;
% Prate = FertApplication.P;
% Krate = FertApplication.K;
% 
% % Ndates = [datenum(2004,3,4)];
% 
% 
% %     Ndates = datenum(1959:1:2015,5,20);
% %     Nrate  = 50*ones(size(Ndates));
% Ndates = [FertDates];
% Nrate = [Nrate];
%     for nn = 1:length(Ndates)
%         ix = find( abs(Ndates(nn) -  datenum(FertN.Time))<0.01);
%         if ~isempty(ix)
%             FertN.Nrate(ix) = Nrate(nn);
%         end
%     end    
%    
% 
% % note - To convert P2O5 to P divide by 2.29.
% %     Pdates = datenum(1959:1:2015,5,20);
% %     Prate  = 50*ones(size(Pdates));
% Pdates = [FertDates];
% Prate = [Prate];
% 
%     for nn = 1:length(Pdates)
%         ix = find( abs(Pdates(nn) -  datenum(FertP.Time))<0.01);
%         if ~isempty(ix)
%             FertP.Prate(ix) = Prate(nn);
%         end
%     end
%     %note -To convert K2O to K divide by 1.21
% %     Kdates = datenum(1959:1:2015,5,20);
% %     Krate  = 50*ones(size(Kdates));
%     Kdates = [FertDates];
%     Krate = [Krate];
%     for nn = 1:length(Kdates)
%         ix = find( abs(Kdates(nn) -  datenum(FertK.Time))<0.01);
%         if ~isempty(ix)
%             FertK.Krate(ix) = Krate(nn);
%         end
%     end

    Upl=0.01; %% Soil production [mm/yr]
    HIST=0;
    %[B_IO]=Biogeochemistry_IO(Zs(ms+1)*1000,Lat,Lon,Upl,HIST,FertN,FertP,FertK,DepN,DepP,DepK);
    [B_IO]=Biogeochemistry_IO(Zs(ms+1),Lat,Lon,Upl);
    B_IO.SC_par=[1 1 1 1];

    B_IO.FertN = FertN;
    B_IO.FertP = FertP;
    B_IO.FertK = FertK;

    PHs=5; 
end
%%%
%%%
TBio_L=1;  %%[ton DM / ha ]
TBio_H=0;  %[ton DM / ha ]
%%%%%%%%%%%%%%%%%
Vx_H=0;  %% [mm/ m2 PFT];
Vl_H=0;  %% [mm/ m2 PFT];
Vx_L=0;   %% [mm/ m2 PFT];
Vl_L=100;   %% [mm/ m2 PFT];
%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Initial Conditions
%%%%%%%%%%%%%%%%%%%% Initial Conditions
SWE(1)=0; %% [mm]
SND(1)=0;
Ts(1)=Ta(1)+2;
Tdamp(1)=15;
Tdp(1,:)= 15*ones(1,ms);
TdpI_H(1,:)=2;
TdpI_L(1,:)=2;
Sdp(1,:)= (0.001*dz).*cv_s.*(Tdp(1,:)); %% [J °C/m^2 K]
%%% Snow_alb = soil_alb initial
snow_alb.dir_vis = 0.2;
snow_alb.dif_vis = 0.2;
snow_alb.dir_nir = 0.2;
snow_alb.dif_nir = 0.2;
In_L(1,:)=0; In_H(1,:)=0;
In_urb(1)=0; In_rock(1)= 0;
In_Litter(1)=0;
SP_wc(1)=0 ; %%[mm]
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%