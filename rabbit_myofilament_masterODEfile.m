function dydt = rabbit_myofilament_masterODEfile(t,y,p)
% This function calls the ode files for EC coupling, CaM reactions, CaMKII
% phosphorylation module, and PKA phosphorylation module.

global JCaCyt JCaSL JCaDyad;

%% Select modules to use

flag_ECC = 1;    % if 0, module clamped
flag_cam = 1;    % if 0, module clamped
flag_CaMKII = 1; % if 0, module clamped
flag_BAR = 1;    % if 0, module clamped

Ca_clamp = 0; % 0 Ca-free, 1 Ca-clamp
% set Ca_clamp = 1 for clamping [Ca] to the initial value.

%% Collect params and ICs for each module

ny_ECC = 83; % 83 state vars in ECC module
ny_cam = 15; % 15 state vars in each CaM module
ny_CaMKII = 6; % 6 state vars in CaMKII module
ny_BAR = 30+6; % 30+6 state vars in BAR module

% Allocate ICs for each moduel
% Ca_j is y(36), Ca_sl is y(37), Ca_cytosol is y(38)
% y(54) -> y(59) are state transitions for myofilament model
% y(60) -> y(65) are state transitions for mode 1 junctional LCCs
% y(66) -> y(71) are state transitoins for mode 2 junctional LCCs
% y(72) -> y(77) are state transitions for mode 1 sarcolemmal LCCs
% y(78) -> y(83) are state transitions for mode 2 sarcolemmal LCCs

y_ecc = y(1:ny_ECC); 
y_camDyad = y(ny_ECC+1:ny_ECC+ny_cam);
y_camSL = y(ny_ECC+ny_cam+1:ny_ECC+2*ny_cam);
y_camCyt = y(ny_ECC+2*ny_cam+1:ny_ECC+3*ny_cam);
y_CaMKII = y(ny_ECC+3*ny_cam+1:ny_ECC+3*ny_cam+ny_CaMKII);
y_BAR = y(ny_ECC+3*ny_cam+ny_CaMKII+1:ny_ECC+3*ny_cam+ny_CaMKII+ny_BAR);

% break up parameters from master function into modules
paramsCell=mat2cell(p,ones(size(p,1),1),ones(size(p,2),1));
[cycleLength,prot_input_par,CaMtotDyad,BtotDyad,CaMKIItotDyad,CaNtotDyad,PP1totDyad,...
    CaMtotSL,BtotSL,CaMKIItotSL,CaNtotSL,PP1totSL,...
    CaMtotCyt,BtotCyt,CaMKIItotCyt,CaNtotCyt,PP1totCyt...
    LCCtotDyad,RyRtot,PP1_dyad,PP2A_dyad,OA,PLBtot,LCCtotSL,PP1_SL...
    Ligtot,LCCtotBA,RyRtotBA,PLBtotBA,TnItotBA,IKstotBA,ICFTRtotBA,PP1_PLBtot...
    PLMtotBA,MyototBA,IKrtotBA,IClCatotBA,CKIIOE]=paramsCell{:};

K = 135; % [mM]
Mg = 1;  % [mM]
%% Distribute parameters by module

% CaM module
CaDyad = y(36)*1e3; % from ECC model, *** Converting from [mM] to [uM] ***
compart_dyad = 2;
% ** NOTE: Btotdyad being sent to the dyad camODEfile is set to zero, but is used below for transfer between SL and dyad
pCaMDyad = [K, Mg, CaMtotDyad, 0, CaMKIItotDyad, CaNtotDyad, PP1totDyad, CaDyad, cycleLength, compart_dyad];
CaSL = y(37)*1e3; % from ECC model, *** Converting from [mM] to [uM] ***
compartSL = 1;
pCaMSL = [K, Mg, CaMtotSL, BtotSL, CaMKIItotSL, CaNtotSL, PP1totSL, CaSL, cycleLength, compartSL];
CaCyt = y(38)*1e3; % from ECC model, *** Converting from [mM] to [uM] ***
compartCyt = 0;
pCaMCyt = [K, Mg, CaMtotCyt, BtotCyt, CaMKIItotCyt, CaNtotCyt, PP1totCyt, CaCyt, cycleLength, compartCyt];

% CaMKII phosphorylation module 
CaMKIIact_Dyad = CaMKIItotDyad.*(y(ny_ECC+8)+y(ny_ECC+9)+y(ny_ECC+10)+y(ny_ECC+11)); % Multiply total by fraction
CaMKIIact_SL = CaMKIItotSL.*(y(ny_ECC+ny_cam+8)+y(ny_ECC+ny_cam+9)+y(ny_ECC+ny_cam+10)+y(ny_ECC+ny_cam+11));
PP1_PLB_avail = y(ny_ECC+3*ny_cam+ny_CaMKII+22)./PP1_PLBtot + .0091;  % Active PP1 near PLB / total PP1 conc + basal value
pCaMKII = [CaMKIIact_Dyad,LCCtotDyad,RyRtot,PP1_dyad,PP2A_dyad,OA,PLBtot,...
           CaMKIIact_SL,LCCtotSL,PP1_SL,...
           PP1_PLB_avail];

LCC_CKdyadp = y(ny_ECC+3*ny_cam+2)./LCCtotDyad;   % fractional CaMKII-dependent LCC dyad phosphorylation
RyR_CKp = y(ny_ECC+3*ny_cam+4)./RyRtot;           % fractional CaMKII-dependent RyR phosphorylation
PLB_CKp = y(ny_ECC+3*ny_cam+5)./PLBtot;           % fractional CaMKII-dependent PLB phosphorylation

% BAR (PKA phosphorylation) module
pBAR = [Ligtot,LCCtotBA,RyRtotBA,PLBtotBA,TnItotBA,IKstotBA,ICFTRtotBA,PP1_PLBtot,PLMtotBA,MyototBA,IKrtotBA,IClCatotBA];
LCCa_PKAp = y(ny_ECC+3*ny_cam+ny_CaMKII+23)./LCCtotBA;
LCCb_PKAp = y(ny_ECC+3*ny_cam+ny_CaMKII+24)./LCCtotBA;
PLB_PKAn = (PLBtotBA - y(ny_ECC+3*ny_cam+ny_CaMKII+19))./PLBtotBA; % non-phosphorylated PLB targets
RyR_PKAp = y(ny_ECC+3*ny_cam+ny_CaMKII+25)./RyRtotBA;
TnI_PKAp = y(ny_ECC+3*ny_cam+ny_CaMKII+26)./TnItotBA;
IKs_PKAp = y(ny_ECC+3*ny_cam+ny_CaMKII+29)./IKstotBA;
ICFTR_PKAp = y(ny_ECC+3*ny_cam+ny_CaMKII+30)./ICFTRtotBA;
PLM_PKAp = y(ny_ECC+3*ny_cam+ny_CaMKII+31)./PLMtotBA;
Myo_PKAp = y(ny_ECC+3*ny_cam+ny_CaMKII+32)./MyototBA;
IKr_PKAp = y(ny_ECC+3*ny_cam+ny_CaMKII+35)./IKrtotBA;
IClCa_PKAp = y(ny_ECC+3*ny_cam+ny_CaMKII+35)./IClCatotBA;

% ECC module
pECC = [cycleLength,LCC_CKdyadp,RyR_CKp,PLB_CKp,...
        LCCa_PKAp,LCCb_PKAp,PLB_PKAn,RyR_PKAp,TnI_PKAp,IKs_PKAp,ICFTR_PKAp,...
        PLM_PKAp,Myo_PKAp,IKr_PKAp,IClCa_PKAp,CKIIOE,prot_input_par];
%% Solve dydt in each module

if flag_ECC==0,
    dydt_ecc = zeros(1,length(y_ecc))';
else
    dydt_ecc = rabbit_myofilament_eccODEfile(t,y_ecc,pECC);
end

if flag_cam==0,
    dydt_camDyad = zeros(1,length(y_camDyad))';
	dydt_camSL = zeros(1,length(y_camSL))';
	dydt_camCyt = zeros(1,length(y_camCyt))';
    JCaDyad = 0; JCaSL = 0; JCaCyt = 0;
else
    dydt_camDyad = rabbit_myofilament_camODEfile(t,y_camDyad,pCaMDyad);
    dydt_camSL = rabbit_myofilament_camODEfile(t,y_camSL,pCaMSL);
    dydt_camCyt = rabbit_myofilament_camODEfile(t,y_camCyt,pCaMCyt);
end

if flag_CaMKII==0,
    dydt_CaMKIIDyad = zeros(1,length(y_CaMKII))';
else
    dydt_CaMKIIDyad = rabbit_myofilament_camkiiODEfile(t,y_CaMKII,pCaMKII);
end

if flag_BAR==0,
    dydt_BAR = zeros(1,length(y_BAR))';
else
    dydt_BAR = rabbit_myofilament_barODEfile(t,y_BAR,pBAR);
end

% incorporate Ca buffering from CaM, convert JCaCyt from uM/msec to mM/msec
if Ca_clamp == 0,
    dydt_ecc(36) = dydt_ecc(36) + 1e-3*JCaDyad;     % 0 with Ca-clamp
    dydt_ecc(37) = dydt_ecc(37) + 1e-3*JCaSL;       % 0 with Ca-clamp
    dydt_ecc(38) = dydt_ecc(38) + 1e-3*JCaCyt;      % 0 with Ca-clamp
end
    
% incorporate CaM diffusion between compartments
Vmyo = 2.1454e-11;      % [L]
Vdyad = 1.7790e-014;    % [L]
VSL = 6.6013e-013;      % [L]
kDyadSL = 3.6363e-16;	% [L/msec]
kSLmyo = 8.587e-15;     % [L/msec]
k0Boff = 0.0014;        % [s^-1] 
k0Bon = k0Boff/0.2;     % [uM^-1 s^-1] kon = koff/Kd
k2Boff = k0Boff/100;    % [s^-1] 
k2Bon = k0Bon;          % [uM^-1 s^-1]
k4Boff = k2Boff;        % [s^-1]
k4Bon = k0Bon;          % [uM^-1 s^-1]
CaMtotDyad = sum(y_camDyad(1:6))+CaMKIItotDyad*sum(y_camDyad(7:10))+sum(y_camDyad(13:15));
Bdyad = BtotDyad - CaMtotDyad; % [uM dyad]
J_cam_dyadSL = 1e-3*(k0Boff*y_camDyad(1) - k0Bon*Bdyad*y_camSL(1)); % [uM/msec dyad]
J_ca2cam_dyadSL = 1e-3*(k2Boff*y_camDyad(2) - k2Bon*Bdyad*y_camSL(2)); % [uM/msec dyad]
J_ca4cam_dyadSL = 1e-3*(k2Boff*y_camDyad(3) - k4Bon*Bdyad*y_camSL(3)); % [uM/msec dyad]
J_cam_SLmyo = kSLmyo*(y_camSL(1)-y_camCyt(1)); % [umol/msec]
J_ca2cam_SLmyo = kSLmyo*(y_camSL(2)-y_camCyt(2)); % [umol/msec]
J_ca4cam_SLmyo = kSLmyo*(y_camSL(3)-y_camCyt(3)); % [umol/msec]
dydt_camDyad(1) = dydt_camDyad(1) - J_cam_dyadSL;
dydt_camDyad(2) = dydt_camDyad(2) - J_ca2cam_dyadSL;
dydt_camDyad(3) = dydt_camDyad(3) - J_ca4cam_dyadSL;
dydt_camSL(1) = dydt_camSL(1) + J_cam_dyadSL*Vdyad/VSL - J_cam_SLmyo/VSL;
dydt_camSL(2) = dydt_camSL(2) + J_ca2cam_dyadSL*Vdyad/VSL - J_ca2cam_SLmyo/VSL;
dydt_camSL(3) = dydt_camSL(3) + J_ca4cam_dyadSL*Vdyad/VSL - J_ca4cam_SLmyo/VSL;
dydt_camCyt(1) = dydt_camCyt(1) + J_cam_SLmyo/Vmyo;
dydt_camCyt(2) = dydt_camCyt(2) + J_ca2cam_SLmyo/Vmyo;
dydt_camCyt(3) = dydt_camCyt(3) + J_ca4cam_SLmyo/Vmyo;
%% Collect all dydt terms

dydt = [dydt_ecc; dydt_camDyad; dydt_camSL; dydt_camCyt; dydt_CaMKIIDyad; dydt_BAR];