function ydot = rabbit_myofilament_barODEfile(t,y,pin)
% This function computes the PKA-dependent phosphorylation profiles for
% LTCC, RyR, PLB, PP1, PLM, IKs, ICFTTR, IKr, ICLCa, and myofilament targets. 

%% Assign passed in params

Ltot = pin(1);
LCCtot = pin(2);
RyRtot = pin(3);
PLBtot = pin(4);
TnItot = pin(5);
IKstot = pin(6);
ICFTRtot = pin(7);
PP1_PLBtot = pin(8);
PLMtot = pin(9);
myotot = pin(10);
IKrtot = pin(11);
IClCatot = pin(12);
%% ----- Signaling model parameters -------

% b-AR/Gs module
Ltot_bar = Ltot;            % Ltot_bar      [uM] ** apply agonist concentration here **
sumb1AR = 0.028;            % sumb1AR       [uM] % not used
Gstot = 3.83;               % Gstot         [uM]
Kl = 0.285;                 % Kl            [uM]
Kr = 0.062;                 % Kr            [uM]
Kc = 33.0;                  % Kc            [uM]
k_barkp = 1.1e-3;           % k_barkp       [1/sec]
k_barkm = 2.2e-3;           % k_barkm       [1/sec]
k_pkap = 3.6e-3;            % k_pkap        [1/sec/uM]
k_pkam = 2.2e-3;            % k_pkam        [1/sec]
k_gact = 16.0;              % k_gact        [1/sec]
k_hyd = 0.8;                % k_hyd         [1/sec]
k_reassoc = 1.21e3;         % k_reassoc     [1/sec/uM]
% cAMP module
ACtot = 0.047;              % ACtot         [uM]
ATP = 5.0e3;                % ATP           [uM]
PDE3tot = 0.036;            % PDE3tot       [uM] % Changed from .060 to .036
PDE4tot = 0.036;            % PDE4tot       [uM]
IBMXtot = 0.0;              % IBMXtot       [uM]
Fsktot = 0.0;               % Fsktot        [uM] (10 uM when used)
k_ac_basal = 0.2;           % k_ac_basal    [1/sec]
k_ac_gsa = 8.5;             % k_ac_gsa      [1/sec]
k_ac_fsk = 7.3;             % k_ac_fsk      [1/sec]
Km_basal = 1.03e3;          % Km_basal      [uM]
Km_gsa = 315.0;             % Km_gsa        [uM]
Km_fsk = 860.0;             % Km_fsk        [uM]
Kgsa = 0.4;                 % Kgsa          [uM]
Kfsk = 44.0;                % Kfsk          [uM]
k_pde3 = 3.5;               % k_pde3        [1/sec]
Km_pde3 = 0.15;             % Km_pde3       [uM]
k_pde4 = 5.0;               % k_pde4        [1/sec]
Km_pde4 = 1.3;              % Km_pde4       [uM]
Ki_ibmx = 30.0;             % Ki_ibmx       [uM]
% PKA module
PKAItot = 0.46;             % PKAItot       [uM]
PKAIItot = 0.084;           % PKAIItot      [uM]
PKItot = 0.18;              % PKItot        [uM]
Ka = 9.14;                  % Ka            [uM]
Kb = 1.64;                  % Kb            [uM]
Kd = 4.375;                 % Kd            [uM]
Ki_pki = 0.2e-3;            % Ki_pki        [uM]
% PLB & PP1 module
epsilon = 10;               % epsilon       [none]
PLBtot_bar = PLBtot;        % PLBtot_bar    [uM]
PP1_PLBtot_bar = PP1_PLBtot; % PP1_PLBtot_bar [uM]
Inhib1tot = 0.3;            % Inhib1tot     [uM]
k_pka_plb = 54;             % k_pka_plb     [1/sec]
Km_pka_plb = 21;            % Km_pka_plb    [uM]
k_pp1_plb = 8.5;            % k_pp1_plb     [1/sec]
Km_pp1_plb = 7.0;           % Km_pp1_plb    [uM]
k_pka_i1 = 60;              % k_pka_i1      [1/sec]
Km_pka_i1 = 1.0;            % Km_pka_i1     [uM]
Vmax_pp2a_i1 = 14.0;        % Vmax_pp2a_i1  [uM/sec]
Km_pp2a_i1 = 1.0;           % Km_pp2a_i1    [uM]
Ki_inhib1 = 1.0e-3;         % Ki_inhib1     [uM]
% LCC module
LCCtot_bar = LCCtot;        % LCCtot_bar    [uM]
PKAIIlcctot = 0.025;        % PKAIIlcctot   [uM]
PP1lcctot = 0.025;          % PP1lcctot     [uM]
PP2Alcctot = 0.025;         % PP2Alcctot    [uM]
k_pka_lcc = 54;             % k_pka_lcc     [1/sec]
Km_pka_lcc = 21;            % Km_pka_lcc    [uM]
k_pp1_lcc = 8.52;           % k_pp1_lcc     [1/sec]
Km_pp1_lcc = 3;             % Km_pp1_lcc    [uM]
k_pp2a_lcc = 10.1;          % k_pp2a_lcc    [1/sec]
Km_pp2a_lcc = 3;            % Km_pp2a_lcc   [uM]
% RyR module
RyRtot_bar = RyRtot;        % RyRtot_bar    [uM]
PKAIIryrtot = 0.034;        % PKAIIryrtot   [uM]
PP1ryr = 0.034;             % PP1ryr        [uM]
PP2Aryr = 0.034;            % PP2Aryr       [uM]
kcat_pka_ryr = 54;          % kcat_pka_ryr  [1/sec]
Km_pka_ryr = 21;            % Km_pka_ryr    [uM]
kcat_pp1_ryr = 8.52;        % kcat_pp1_ryr  [1/sec]
Km_pp1_ryr = 7;             % Km_pp1_ryr    [uM]
kcat_pp2a_ryr = 10.1;       % kcat_pp2a_ryr [1/sec]
Km_pp2a_ryr = 4.1;          % Km_pp2a_ryr   [uM]
% TnI module
TnItot_bar = TnItot;        % TnItot_bar    [uM]
PP2Atni = 0.67;             % PP2Atni       [uM]
kcat_pka_tni = 54;          % kcat_pka_tni  [1/sec]
Km_pka_tni = 21;            % Km_pka_tni    [uM]
kcat_pp2a_tni = 10.1;       % kcat_pp2a_tni [1/sec]
Km_pp2a_tni = 4.1;          % Km_pp2a_tni   [uM]
% IKs module
IKstot_bar = IKstot;        % IKstot_bar    [uM]
Yotiao_tot = 0.025;         % Yotiao_tot    [uM]
K_yotiao = 0.1e-3;          % K_yotiao      [uM] ** apply G589D mutation here **
PKAII_ikstot = 0.025;       % PKAII_ikstot  [uM]
PP1_ikstot = 0.025;         % PP1_ikstot    [uM]
k_pka_iks = 1.87; %54;      % k_pka_iks     [1/sec] % adjusted as in Xie et al 2013
Km_pka_iks = 21;            % Km_pka_iks    [uM]
k_pp1_iks = 0.19; %8.52;    % k_pp1_iks     [1/sec] % adjusted as in Xie et al 2013
Km_pp1_iks = 7;             % Km_pp1_iks    [uM]
% CFTR module - Added 04/30/10 by Anthony Soltis
ICFTRtot_bar = ICFTRtot;    % ICFTRtot_bar  [uM]
PKAII_CFTRtot = 0.025;      % PKAII_CFTRtot [uM]
PP1_CFTRtot = 0.025;        % PP1_CFTRtot   [uM]
k_pka_CFTR = 54;            % k_pka_CFTR    [1/sec]
Km_pka_CFTR = 8.5;          % Km_pka_CFTR   [uM]
k_pp1_CFTR = 8.52;          % k_pp1_CFTR    [1/sec]
Km_pp1_CFTR = 7;            % Km_pp1_CFTR   [uM]
% PLM module (from PLB)
PLMtot_bar = PLMtot;        % PLMtot_bar    [uM]
kcat_pka_plm = 54;          % kcat_pka_plm  [1/sec]
Km_pka_plm = 21;            % Km_pka_plm    [uM]
kcat_pp2a_plm = 8.5;        % kcat_pp2a_plm [1/sec]
Km_pp2a_plm = 7.0;          % Km_pp2a_plm   [uM]
% Myofilament module (from TnI)
MYOtot_bar = myotot;        % MYOtot_bar    [uM]
PP2Amyo = 0.67;             % PP2Amyo       [uM]
kcat_pka_myo = 54;          % kcat_pka_myo  [1/sec]
Km_pka_myo = 21;            % Km_pka_myo    [uM]
kcat_pp2a_myo = 10.1;       % kcat_pp2a_myo [1/sec]
Km_pp2a_myo = 4.1;          % Km_pp2a_myo   [uM]
% IKr module (from IKs)
IKrtot_bar = IKrtot;        % IKrtot_bar    [uM]
Yotiaor_tot = 0.025;        % Yotiaor_tot   [uM]
K_yotiaor = 0.1e-3;         % K_yotiaor     [uM]
PKAII_ikrtot = 0.025;       % PKAII_ikrtot  [uM]
PP1_ikrtot = 0.025;         % PP1_ikrtot    [uM]
k_pka_ikr = 1.87; %54;      % k_pka_ikr     [1/sec] % adjusted as in Xie et al 2013
Km_pka_ikr = 21;            % Km_pka_ikr    [uM]
k_pp1_ikr = 0.19; %8.52;    % k_pp1_ikr     [1/sec] % adjusted as in Xie et al 2013
Km_pp1_ikr = 7;             % Km_pp1_ikr    [uM]
% IClCa module (from CFTR)
IClCatot_bar = IClCatot;	% IClCatot_bar  [uM]
PKAII_ClCatot = 0.025;      % PKAII_ClCatot [uM]
PP1_ClCatot = 0.025;        % PP1_ClCatot   [uM]
k_pka_ClCa = 54;            % k_pka_ClCa    [1/sec]
Km_pka_ClCa = 8.5;          % Km_pka_ClCa   [uM]
k_pp1_ClCa = 8.52;          % k_pp1_ClCa    [1/sec]
Km_pp1_ClCa = 7;            % Km_pp1_ClCa   [uM]
%% ----------- SIGNALING MODEL -----------

ydot = zeros(size(y));
%% b-AR module
LR = y(1)*y(2)/Kl;
LRG = LR*y(3)/Kr;
RG = y(2)*y(3)/Kc;
BARKDESENS = k_barkp*(LR+LRG);
BARKRESENS = k_barkm*y(5);
PKADESENS = k_pkap*y(17)*y(4);  
PKARESENS = k_pkam*y(6);
GACT = k_gact*(RG+LRG);
HYD = k_hyd*y(7);
REASSOC = k_reassoc*y(8)*y(9);
ydot(1) = Ltot_bar-LR-LRG-y(1);
ydot(2) = y(4)-LR-LRG-RG-y(2);
ydot(3) = Gstot-LRG-RG-y(3);
ydot(4) = (BARKRESENS-BARKDESENS)+(PKARESENS-PKADESENS);
ydot(5) = BARKDESENS-BARKRESENS;
ydot(6) = PKADESENS-PKARESENS;
ydot(7) = GACT-HYD;
ydot(8) = HYD-REASSOC;
ydot(9) = GACT-REASSOC;
% end b-AR module

%% cAMP module
Gsa_gtp_AC = y(10)*y(12)/Kgsa;
Fsk_AC = y(11)*y(12)/Kfsk;
AC_ACT_BASAL = k_ac_basal*y(12)*ATP/(Km_basal+ATP);	    
AC_ACT_GSA = k_ac_gsa*Gsa_gtp_AC*ATP/(Km_gsa+ATP); 
AC_ACT_FSK = k_ac_fsk*Fsk_AC*ATP/(Km_fsk+ATP);	   
% PDE3_ACT = k_pde3*y(13)*y(16)/(Km_pde3+y(16));	
% PDE4_ACT = k_pde4*y(13)*y(16)/(Km_pde4+y(16));	
PDE3_ACT = k_pde3*PDE3tot*y(16)/(Km_pde3*(1+IBMXtot/Ki_ibmx)+y(16));	% new PDE3 term w IBMX
PDE4_ACT = k_pde4*PDE4tot*y(16)/(Km_pde4*(1+IBMXtot/Ki_ibmx)+y(16));	% new PDE4 term w IBMX
PDE_IBMX = y(13)*y(14)/Ki_ibmx; % not used
ydot(10) = y(7)-Gsa_gtp_AC-y(10);
ydot(11) = Fsktot-Fsk_AC-y(11);
ydot(12) = ACtot-Gsa_gtp_AC-y(12);  % note: assumes Fsk = 0.  Change Gsa_gtp_AC to Fsk_AC for Forskolin.
% ydot(13) = PDE4tot-PDE_IBMX-y(13);
% ydot(14) = IBMXtot-PDE_IBMX-y(14);
ydot(13) = 0;
ydot(14) = 0;
ydot(15) = AC_ACT_BASAL+AC_ACT_GSA+AC_ACT_FSK-PDE3_ACT-PDE4_ACT;
% end cAMP module

%% PKA module
PKI = PKItot*Ki_pki/(Ki_pki+y(17)+y(18));
A2RC_I = (y(17)/Kd)*y(17)*(1+PKI/Ki_pki);
A2R_I = y(17)*(1+PKI/Ki_pki);
A2RC_II = (y(18)/Kd)*y(18)*(1+PKI/Ki_pki);
A2R_II = y(18)*(1+PKI/Ki_pki);
ARC_I = (Ka/y(16))*A2RC_I;
ARC_II = (Ka/y(16))*A2RC_II;
ydot(16) = y(15)-(ARC_I+2*A2RC_I+2*A2R_I)-(ARC_II+2*A2RC_II+2*A2R_II)-y(16);
PKAtemp = Ka*Kb/Kd+Ka*y(16)/Kd+y(16)^2/Kd;
ydot(17) = 2*PKAItot*y(16)^2-y(17)*(1+PKI/Ki_pki)*(PKAtemp*y(17)+y(16)^2);
ydot(18) = 2*PKAIItot*y(16)^2-y(18)*(1+PKI/Ki_pki)*(PKAtemp*y(18)+y(16)^2);
% end PKA module

%% PLB & PP1 module
PLBn = PLBtot_bar-y(19);  % Non-phos = tot - phos
PLB_PHOSPH = k_pka_plb*y(17)*PLBn/(Km_pka_plb+PLBn);
PLB_DEPHOSPH = k_pp1_plb*y(22)*y(19)/(Km_pp1_plb+y(19));
ydot(19) = PLB_PHOSPH-PLB_DEPHOSPH;
 
Inhib1n = Inhib1tot-y(20);  % Non-phos = tot - phos
Inhib1p_PP1 = y(21)*y(22)/Ki_inhib1;
Inhib1_PHOSPH = k_pka_i1*y(17)*Inhib1n/(Km_pka_i1+Inhib1n); 
Inhib1_DEPHOSPH = Vmax_pp2a_i1*y(20)/(Km_pp2a_i1+y(20));
ydot(20) = Inhib1_PHOSPH-Inhib1_DEPHOSPH;
ydot(21) = y(20)-Inhib1p_PP1-y(21);
ydot(22) = PP1_PLBtot_bar-Inhib1p_PP1-y(22);
% end PLB & PP1 module

%% LTCC module
PKAClcc = (PKAIIlcctot/PKAIItot)*y(18);

LCCan = LCCtot_bar-y(23);  % Non-phos = tot - phos
LCCa_PHOSPH = epsilon*k_pka_lcc*PKAClcc*LCCan/(Km_pka_lcc + epsilon*LCCan);
LCCa_DEPHOSPH = epsilon*k_pp2a_lcc*PP2Alcctot*y(23)/(Km_pp2a_lcc+epsilon*y(23));
ydot(23) = LCCa_PHOSPH - LCCa_DEPHOSPH;
 
LCCbn = LCCtot_bar-y(24);  % Non-phos = tot - phos
LCCb_PHOSPH = epsilon*k_pka_lcc*PKAClcc*LCCbn/(Km_pka_lcc+epsilon*LCCbn);   
LCCb_DEPHOSPH = epsilon*k_pp1_lcc*PP1lcctot*y(24)/(Km_pp1_lcc+epsilon*y(24));
ydot(24) = LCCb_PHOSPH-LCCb_DEPHOSPH;
% end LCC module

%% RyR module
PKACryr = (PKAIIryrtot/PKAIItot)*y(18);
RyRn = RyRtot_bar-y(25);  % Non-phos = tot - phos
RyRPHOSPH = epsilon*kcat_pka_ryr*PKACryr*RyRn/(Km_pka_ryr+epsilon*RyRn);
RyRDEPHOSPH1 = epsilon*kcat_pp1_ryr*PP1ryr*y(25)/(Km_pp1_ryr+epsilon*y(25));
RyRDEPHOSPH2A = epsilon*kcat_pp2a_ryr*PP2Aryr*y(25)/(Km_pp2a_ryr+epsilon*y(25));
ydot(25) = RyRPHOSPH-RyRDEPHOSPH1-RyRDEPHOSPH2A;
% end RyR module

%% TnI module
TnIn = TnItot_bar-y(26);  % Non-phos = tot - phos
TnIPHOSPH = kcat_pka_tni*y(17)*TnIn/(Km_pka_tni+TnIn);
TnIDEPHOSPH = kcat_pp2a_tni*PP2Atni*y(26)/(Km_pp2a_tni+y(26));
ydot(26) = TnIPHOSPH-TnIDEPHOSPH;
% end TnI module

%% IKs module
IksYot = y(27)*y(28)/K_yotiao;             % [uM]
ydot(27) = IKstot_bar - IksYot - y(27);    % [uM]
ydot(28) = Yotiao_tot - IksYot - y(28);    % [uM]
PKACiks = (IksYot/IKstot_bar)*(PKAII_ikstot/PKAIItot)*y(18);
PP1iks = (IksYot/IKstot_bar)*PP1_ikstot;
Iksn = IKstot_bar-y(29);  % Non-phos = tot - phos
IKS_PHOSPH = epsilon*k_pka_iks*PKACiks*Iksn/(Km_pka_iks+epsilon*Iksn);
IKS_DEPHOSPH = epsilon*k_pp1_iks*PP1iks*y(29)/(Km_pp1_iks+epsilon*y(29));
ydot(29) = IKS_PHOSPH-IKS_DEPHOSPH;
% end IKs module

%% CFTR module (included 04/30/10)
CFTRn = ICFTRtot_bar - y(30);  % Non-phos = tot - phos
PKAC_CFTR = (PKAII_CFTRtot/PKAIItot)*y(18);    % (PKACFTRtot/PKAIItot)*PKAIIact
CFTRphos = epsilon*CFTRn*PKAC_CFTR*k_pka_CFTR/(Km_pka_CFTR+epsilon*CFTRn);
CFTRdephos = PP1_CFTRtot*k_pp1_CFTR*epsilon*y(30)/(Km_pp1_CFTR+epsilon*y(30));
ydot(30) = CFTRphos - CFTRdephos;
% end CFTR module

%% PLM module (from PLB)
PLMn = PLMtot_bar-y(31);  % Non-phos = tot - phos
PLM_PHOSPH = kcat_pka_plm*y(17)*PLMn/(Km_pka_plm+PLMn);
PLM_DEPHOSPH = kcat_pp2a_plm*y(22)*y(31)/(Km_pp2a_plm+y(31));
ydot(31) = PLM_PHOSPH-PLM_DEPHOSPH;
% end PLM module

%% Myofilament module (from TnI)
Myon = MYOtot_bar-y(32);  % Non-phos = tot - phos
MyoPHOSPH = kcat_pka_myo*y(17)*Myon/(Km_pka_myo+Myon);
MyoDEPHOSPH = kcat_pp2a_myo*PP2Amyo*y(32)/(Km_pp2a_myo+y(32));
ydot(32) = MyoPHOSPH-MyoDEPHOSPH;
% end myofilament module

%% IKr module (from IKs)
IkrYot = y(33)*y(34)/K_yotiaor;             % [uM]
ydot(33) = IKrtot_bar - IkrYot - y(33);     % [uM]
ydot(34) = Yotiaor_tot - IkrYot - y(34);    % [uM]
PKACikr = (IkrYot/IKrtot_bar)*(PKAII_ikrtot/PKAIItot)*y(18);
PP1ikr = (IkrYot/IKrtot_bar)*PP1_ikrtot;
Ikrn = IKrtot_bar-y(35);  % Non-phos = tot - phos
IKR_PHOSPH = epsilon*k_pka_ikr*PKACikr*Ikrn/(Km_pka_ikr+epsilon*Ikrn);
IKR_DEPHOSPH = epsilon*k_pp1_ikr*PP1ikr*y(35)/(Km_pp1_ikr+epsilon*y(35));
ydot(35) = IKR_PHOSPH-IKR_DEPHOSPH;
% end IKr module

%% ICl(Ca) module
ClCan = IClCatot_bar - y(36);  % Non-phos = tot - phos
PKAC_ClCa = (PKAII_ClCatot/PKAIItot)*y(18);    % (PKACFTRtot/PKAIItot)*PKAIIact
ClCaphos = epsilon*ClCan*PKAC_ClCa*k_pka_ClCa/(Km_pka_ClCa+epsilon*ClCan);
ClCadephos = PP1_ClCatot*k_pp1_ClCa*epsilon*y(36)/(Km_pp1_ClCa+epsilon*y(36));
ydot(36) = ClCaphos - ClCadephos;
% end ICl(Ca) module

%%
%ydot(19) = 0; % PLB
%ydot(23) = 0; ydot(24) = 0; % LCCa and LCCb
%ydot(25) = 0; % RyR
% %ydot(26) = 0; % TnI % not used
%ydot(29) = 0; % IKs
%ydot(30) = 0; % ICFTR
%ydot(31) = 0; % PLM
%ydot(35) = 0; % IKr
%ydot(36) = 0; % IClCa

%ydot(32) = 0; % Myofilament

%% Gather odes
% Need to convert all ydot terms that are ODEs (not DAEs) to miliseconds
odes = [4,5,6,7,8,9,13,14,15,19,20,23,24,25,26,29,30,31,32,35,36];
ydot(odes) = ydot(odes).*1e-3;