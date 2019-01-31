function dYdt=cadynamicwithIP3_buff_ncon(t,Y)
%%% this file is based on the model in the paper by Jungmin han
%%%  and Vipul Periwal " A mathematical model for calcium dynamics..."
%%% Y =[C_cyto;C_mam; C_mito;ADPM;ADPC;NADPHM;VM; P; P_mam; h42;nh42; C_erl;pX;pXn] \
%%% Y =[1    ;2      ; 3    ;4   ; 5  ;6      ;7; 8; 9;     10;  11;   12;   13;14] 

%%
C = Y(1); %% cyto. concentraiton
Cnd = Y(2); %% MAM concentration
Cm = Y(3);  %%% mito. concentration
ADPm = Y(4);
ADPc = Y(5);
NADHm = Y(6);
Vm = Y(7);
P = Y(8);
Pnd = Y(9);
h42 = Y(10);
nh42 = Y(11);
Cer = Y(12);
Xp = Y(13);%%% pump protent on the cytosol side
Xpnd = Y(14);%%% pump protent on the MAM side
%%%
%% Parameters ###

Ps=0.3;

RS1 = 0.25;
RS2 = 0.15;
kipr=0.3;
nkipr=0.3;
Vmcu=0.00005;
nVmcu=0.00005;
Vs=30;
nVs=20;
Vncx=0.5;
nVncx=0.5;

RV1 = 2000;
RV2 = 10;
RV3 = 15;

% fc = 1;
% fm = 1;
% fn = 1;
% fer = 1;

fc = 0.01;
fer = 0.01;
fm =  1;
fn = 1;

Amt =  15000;
Act = 2500;
a1 = 10;
a2 = 3.43;
Cpc = 1.8;


NADt=250;

Dc=0.1;
Dp=1;
tp=1;
Dhyp=1;

kleak=0.001;
nkleak=0.001;

Ks=0.35;
Ke=0.05;

p2=0.016;

K1=19;
K2=0.38;
L=50;
p1=0.1;


kgly=450;
q1=1;
q2=0.1;

ko=600;
q3=100;
q4=177;
q5=5;

Vagc=180;
Kagc=0.14;
p4=0.01;

Vant=5000;
ac=0.111;
am=0.139;
F=96480;
R=8315;
Temp=310.16;

Vf=35000;
q6=10000;
q7=190;
q8=8.5;

Vhyd=150;
Khyd = 1000;

q9=2;
q10=-30;


%%%% parameters for buffer 
pkm2 = 0.97;  %% k-2
pK3s = 1.11111E-5;%% K_3^2
pk4 = 0.4; %% k_4
pPt = 15; % P_t
pk2 = 0.6;%%k_2
pkm4 = 1.2E-3; %%k-4
pK1s =0.7; %% K_1^2
gammap = 10;
Yp = RV2*(pPt-Xp-Xpnd/RV1); %% total protein in ER
fcd = 1/(1+4*C*pK1s/(pK1s+C^2)^2*Xp); %% buffer coefficient of cyto
fcnd = 1/(1+4*Cnd*pK1s/(pK1s+Cnd^2)^2*Xpnd); %% buffer coefficient of MAM
ferd = 1/(1+4*Cer*pK3s/(pK3s*Cer^2+1)^2*Yp); %% buffer coefficient of ER

% figure(4)
% plot(t,fcd,'r*');
% hold on
% plot(t,ferd,'bo');

Jcb = -2*C^2*pK1s*(pk2-pkm4)/(pK1s+C^2)^2*Xp-2*(C^2*pk4-pkm2*pK3s*pK1s*Cer^2)...
       /(1+pK3s*Cer^2)/(pK1s+C^2)*(1-RS1)*(pPt-Xp-Xpnd/RV1); %%SERCA cyto to ER flux
   
Jerb = -2*(pK3s*Cer^2*pkm4*pK1s-pk2*C^2)/(1+pK3s*Cer^2)/(pK1s+C^2)*Xp*RV2...
       -2*(pK3s*Cer^2*pkm4*pK1s-pk2*Cnd^2)/(1+pK3s*Cer^2)/(pK1s+Cnd^2)*Xpnd*RV2/RV1...
       +2*(pK3s*Cer^2*(pk4-pkm2)/(1+pK3s*Cer^2)^2)*(pPt-Xp-Xpnd/RV1)*RV2;%%SERCA ER to Cyto flux
   
   

Jcbnd = -2*Cnd^2*pK1s*(pk2-pkm4)/(pK1s+Cnd^2)^2*Xpnd-2*(Cnd^2*pk4-pkm2*pK3s*pK1s*Cer^2)...
        /(1+pK3s*Cer^2)/(pK1s+Cnd^2)*RS1*RV1*(pPt-Xp-Xpnd/RV1); %%SERCA MAM to ER flux
% Jerbnd = -2*gammap*(pK3s*Cer^2*pkm4*pK1s-pk2*Cnd^2)/(1+pK3s*Cer^2)/(pK1s+Cnd^2)*Xpnd...
%         +2*gammap*(pK3s*Cer^2*(pk4-pkm2)/(1+pK3s*Cer^2)^2)*(pPt-Xpnd);%%SERCA ER to MAM flux

% 
% figure(7)
% plot(t,Jcb,'r*');
% hold on
% plot(t,Jerb,'bo');
% 
% figure(6)
% plot(t,Jcbnd,'r*');
% hold on
% plot(t,Jerbnd,'bo');
%% 
%Cer = RV2*fer*(Ct-C/fc-Cm/(fm*RV3)-Cnd/(fn*RV1));
FRT = F/(R*Temp);


NADm = NADt - NADHm;
ATPm = Amt - ADPm;
ATPc = Act - ADPc; 
%% IPR model ###
 

%Cp0=700;
Cp0=100;
q26=10500;
q62=4010;
phi=q26/(q26+q62);

q42a=1.8*(P)^2/((P)^2+0.34);
V42=110*(P)^2/((P)^2+0.01);

k42=0.49+0.543*(P)^3/((P)^3+64);
km42=0.41+25*(P)^3/((P)^3+274.6);

q24a=1+5/((P)^2+0.25);
V24=62+880/((P)^2+4);

k24=0.35;
km24=80;

%Cp=Cp0*(Cer/680);
Cp=Cp0*(Cer/400);
mi24=Cp^3/(Cp^3+k24^3);
hi24=km24^2/(Cp^2+km24^2);

mi42=C^3/(C^3+k42^3);
hi42=km42^3/(C^3+km42^3);

min42=Cnd^3/(Cnd^3+k42^3);
hin42=km42^3/(Cnd^3+km42^3);

q42=q42a+V42*mi42*h42;
nq42=q42a+V42*min42*nh42;
q24=q24a+V24*(1-mi24*hi24);

D = q42*(q62+q26)/(q42*q62+q42*q26+q24*q62);
nD = nq42*(q62+q26)/(nq42*q62+nq42*q26+q24*q62);
 
L42=0.2;
LH=10;
Lm42=100;
Lh42=(1-D)*L42+D*LH;
nLh42=(1-nD)*L42+nD*LH;

%%
Oipr=phi*D;
nOipr=phi*nD;
 

 
 
 
%% IP_3 dynmaics 
V_delta = 0.05;
K_PLC = 0.1;
k_delta = 1.5;
r_5p = 0.04;
v_3k =  2;
K_D = 0.7;
K_3 = 1;
 
J_PLC = V_delta/(1+Y(8)/k_delta).*Y(1).^2/(Y(1).^2+K_PLC^2);
%J_PLC =0.05.*Y(1).^2/(Y(1).^2+0.3^2); %%% modelf of JTB 251
J_5P = r_5p*Y(8);
J_3K =  v_3k*(Y(1)^4/(Y(1)^4+K_D^4))*Y(8)/(Y(8)+K_3);


%%  Calcium fluxes ###
Kout = 0.5;
Jout = Kout*C;
Jipr	= kipr*Oipr*(Cer-C);
Jleak	= kleak*(Cer-C);
% Jserca	= Vs*(C^2/(Ks^2+C^2))*(ATPc/(Ke+ATPc))*(1-tanh ((Cer-280)/200/0.1))/2; 
Ksercaf =0.26;
nKsercaf =2;
Ksercar =250;
Jserca	= 137*((C/Ksercaf)^0.75-(Cer/Ksercar)^0.75)/(1+(C/Ksercaf)^0.75+(Cer/Ksercar)^0.75)*(ATPc/(Ke+ATPc)); 
% figure(5)
% plot(t,Jserca,'r*');
% hold on

Jncx	= Vncx*(Cm/C)*exp(p2*Vm);
Jmcu	= Vmcu*(C/K1)*(1+C/K1)^3*exp(p1*Vm)/((1+C/K1)^4+(L/(1+C/K2)^2.8));

% mcu_g = 0.0046875; 
% mcu_km = 19;
% mcu_N = 200;
% mcu_p =  0.9;
% Jmcu    = mcu_p*mcu_N*mcu_g*(Y(7)-12.9*log(Y(1)/Y(3)))/(1+mcu_km/Y(1));%% PNAS Mitochondrial calcium uptake

Jdiff	= Dc*(Cnd-C);


nJipr	= nkipr*nOipr*(Cer-Cnd);
nJleak	= nkleak*(Cer-Cnd);
nJserca	= nVs*(Cnd^2/(Ks^2+Cnd^2))*(ATPc/(Ke+ATPc))*(1-tanh ((Cer-280)/200/0.1))/2;
% nJserca	= 137*((Cnd/nKsercaf)^0.75-(Cer/Ksercar)^0.75)/(1+(Cnd/nKsercaf)^0.75+(Cer/Ksercar)^0.75)*(ATPc/(Ke+ATPc)); 
%  figure(6)
% plot(t,nJserca,'r*');
% hold on
nJmcu	= nVmcu*(Cnd/K1)*(1+Cnd/K1)^3*exp(p1*Vm)/((1+Cnd/K1)^4+(L/(1+Cnd/K2)^2.8));
nJncx	= nVncx*(Cm/Cnd)*exp(p2*Vm);
%% ### Metabolic pathways ###
Ipdh	= kgly*(1/(q1+NADHm/NADm))*(Cm/(q2+Cm));
Iagc	= Vagc*(C/(Kagc+C))*(q2/(q2+Cm))*exp(p4*Vm);
Io	= ko*(NADHm/(q3+NADHm))/(1+exp((Vm-q4)/q5));
Jant	= Vant*(1-(ac*ATPc*ADPm/(am*ADPc*ATPm))*exp(-FRT*Vm))/((1+ac*(ATPc/ADPm)*exp(-0.5*FRT*Vm))*(1+ADPm/(am*ATPm)));
If1	= Vf*(q6/(q6+ATPm))/(1+exp((q7-Vm)/q8));
Ihyd	= Jserca*(1-RS1)/2+nJserca*RS1/2+Vhyd*(ATPc/(ATPc+Khyd));

Jhleak	= q9*Vm+q10;
%% Xp dynamics flux
JXp = (pkm2*pK3s*Cer^2+pk4)/(1+pK3s*Cer^2)*(pPt-Xp)-(pk2*C^2+pkm4*pK1s)/(pK1s+C^2)*Xp;
JXpnd = (pkm2*pK3s*Cer^2+pk4)/(1+pK3s*Cer^2)*(pPt-Xpnd)-(pk2*Cnd^2+pkm4*pK1s)/(pK1s+Cnd^2)*Xpnd;
 %% ODEs 
% dCdt = fcd*(Jipr*(1-RS1) - Jserca*(1-RS1) + Jleak*(1-RS1) + Jncx*(1-RS2)/RV3 - Jmcu*(1-RS2)/RV3 + Jdiff-Jout+Jcb);
dCdt = fcd*(Jipr*(1-RS1) +Jcb + Jleak*(1-RS1) + Jncx*(1-RS2)/RV3 - Jmcu*(1-RS2)/RV3 + Jdiff-Jout);

% dCnddt = fn*RV1*(nJipr*RS1 + Jserca*RS1 - Jdiff + nJleak*RS1 + nJncx*RS2/RV3 - nJmcu*RS2/RV3);
dCnddt = fcnd*RV1*(nJipr*RS1 + Jcbnd - Jdiff + nJleak*RS1 + nJncx*RS2/RV3 - nJmcu*RS2/RV3);

dCmdt =  0;% fm*(nJmcu*RS2 - nJncx*RS2 + Jmcu*(1-RS2) - Jncx*(1-RS2));

dADPmdt = 0;%Jant - If1;
dADPcdt = 0;%Ihyd - Jant/RV3;
dNADHmdt = 0;%Ipdh - Io + Iagc;
dVmdt = (a1*Io - a2*If1 - Jant - Jhleak - Jncx*(1-RS2) - nJncx*RS2 - 2*(Jmcu*(1-RS2) + nJmcu*RS2) - Iagc)/Cpc;

dPdt=  J_PLC-J_5P-J_3K;  

dPnddt = RV1*Dp*(P-Pnd);

dh42dt = Lh42*(hi42-h42);
dnh42dt = nLh42*(hin42-nh42);

% dCerdt = ferd*RV2*(-nJipr*RS1 + nJserca*RS1 - nJleak*RS1-Jipr*(1-RS1) + Jserca*(1-RS1)-Jleak*(1-RS1))+ferd*Jerb;
dCerdt = ferd*RV2*(-nJipr*RS1  - nJleak*RS1-Jipr*(1-RS1) -Jleak*(1-RS1))+ferd*(Jerb);
dXpdt = JXp;
dXpnddt = JXpnd;
%% 
dYdt = [dCdt;dCnddt;dCmdt;dADPmdt;dADPcdt;dNADHmdt;dVmdt;dPdt;dPnddt;dh42dt;dnh42dt;dCerdt;dXpdt;dXpnddt];
end