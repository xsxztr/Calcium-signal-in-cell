
clc
clear all
%% 
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

Ct=2000;
NADt=250;

C0 = 0.1;%0.1;
Cm0 = 1;%0.1;
Cnd0  = 0.1;
ADPM0 = 7400;%6100;
ADPC0 = 2000;%1150;
NADPHM0 = 110;%90;
VM0= 160;%154;
P0 = 0.1; %0.1;
P_MAM0 = 0.2;%0;
h420=0.5;%0.99863;
nh420 = 0.3;%0.9975;
Cer0 = 197.9993; %RV2*fer*(Ct-C0/fc-Cm0/(fm*RV3)-Cnd0/(fn*RV1));

%close all
%%% Y =[C_cyto;C_mam; C_mito;ADPM;ADPC;NADPHM;VM; P; P_mam; h42;nh42; ] 
%Y0 = [0.1; 0.1;0.1;7;1.5;0.1;0;0;0;0;0];
options= odeset('MaxStep',0.01);
%Y0 =[0.04561;0.055717;0.00013841;9622;2311.7;54.196;154.9;0;0;0.99863;0.9975];
%Y0 =[0.2;0.1;0.1;6100;1150;90;154;0.1;0;0.99863;0.9975];
Y0 = [C0;Cnd0;Cm0;ADPM0;ADPC0;NADPHM0;VM0;P0;P_MAM0;h420;nh420;Cer0];
% C0 = Y0(1);
% Cnd0 = Y0(2);
% Cm0 = Y0(3);
 %% 
[t, Y] = ode15s(@(t,Y)cadynamicwithIP3_ncon(t,Y),[0,3600],Y0);
C = Y(:,1);
Cnd = Y(:,2);
Cm = Y(:,3);

Cer = Y(:,12);
t = t/60;
figure(1)
yyaxis left
plot(t,Y(:,1),'b-');
hold on
plot(t,Y(:,2),'k-');
plot(t,Y(:,3),'g-');
 
xlabel('Time(min)')
ylabel('[Ca^{2+}]_{cyto}(\mu M), Ca^{2+}]_{MAM},Ca^{2+}]_{mito}(\mu M)   ')
hold off

yyaxis right
plot(t,Cer,'r-');
ylabel('[Ca^{2+}]_{ER}(\mu M)')
legend('[Ca^{2+}]_{cyto}','[Ca^{2+}]_{MAM}','Ca^{2+}]_{mito}(\mu M)', '[Ca^{2+}]_{ER}')
% figure(3)
% plot(t,Y(:,2),'k')
% xlabel('Time')
% ylabel('[Ca^{2+}]_{mt} (\mu M)')

figure(2)
plot(t,Y(:,8),'k');
xlabel('Time')
ylabel('IP_3 (\mu M)')

figure(3)
plot(t,Y(:,7),'k');
xlabel('Time')
ylabel('Potential (mV)');

Ctall = Cer/(RV2*fer)+C/fc+Cm/(fm*RV3)+Cnd/(fn*RV1);
figure(4)
plot(t,Ctall);
xlabel('Time')
ylabel('Total Ca^{2} (\mu M)');