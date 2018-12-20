
clc
clear all
%% 
RV1 = 2000;
RV2 = 10;
RV3 = 15;

fc = 1;
fm = 1;
fn = 1;
fer = 1;

Ct=50;
NADt=250;




%%%



%close all
%%% Y =[C_cyto;C_mam; C_mito;ADPM;ADPC;NADPHM;VM; P; P_mam; h42;nh42; ] 
%Y0 = [0.1; 0.1;0.1;7;1.5;0.1;0;0;0;0;0];
options= odeset('MaxStep',0.01);
Y0 =[0.04561;0.055717;0.00013841;9622;2311.7;54.196;154.9;0;0;0.99863;0.9975];
%[t, Y] = ode15s(@(t,Y)cadynamicwithIP3(t,Y),[0,100],Y0,options);
%[t, Y] = ode15s(@(t,Y)cadynamicwithIP3(t,Y),[0,100],Y0);
[t, Y] = ode15s(@(t,Y)cadynamicwithIP3_dyn(t,Y),[0,1800],Y0);
C = Y(:,1);
Cnd = Y(:,2);
Cm = Y(:,3);

Cer = RV2*fer*(Ct-C/fc-Cm/(fm*RV3)-Cnd/(fn*RV1));
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
legend('[Ca^{2+}]_{cyto}','[Ca^{2+}]_{mi}','Ca^{2+}]_{mito}(\mu M)', '[Ca^{2+}]_{ER}')
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
ylabel('Potential (mV)')
