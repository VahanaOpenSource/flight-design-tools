%Script to Demo using the function "runPolars", an XFOIL interface
%ZRH March '16
clc;clear;close all;fclose all;

%% demonstrate basic usage for 1 Mach, varying alpha and re#
airfoil = '4412';             %use a string of 4 or 5 digits for NACA sections
alpha   = -10:.25:20;         %deg 
Re      = [.2 1 3]*1e6;       %Reynolds Number(s)
Mach    = 0.1;                %Mach # 
[basicPolarName, P1] = runPolars(airfoil,'aSweep',alpha,'rSweep',Re,'Mach',0.15);

%% plot a 2x2 grid of polar results
subplot(2,2,1)
plot(P1.CD,P1.CL,'.-')
grid on; xlabel('C_D'); ylabel('C_L');
title(['NACA ' airfoil]);
legend(num2str(Re')); 

subplot(2,2,2)
plot(P1.alpha,P1.CL,'.-')
grid on; xlabel('\alpha [deg]'); ylabel('C_L');

subplot(2,2,3)
plot(P1.alpha,P1.CM,'.-')
grid on; xlabel('\alpha [deg]'); ylabel('C_M');

subplot(2,2,4)
plot(P1.alpha,P1.CL./P1.CD,'.-')
grid on; xlabel('\alpha [deg]'); ylabel('L/D');

%% identify failed convergence points where runPolars used interpolation
figure; hold on;
for ii = 1:length(P1.Re)
    xfi = ~P1.conFlag(:,ii);  %indices of failed xfoil convergence for this alpha sweep
    h1(ii) = plot(P1.CD(:,ii),P1.CL(:,ii),'.-');
    plot(P1.CD(xfi,ii),P1.CL(xfi,ii),'ko');
end
grid on; xlabel('CD'); ylabel('CL');
title({'Unconverged Points (circled) for';['NACA ' airfoil]})
legend(h1,num2str(Re'));

