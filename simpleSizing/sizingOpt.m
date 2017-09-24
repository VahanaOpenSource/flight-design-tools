clear;clc
increment = 10;
payload = [100:increment:600];
options = optimoptions(@fmincon);%,'ScaleProblem',true);
scaleFactors = [1000 1 .01 10];
x = [540 1 .01 10];
x = x./scaleFactors;

outputlegend = {...
    'Fixed',...
    'Motors',...
    'Rotors',...
    'Battery',...
    'Wires',...
    'Wing',...
    'Other',...
    'Payload'};

for i=1:length(payload)
    m_payload=payload(i);
    
    lb = [0 0 0 4]; ub = [inf inf inf inf];
    lb = lb./scaleFactors; ub = ub./scaleFactors;
    [x,fval,exitflag] = fmincon(@(x)sizing(x(1),x(2),x(3),x(4),m_payload,0),...
                x,[],[],[],[],lb,ub,...
                @(x)sizing(x(1),x(2),x(3),x(4),m_payload,1),options);
    [weights(i,:), metrics(i,:)] = sizing(x(1),x(2),x(3),x(4),m_payload,2);
   
    if exitflag>0
        m(i)  = x(1)*scaleFactors(1);
        r(i)  = x(2)*scaleFactors(2);
        t(i)  = x(3)*scaleFactors(3);
        AR(i) = x(4)*scaleFactors(4);
    else
        m(i) = nan;
        r(i) = nan;
        t(i) = nan;
        AR(i) = nan;
        weights(i,:) = nan*weights(i,:);
        metrics(i,:) = nan*metrics(i,:);
    end
end
%%
figure(1)
bar(payload,weights,'stacked');
legend(outputlegend,'Location','best','Orientation','horizontal');
title('eVTOL Scaling Trends w/ Usable Payload')
ylabel('Vehicle Gross Mass [kg]')
grid on;
%%
figure(2); clf
subplot(2,2,1); hold on
plot(payload,2*r); ylabel('Prop diameter [m]');
yyaxis right
plot(payload,metrics(:,9)'); ylabel('Disk loading [N/m^2]')
xlabel('Payload [kg]')

subplot(2,2,2); hold on
plot(payload,1000*t); ylabel('Spar cap thickness [mm]');
yyaxis right
plot(payload,metrics(:,1)'); ylabel('Wingspan [m]')
xlabel('Payload [kg]')

subplot(2,2,3);
plot(payload,AR); hold on; ylabel('Aspect ratio')
yyaxis right
plot(payload,metrics(:,10)'*3600); ylabel('C-rate [1/h]')
xlabel('Payload [kg]');

subplot(2,2,4)
plot(payload,payload./m); hold on
plot(payload,weights(:,4)'./m);
legend('Payload','Battery');
xlabel('Payload [kg]'); ylabel('Mass fraction')
ylim([0 1])
