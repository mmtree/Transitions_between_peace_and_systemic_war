clear all; close all; clc;

set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',11)

%plr = parula;

h = plot(1:10,1:10,1:10,2:11,1:10,2:11,3:12,1:10,1:10,1:10,1:10,2:11,1:10,2:11,3:12,1:10);
c = get(h,'Color');

%% Stable
fig = figure('position', [0, 0, 600, 200]); hold on;

alpha = 0.08;
beta = 0.5;
L=6;
lambda0=1.0;

x = -3:0.1:30;

y1 = alpha*x.^2 - beta*x + beta*lambda0;
y2 = L*tanh((alpha*x.^2)/L) - beta*x + beta*lambda0;  %L=6
y3 = L - beta*x + beta*lambda0;
%y3 = 4*tanh((alpha*x.^2)/4) - beta*x + beta*lambda0;

plot(x, zeros(1,length(x)),'color',[0.5 0.5 0.5]);
plot(zeros(1,length(x)),x,'color',[0.5 0.5 0.5]);

plot(x,y2,'linewidth',1.5,'color',c{6});  %blue
plot(x,y1,':','linewidth',2,'color',c{4});   %purple
plot(x,y3,':','linewidth',2,'color',c{4});   %purple
%plot(x,y3,'k');
ylim([-0.8 1.2]);
xlim([-1 15]);

%{
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/hysteresis_stable','-dpdf','-r400')
%}
%% peace to war
fig = figure('position', [0, 0, 300, 120]); hold on;

alpha = 0.08;
beta = 0.5;
lambda0=1.6;

x = -3:0.1:30;

y1 = alpha*x.^2 - beta*x + beta*lambda0;
y2 = L*tanh((alpha*x.^2)/L) - beta*x + beta*lambda0;
y3 = L - beta*x + beta*lambda0;

plot(x, zeros(1,length(x)),'color',[0.5 0.5 0.5]);
plot(zeros(1,length(x)),x,'color',[0.5 0.5 0.5]);

plot(x,y2,'linewidth',1.5,'color',c{6});  %blue
plot(x,y1,':','linewidth',2,'color',c{4});   %purple
plot(x,y3,':','linewidth',2,'color',c{4});   %purple

%plot(x,y3,'k');
ylim([-0.8 1.2]);
xlim([-1 15]);
%{
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/hysteresis_peace2war','-dpdf','-r400')
%}

%% War to Peace
fig = figure('position', [0, 0, 300, 120]); hold on;

alpha = 0.0515;
beta = 0.5;
lambda0=1.6;

x = -3:0.1:30;

y1 = alpha*x.^2 - beta*x + beta*lambda0;
y2 = L*tanh((alpha*x.^2)/L) - beta*x + beta*lambda0;
y3 = L - beta*x + beta*lambda0;

plot(x, zeros(1,length(x)),'color',[0.5 0.5 0.5]);
plot(zeros(1,length(x)),x,'color',[0.5 0.5 0.5]);

plot(x,y2,'linewidth',1.5,'color',c{6});  %blue
plot(x,y1,':','linewidth',2,'color',c{4});   %purple
plot(x,y3,':','linewidth',2,'color',c{4});   %purple

%plot(x,y3,'k');
ylim([-0.8 1.2]);
xlim([-1 15]);
%{
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/hysteresis_war2peace','-dpdf','-r400')
%}


%%

% norm_diff =norm(abs(V) - abs(V_tan))
% norm_A = norm(A)
% 
% norm_diff/norm_A
