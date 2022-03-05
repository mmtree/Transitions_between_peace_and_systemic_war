
% vary the level of postive ties, homog <--> factions and measure the
% leading eigenvalue and phi
close all; clear all; clc;
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');


set(0,'defaultAxesFontSize',10);

clr = get(gca,'colororder');
close all;

%%

sz = 100;  %10
mat = zeros(sz);
d = 0.4;
mat11 = rand(sz/2)<d;
mat11 = triu(mat11) + triu(mat11,1)';
mat22 = rand(sz/2)<d;
mat22 = triu(mat22) + triu(mat22,1)';

fig = figure('position', [0, 0, 150, 150]); hold on; 
fig = figure('position', [0, 0, 150, 150]); hold on;
fig = figure('position', [0, 0, 150, 150]); hold on;
fig = figure('position', [0, 0, 150, 150]); hold on;

pop_vec = 1.0:-0.02:0.0;  % 1:-0.02:0;

pop = pop_vec(1);
mat12_mask = rand(sz/2)<d;
mat12 = 2*(rand(sz/2)<pop)-1;
mat12 = mat12.*mat12_mask;

FPs = [mat11 mat12;
    mat12' mat22];
 

%%% dynamics parameters %%%
beta = 1;
alpha = 0.001; % 0.02
tspan = 0:0.1:20;

u_vec = [];

for pop = pop_vec
    mat12_mask = rand(sz/2)<d;
    mat12 = 2*(rand(sz/2)<pop)-1;
    mat12 = mat12.*mat12_mask;
    
    A = [mat11 mat12;
        mat12' mat22];
    %A = A/sqrt(sz);
    
    if pop == 0.3
        A30=A;
    end
    
    if pop == 0.7
        A70=A;
    end
    
    
    [V,D] = eig(A);
    v1 = V(:,end);
    v2 = V(:,end-1);
    Ap = A>0;
    An = A<0;
    kp = sum(Ap);
    kn = sum(An);
    mp = sum(kp)/2;
    mn = sum(kn)/2;
    Bp = (kp'*kp)/(2*mp);
    Bn = (kn'*kn)/(2*mn);
    Mp = Ap - Bp;
    Mn = An - Bn;
    M = Mp - Mn;
    phi = diag(V'*M*V);
    phi1 = v1'*M*v1;
    phi2 = v2'*M*v2;

    figure(1);
    plot(diag(D)', pop*ones(1,sz),'.','color',[0.5 0.5 0.5]);

    figure(2);
    plot(phi', pop*ones(1,sz),'.','color',[0.5 0.5 0.5]);
    plot(phi1,pop,'.','color','b','markersize',5);
    plot(phi2,pop,'.','color','c','markersize',5);
    
    %%%%%%%%%%%%%%%%%%%%
    %%% run dynamics %%%
    %%%%%%%%%%%%%%%%%%%%
    y01 = reshape(A,[],1);
    [t,y] = ode45(@(t,y) odefcn_control(t,y,alpha,sz,FPs,beta), tspan, y01);
    y_end = reshape(y(end,:),sz,sz);
    A2 = (y_end + y_end')/2;
    
    if pop == 0.3
        Afin30=A2;
    end
    
    if pop == 0.7
        Afin70=A2;
    end
    
    [V,D] = eig(A2);
    v1 = V(:,end);
    v2 = V(:,end-1);
    Ap = A2>0;
    An = A2<0;
    kp = sum(Ap);
    kn = sum(An);
    mp = sum(kp)/2;
    mn = sum(kn)/2;
    Bp = (kp'*kp)/(2*mp);
    Bn = (kn'*kn)/(2*mn);
    Mp = Ap - Bp;
    Mn = An - Bn;
    M = Mp - Mn;
    phi = diag(V'*M*V);
    phi1 = v1'*M*v1;
    phi2 = v2'*M*v2;
    
    u_vec = [u_vec; v1'*v1(1,1)];

    figure(3);
    plot(diag(D)', pop*ones(1,sz),'.','color',[0.5 0.5 0.5]);
    
    figure(4);
    plot(phi', pop*ones(1,sz),'.','color',[0.5 0.5 0.5]);
    plot(phi1,pop,'.','color','b','markersize',5);
    plot(phi2,pop,'.','color','c','markersize',5);
    
end

%%
figure(1); ylim([0 1]);
%{
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/homog_fact_eigs_before','-dpdf','-r800')
%}
figure(2); ylim([0 1]);
figure(3); ylim([0 1]);
%{
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/homog_fact_eigs_after','-dpdf','-r800')
%}
figure(4); ylim([0 1]);

%%
fig = figure('position', [0, 0, 150, 150]); hold on;
imagesc(flip(u_vec));
caxis([-0.01 0.01])
ylab = [0 0.5 1.0];
set(gca, 'YTick', 0:25:50, 'YTickLabel', ylab) %
ylim([0 50]);
%{
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/homog_fact_vectors','-dpdf','-r800')
%}

%% show examples
fig = figure('position', [0, 0, 150, 150]); hold on;
imagesc(flip(A30));
%{
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/homog_fact_A30_before','-dpdf','-r800')
%}
fig = figure('position', [0, 0, 150, 150]); hold on;
imagesc(flip(A70));
%{
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/homog_fact_A70_before','-dpdf','-r800')
%}
fig = figure('position', [0, 0, 150, 150]); hold on;
imagesc(flip(Afin30));
cl = caxis;

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/homog_fact_A30_after','-dpdf','-r800')
%{
fig = figure('position', [0, 0, 150, 150]); hold on;
imagesc(flip(Afin70));
caxis(cl);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/homog_fact_A70_after','-dpdf','-r800')
%}


%% functions
function dydt = odefcn_control(t,y,alpha,sz,FPs,beta,tspan)
    y = reshape(y,sz,sz);
    dydt = zeros(sz,sz);
    struc_bal = 4*tanh(alpha*(y^2)/4);  % values that are too large in absolute value level off
    dydt = stability(y,FPs,beta) + 1*struc_bal;
    dydt = reshape(dydt,[],1);
end

function stab = stability(x,FPs,beta)
    stab = -beta*(x-FPs);
end

