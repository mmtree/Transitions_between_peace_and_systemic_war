%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Make the first, intro figure in the matrix control paper %%%%%%%%
%%%%% show model dynamics, and st.dev %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all; clc;
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',11)


%% increase alpha

sz = 10;  %10
rng(2);
y0 = 0.8*randn(sz);
y0 = triu(y0);
y0 = y0 + y0';


com_struc = [ones(sz/2,1);-1*ones(sz/2,1)];
y0 = y0 + 0.4*com_struc*com_struc';
y0 = y0 - diag(diag(y0));
FPs = y0;
lambda_01 = max(eig(FPs)); 
y01 = reshape(y0,[],1);

tspan = 0:0.01:15;

y_all = [];
% the first time period, not enough alpha to create war
a = 0.05;
[t,y] = ode45(@(t,y) odefcn_control(t,y,a,sz,FPs), tspan, y01);
y_all = [y_all; y];
t_all = tspan;
y_end = reshape(y(end,:),sz,sz);
y_end = (y_end + y_end')/2;
img1 = y_end;

%%
% the second time period, larger alpha = war
y01 = reshape(y_end,[],1);
a = 0.1;
[t,y] = ode45(@(t,y) odefcn_control(t,y,a,sz,FPs), tspan, y01);
y_all = [y_all; y(2:end,:)];
t_all = [t_all t_all(end)+tspan(2:end)];
y_end = reshape(y(end,:),sz,sz);
y_end = (y_end + y_end')/2;
img2 = y_end;

% the third time period, small alpha, still war though
y01 = reshape(y_end,[],1);
a = 0.05;
[t,y] = ode45(@(t,y) odefcn_control(t,y,a,sz,FPs), tspan, y01);
y_all = [y_all; y(2:end,:)];
t_all = [t_all t_all(end)+tspan(2:end)];
y_end = reshape(y(end,:),sz,sz);
y_end = (y_end + y_end')/2;
img3 = y_end;

% the fourth time period, smaller alpha, back to peace
y01 = reshape(y_end,[],1);
a = 0.01;
[t,y] = ode45(@(t,y) odefcn_control(t,y,a,sz,FPs), tspan, y01);
y_all = [y_all; y(2:end,:)];
t_all = [t_all t_all(end)+tspan(2:end)];
y_end = reshape(y(end,:),sz,sz);
y_end = (y_end + y_end')/2;
img4 = y_end;

%% plot our model connectivity timeseries
fig = figure('position', [0, 0, 350, 100]); hold on;
plot(t_all,y_all,'color','k');
%{
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/into_fig_conn_ts','-dpdf','-r400')
%}

%% plot the standard deviation over time
fig = figure('position', [0, 0, 350, 60]); hold on;

st_dev_ts = std(y_all');
plot(t_all,st_dev_ts, 'linewidth',1);
%{
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/into_fig_stdev_ts','-dpdf','-r400')
%}
%% plot image snapshots along the way

% image 1
fig = figure('position', [0, 0, 100, 100]); hold on;
imagesc(fliplr(img1));
pbaspect([1 1 1]);
c=colorbar();
c.FontName = 'Times New Roman';
axis off
%{
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/into_fig_img1','-dpdf','-r400')
%}
% images 2
fig = figure('position', [0, 0, 100, 100]); hold on;
imagesc(fliplr(img2));
pbaspect([1 1 1]);
c=colorbar();
c.FontName = 'Times New Roman';
axis off
%{
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/into_fig_img2','-dpdf','-r400')
%}
% image 3
fig = figure('position', [0, 0, 100, 100]); hold on;
imagesc(fliplr(img3));
pbaspect([1 1 1]);
c=colorbar();
c.FontName = 'Times New Roman';
axis off
%{
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/into_fig_img3','-dpdf','-r400')
%}
% image 4
fig = figure('position', [0, 0, 100, 100]); hold on;
imagesc(fliplr(img4));
pbaspect([1 1 1]);
c=colorbar();
c.FontName = 'Times New Roman';
axis off
%{
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/into_fig_img4','-dpdf','-r400')
%}
%% Marvel et al. dynamics
sz = 10;  %10
rng(2);
y0 = 0.8*randn(sz);
y0 = triu(y0);
y0 = y0 + y0';


com_struc = [ones(sz/2,1);-1*ones(sz/2,1)];
y0 = y0 + 0.4*com_struc*com_struc';
y0 = y0 - diag(diag(y0));
FPs = y0;
lambda_01 = max(eig(FPs)); 
y01 = reshape(y0,[],1);

tspan = 0:0.01:4;


% the first time period, not enough alpha to create war
a = 0.05;
[t,y] = ode45(@(t,y) odefcn_marvel(t,y,a,sz,FPs), tspan, y01);


fig = figure('position', [0, 0, 120, 100]); hold on;
plot(tspan(200:end),y(200:end,1:5:100),'color','k','linewidth',0.5);
ylim([-30,30]);
%{
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/into_fig_conn_ts_marvel','-dpdf','-r400')
%}


%% plot the standard deviation as a function of alpha

h = plot(1:10,1:10,1:10,2:11,1:10,2:11,3:12,1:10,1:10,1:10,1:10,2:11,1:10,2:11,3:12,1:10);
clr = get(h,'Color');

%%%
sz = 10;  %10
rng(2);
y0 = 0.8*randn(sz);
y0 = triu(y0);
y0 = y0 + y0';
com_struc = [ones(sz/2,1);-1*ones(sz/2,1)];
y0 = y0 + 0.4*com_struc*com_struc';
y0 = y0 - diag(diag(y0));
FPs = y0;
lambda_01 = max(eig(FPs)); 
y01 = reshape(y0,[],1);
tspan = 0:0.01:40;
%%%

a_vec = 0:0.001:0.08; %0.1 %0.2   %0.01*ones(1,3); % %


% peace to war
lambda1_vec = [];
std_vec = [];
for a = a_vec
    [t,y] = ode45(@(t,y) odefcn_control(t,y,a,sz,FPs), tspan, y01);
    y_end = reshape(y(end,:),sz,sz);
    y_end = (y_end + y_end')/2;
    lambda1_vec = [lambda1_vec max(eig(y_end))];
    std_vec = [std_vec std(reshape(y_end,[],1))];
    %disp(a);
end

% war to peace
y01 = reshape(y_end,[],1);
lambda1_w2p_vec = [];
std_w2p_vec = [];
a_vec_back = 0.08:-0.001:0;
for a = a_vec_back
    [t,y] = ode45(@(t,y) odefcn_control(t,y,a,sz,FPs), tspan, y01);
    y_end = reshape(y(end,:),sz,sz);
    y_end = (y_end + y_end')/2;
    lambda1_w2p_vec = [lambda1_w2p_vec max(eig(y_end))];
    std_w2p_vec = [std_w2p_vec std(reshape(y_end,[],1))];
    %disp(a);
end

%% analytical transitions in plot 1(c), std for alpha
ap_star = 1/(4*max(eig(FPs)));
aw_star = 1/max(eig(reshape(y01,sz,sz)));

L=8; % check this is true in function
beta = 1; % check this is true in function
lam = max(eig(FPs));
N=sz;
aw_star_better = (L*N*beta^2)/((L*N+lam*beta)^2); %% needs to be fixed

%%
fig = figure('position', [0, 0, 150, 90]); hold on;
plot(a_vec,std_vec,'color','k','linewidth',1);
plot(a_vec_back,std_w2p_vec,'color','k','linewidth',1);

plot(a_vec_back,std_w2p_vec,'.','color',clr{3},'MarkerSize',8,'linewidth',0.5);
plot(a_vec,std_vec,'.','color',clr{4},'MarkerSize',8,'linewidth',0.5);

plot(a_vec(1:18),std_vec(1:18),'.','color',clr{6},'MarkerSize',8,'linewidth',0.5);
plot(a_vec_back(1:27),std_w2p_vec(1:27),'.','color',clr{6},'MarkerSize',8,'linewidth',0.5);

plot(ap_star*[1 1], [0 10],':','color',0.2*[1 1 1],'linewidth',1);
plot(aw_star*[1 1], [0 10],':','color',0.2*[1 1 1],'linewidth',1);

xticks(0:0.02:0.08);

xlim([0 0.08]);
ylim([0 10]);
%{
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/intro_fig_hystersis_std','-dpdf','-r800')
%}


%% functions

function dydt = odefcn_control(t,y,a,sz,FPs)
    y = reshape(y,sz,sz);
    dydt = zeros(sz,sz);
    
    
    struc_bal = 8*tanh(a*(y^2)/8);  % values that are too large in absolute value level off
    
    dydt = stability(y,FPs) + struc_bal;  %don't want a square root --- too much
   
    
    %dydt = dydt - diag(diag(dydt));
    dydt = reshape(dydt,[],1);
end

function stab = stability(x,FPs)
    stab = -1*(x-FPs);
    %stab = -1*(x-FPs).^3;
end

function dydt = odefcn_marvel(t,y,a,sz,FPs)
    y = reshape(y,sz,sz);
    dydt = zeros(sz,sz);
    
    
    struc_bal = a*(y^2);  % values that are too large in absolute value level off
    
    dydt = struc_bal;  %don't want a square root --- too much
   
    
    %dydt = dydt - diag(diag(dydt));
    dydt = reshape(dydt,[],1);
end