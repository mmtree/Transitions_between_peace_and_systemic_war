close all; clear all; clc;
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
%%
%{
fig = figure('position', [0, 0, 300, 300]); hold on;
x = -4:0.01:4;
y =  -x.^2.*sin(10*x/(pi));    %-x.^2.*sin(10*x/(pi));
y2 =  x.^2;                     %x.^2;
y3 = -x.^2;
y4 = -2*x.*(x+2).*(x-2);
plot(x,y);
plot(x,y2);
plot(x,y3);
plot(x,y4);
%}

%%
%{
figure();
x = -10:0.01:10;
y =30*tanh(x/30);
plot(x,y);
%}

%%


%%
c = 0; %control signal strength
sz = 10;  %10
rng(2);
y0 = 0.8*randn(sz);
y0 = triu(y0);
y0 = y0 + y0';

com_struc = [ones(sz/2,1);-1*ones(sz/2,1)];
y0 = y0 + 0.8*com_struc*com_struc';
y0 = y0 - diag(diag(y0));
FPs = y0;
lambda_01 = max(eig(FPs)); 
y01 = reshape(y0,[],1);

tspan = [0 400];
a_vec = 0:0.01:0.1; %0.1 %0.2   %0.01*ones(1,3); % %
eigs_mat = [];
i=1; j = 1;

%%%%%%%%% 
% balance variables
%%%%%%%%%%%%
Bsz = [];
Bsz_null = NaN(100,length(a_vec));
alpha_bal = 2;
ii = 1;

v_inner_1_vec = [];
v_inner_2_vec = [];
V_inner = [];

for a = a_vec
    [t,y] = ode45(@(t,y) odefcn_control(t,y,a,sz,i,j,c,FPs), tspan, y01);
    y_end = reshape(y(end,:),sz,sz);
    y_end = (y_end + y_end')/2;
    eigs_mat = [eigs_mat; eig(y_end)'];
    disp(a);
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%% compute balance %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    y_end = y_end - diag(diag(y_end));
    P = y_end>0.0;
    N = y_end<-0.0; 
    lambda_star = max(max(eig(P+N)),max(eig(P-N)));
    lambda_p = max(eig(double(P)));
    Bsz = [Bsz (1/4)*log(det(alpha_bal*lambda_star*eye(sz) - (P-N))/det(alpha_bal*lambda_star*eye(sz)-(P+N)))];
    num_its = 100;
    for jj = 1:num_its
        num_pos = sum(sum(P))/2; %only count half the matrix
        num_neg = sum(sum(N))/2; %only count half the matrix
        idx_pos = find(triu(P)>0);
        idx_neg = find(triu(N)>0);
        idx_edge = [idx_pos; idx_neg];
        idx_edge_scramble = idx_edge(randperm(length(idx_edge)));
        null_mat = zeros(size(P));
        null_mat(idx_edge_scramble(1:num_pos)) = 1;
        null_mat(idx_edge_scramble(num_pos+1:num_pos+num_neg)) = -1;
        null_mat = null_mat + null_mat';
        null_mat_P = null_mat>0;
        null_mat_N = null_mat<0;
        lambda_star_null = max(max(eig(null_mat_P+null_mat_N)),max(eig(null_mat_P-null_mat_N)));
        lambda_p = max(eig(double(null_mat_P)));
        Bsz_null(jj,ii) = (1/4)*log(det(alpha_bal*lambda_star_null*eye(sz) - (null_mat_P-null_mat_N))/det(alpha_bal*lambda_star_null*eye(sz)-(null_mat_P+null_mat_N)));
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%  compute phi   %%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    P = y_end.*(y_end>0);
    N = -1*y_end.*(y_end<0);

    kp = sum(P);
    kn = sum(N);
    mp = sum(sum(P));
    mn = sum(sum(N)); 
    if mp>0
        Nulp = (kp'*kp)/mp;
        ModP = P-Nulp;
    else
        ModP = P;
    end
    if mn>0
        Nuln = (kn'*kn)/mn;
        ModN = N-Nuln;
    else
        ModN = N;
    end
    
    Mod = ModP - ModN;
    [V1,D1] = eig(Mod);
    A = P-N;
    [Va,Da] = eig(A);
    
    V_inner = [V_inner diag(Va'*Mod*Va)];
    v_inner_1 = Va(:,end)'*Mod*Va(:,end);
    v_inner_2 = Va(:,end-1)'*Mod*Va(:,end-1);
    v_inner_1_vec = [v_inner_1_vec v_inner_1];
    v_inner_2_vec = [v_inner_2_vec v_inner_2];
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    ii = ii+1;

end

a_star = 1/(4*max(eig(FPs)));
%%% balance %%%%
Bsz_null_avg = mean(Bsz_null);
eta = Bsz./Bsz_null_avg;
eta_null = Bsz_null./Bsz_null_avg;
phi = eta-1;
phi_null = eta_null-1;

fig = figure('position', [0, 0, 160, 140]); hold on;
p = plot(eta_null(1:50,:)',a_vec,'o','Markersize',4,'color',[0.6 0.6 0.6 0.8]); 
%title(['Strong balance, ', num2str(UnitsPer),' mo. intervals, WITH state']);
plot(eta,a_vec,'.','markersize',20,'color','b'); 
plot([0 2],a_star*[1 1],'r:','linewidth',2);
%ylabel('$\alpha$'); xlabel('$\eta$');
%{
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/bal_alpha_1','-dpdf','-r100')  %bal_alpha_1
%}
%%%%%%%%%%%%%%%%


fig = figure('position', [0, 0, 300, 140]); hold on;
for i = 1:size(eigs_mat,1)
    scatter(eigs_mat(i,:),a_vec(i)*ones(sz,1),'ko');
end

lambda_1 = (1-sqrt(1-4.*a_vec*lambda_01))./(2*a_vec);

plot(lambda_1,a_vec,'linewidth',2);
a_star = 1/(4*max(eig(FPs)));
plot([-5 25],a_star*[1 1],'r:','linewidth',2);
%ylabel('$\alpha$'); xlabel('$\lambda(X)$')
%{
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/eigvals_alpha_1','-dpdf','-r100'); %eigvals_alpha_1
%}

fig = figure('position', [0, 0, 200, 140]); hold on;
for i = 1:size(V_inner,2)
    scatter(V_inner(:,i),a_vec(i)*ones(sz,1),15,[0.5 0.5 0.5],'filled');
    scatter(v_inner_1_vec(i),a_vec(i),'bo','linewidth',0.5);
    scatter(v_inner_2_vec(i),a_vec(i),'ro','linewidth',0.5);
end
plot([-5 25],a_star*[1 1],'r:','linewidth',2);


%%
rng(2);
y0 = 0.8*randn(sz);
y0 = triu(y0);
y0 = y0 + y0';
y0 = y0 - diag(diag(y0));
FPs = y0;
y01 = reshape(y0,[],1);

tspan = [0 200];
i=2; j=3;
    [t,y] = ode45(@(t,y) odefcn_control(t,y,a,sz,i,j,c,FPs), tspan, y01);
        y_end = reshape(y(end,:),sz,sz);
        y_end = (y_end + y_end')/2;
        eigs_mat = [eigs_mat; eig(y_end)'];
        disp(a);
%}
%%  timeseries edges 1
rng(2);
y0 = 0.8*randn(sz);
y0 = triu(y0);
y0 = y0 + y0';
y0 = y0 - diag(diag(y0));
FPs = y0;
y01 = reshape(y0,[],1);

c = 1.8; %control signal strength
a = 0.1;  %0.07
eigs_mat = [];
tspan = [0 200];

%
i=1; j=5;
[t,y] = ode45(@(t,y) odefcn_control(t,y,a,sz,i,j,c,FPs), tspan, y01);

%
y_tmp = reshape(y(240,:),sz,sz);
y_tmp = (y_tmp + y_tmp')/2;
max(eig(y_tmp))

1/(4*max(eig(y_tmp)))


for i = 1:10
    for j = i:10
        [t,y] = ode45(@(t,y) odefcn_control(t,y,a,sz,i,j,c,FPs), tspan, y01);
        y_end = reshape(y(end,:),sz,sz);
        y_end = (y_end + y_end')/2;
        eigs_mat = [eigs_mat; eig(y_end)'];
        disp(a);
    end
end

%
fig = figure('position', [0, 0, 140, 200]); hold on;
plot(y(1:200,:),'color','k');
ylim([-5,5]);
box off
set(gca,'xticklabel',{[]},'yticklabel',{[]})

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/ts_1','-dpdf','-r100')

%% timeseries edge 2
rng(2);
y0 = 0.8*randn(sz);
y0 = triu(y0);
y0 = y0 + y0';
y0 = y0 - diag(diag(y0));
FPs = y0;
y01 = reshape(y0,[],1);

c = 1.8; %control signal strength
a = 0.07;  %0.07
eigs_mat = [];
tspan = [0 200];

%
i=1; j=5;
[t,y] = ode45(@(t,y) odefcn_control(t,y,a,sz,i,j,c,FPs), tspan, y01);

%
y_tmp = reshape(y(240,:),sz,sz);
y_tmp = (y_tmp + y_tmp')/2;
max(eig(y_tmp))

1/(4*max(eig(y_tmp)))



for i = 1:10
    for j = i:10
        [t,y] = ode45(@(t,y) odefcn_control(t,y,a,sz,i,j,c,FPs), tspan, y01);
        y_end = reshape(y(end,:),sz,sz);
        y_end = (y_end + y_end')/2;
        eigs_mat = [eigs_mat; eig(y_end)'];
        disp(a);
    end
end


%
fig = figure('position', [0, 0, 140, 200]); hold on;
plot(y(1:200,:),'color','k');
ylim([-5,5]);
box off
set(gca,'xticklabel',{[]},'yticklabel',{[]})

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/ts_2','-dpdf','-r100')

%% time series edges 3
rng(2);
y0 = 0.8*randn(sz);
y0 = triu(y0);
y0 = y0 + y0';
y0 = y0 - diag(diag(y0));
FPs = y0;
y01 = reshape(y0,[],1);

c = 1.8; %control signal strength
a = 0.0;  %0.07
eigs_mat = [];
tspan = [0 200];

%
i=1; j=5;
[t,y] = ode45(@(t,y) odefcn_control(t,y,a,sz,i,j,c,FPs), tspan, y01);

%
y_tmp = reshape(y(240,:),sz,sz);
y_tmp = (y_tmp + y_tmp')/2;
max(eig(y_tmp))

1/(4*max(eig(y_tmp)))


for i = 1:10
    for j = i:10
        [t,y] = ode45(@(t,y) odefcn_control(t,y,a,sz,i,j,c,FPs), tspan, y01);
        y_end = reshape(y(end,:),sz,sz);
        y_end = (y_end + y_end')/2;
        eigs_mat = [eigs_mat; eig(y_end)'];
        disp(a);
    end
end


%
fig = figure('position', [0, 0, 140, 200]); hold on;
plot(y(:,:),'color','k');
ylim([-5,5]);
box off
set(gca,'xticklabel',{[]},'yticklabel',{[]})

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/ts_3','-dpdf','-r100')

%%  timeseries edges 4
rng(2);
y0 = 0.8*randn(sz);
y0 = triu(y0);
y0 = y0 + y0';
com_struc = [ones(sz/2,1);-1*ones(sz/2,1)];
y0 = y0 + 0.8*com_struc*com_struc';
y0 = y0 - diag(diag(y0));
FPs = y0;
y01 = reshape(y0,[],1);

c = 1.8; %control signal strength
a = 0.1;  %0.07
eigs_mat = [];
tspan = [0 200];

%
i=1; j=5;
[t,y] = ode45(@(t,y) odefcn_control(t,y,a,sz,i,j,c,FPs), tspan, y01);

%
y_tmp = reshape(y(240,:),sz,sz);
y_tmp = (y_tmp + y_tmp')/2;
max(eig(y_tmp))

1/(4*max(eig(y_tmp)))


for i = 1:10
    for j = i:10
        [t,y] = ode45(@(t,y) odefcn_control(t,y,a,sz,i,j,c,FPs), tspan, y01);
        y_end = reshape(y(end,:),sz,sz);
        y_end = (y_end + y_end')/2;
        eigs_mat = [eigs_mat; eig(y_end)'];
        disp(a);
    end
end

%
fig = figure('position', [0, 0, 140, 200]); hold on;
plot(y(1:200,:),'color','k');
ylim([-5.5,5.5]);
box off
set(gca,'xticklabel',{[]},'yticklabel',{[]})

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/ts_4','-dpdf','-r100')

%%  timeseries edges 5
rng(2);
y0 = 0.8*randn(sz);
y0 = triu(y0);
y0 = y0 + y0';
com_struc = [ones(sz/2,1);-1*ones(sz/2,1)];
y0 = y0 + 0.8*com_struc*com_struc';
y0 = y0 - diag(diag(y0));
FPs = y0;
y01 = reshape(y0,[],1);

c = 1.8; %control signal strength
a = 0.03;  %0.07
eigs_mat = [];
tspan = [0 200];

%
i=1; j=5;
[t,y] = ode45(@(t,y) odefcn_control(t,y,a,sz,i,j,c,FPs), tspan, y01);

%
y_tmp = reshape(y(240,:),sz,sz);
y_tmp = (y_tmp + y_tmp')/2;
max(eig(y_tmp))

1/(4*max(eig(y_tmp)))


for i = 1:10
    for j = i:10
        [t,y] = ode45(@(t,y) odefcn_control(t,y,a,sz,i,j,c,FPs), tspan, y01);
        y_end = reshape(y(end,:),sz,sz);
        y_end = (y_end + y_end')/2;
        eigs_mat = [eigs_mat; eig(y_end)'];
        disp(a);
    end
end

%
fig = figure('position', [0, 0, 140, 200]); hold on;
plot(y(1:200,:),'color','k');
ylim([-5.5,5.5]);
box off
set(gca,'xticklabel',{[]},'yticklabel',{[]})

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/ts_5','-dpdf','-r100')

%%  timeseries edges 6
rng(2);
y0 = 0.8*randn(sz);
y0 = triu(y0);
y0 = y0 + y0';
com_struc = [ones(sz/2,1);-1*ones(sz/2,1)];
y0 = y0 + 0.8*com_struc*com_struc';
y0 = y0 - diag(diag(y0));
FPs = y0;
y01 = reshape(y0,[],1);

c = 1.8; %control signal strength
a = 0;  %0.07
eigs_mat = [];
tspan = [0 200];

%
i=1; j=5;
[t,y] = ode45(@(t,y) odefcn_control(t,y,a,sz,i,j,c,FPs), tspan, y01);

%
y_tmp = reshape(y(240,:),sz,sz);
y_tmp = (y_tmp + y_tmp')/2;
max(eig(y_tmp))

1/(4*max(eig(y_tmp)))


for i = 1:10
    for j = i:10
        [t,y] = ode45(@(t,y) odefcn_control(t,y,a,sz,i,j,c,FPs), tspan, y01);
        y_end = reshape(y(end,:),sz,sz);
        y_end = (y_end + y_end')/2;
        eigs_mat = [eigs_mat; eig(y_end)'];
        disp(a);
    end
end

%
fig = figure('position', [0, 0, 140, 200]); hold on;
plot(y(:,:),'color','k');
ylim([-5.5,5.5]);
box off
set(gca,'xticklabel',{[]},'yticklabel',{[]})

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/ts_6','-dpdf','-r100')


%%
fig = figure('position', [0, 0, 400, 500]); hold on;

subplot(3,2,1);
imagesc(y0); title('IC');
pbaspect([1 1 1]);

subplot(3,2,2);
imagesc(y_end()); title('final');
pbaspect([1 1 1]);

subplot(3,2,[3 4]);
plot(y(1:end,:),'color','k'); title('edge timeseries');
box off
subplot(3,2,5);
plot(eig(y0),ones(sz,1),'bo'); title('IC eigs');

subplot(3,2,6);
plot(eig(y_end),ones(sz,1),'bo'); title('final eigs');
%{
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/IC_final','-dpdf','-r100')
%}
%% functions

function dydt = odefcn_control(t,y,a,sz,i,j,c,FPs)
    y = reshape(y,sz,sz);
    dydt = zeros(sz,sz);
    
    
    struc_bal = 2*tanh(a*(y^2)/2);  % values that are too large in absolute value level off
    
    dydt = stability(y,FPs) + heaviside(t-4)*struc_bal;  %don't want a square root --- too much
    
    if i ~=j
        dydt(i,j) = dydt(i,j) - c*(heaviside(t-7) - heaviside(t-100));
        dydt(j,i) = dydt(i,j);
    end
    
    %dydt = dydt - diag(diag(dydt));
    dydt = reshape(dydt,[],1);
end

function stab = stability(x,FPs)
    stab = -1*(x-FPs);
    %stab = -1*(x-FPs).^3;
end

function out = heaviside(in)
if in<=0
    out = 0;
else
    out = 1;
end
end