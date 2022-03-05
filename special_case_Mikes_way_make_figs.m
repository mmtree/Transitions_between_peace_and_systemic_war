



% make figure with intersecting beta and alpha plot, special case fig
close all; clear all; clc;
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');


set(0,'defaultAxesFontSize',10);

clr = get(gca,'colororder');

beta1 = 1;
beta2 = 2;
beta3 = 1/2;

alpha = 0.03;
L = 6;
N = 10;
mu = 0.5;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plot opposing forces %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure('position', [0, 0, 500, 200]); hold on;
x_vec = 0:0.01:14;
dxdt_alpha = L*tanh((alpha*N*x_vec.^2)/L);
dxdt_beta1 = -beta1*(x_vec-mu);
dxdt_beta2 = -beta2*(x_vec-mu);
dxdt_beta3 = -beta3*(x_vec-mu);

approx_1 = L*ones(1,length(x_vec));
approx_2 =alpha*N*x_vec.^2;
plot(x_vec,dxdt_alpha,'color',clr(4,:),'linewidth',1.5);
plot(x_vec(450:end),approx_1(450:end),'color',clr(4,:),'linewidth',1.8,'linestyle','--');
plot(x_vec(1:452),approx_2(1:452),'color',clr(4,:),'linewidth',1.8,'linestyle','--');

plot(x_vec,-1*dxdt_beta1,'color',clr(6,:),'linewidth',2,'linestyle',':');
plot(x_vec,-1*dxdt_beta2,'color',clr(6,:),'linewidth',2,'linestyle',':');
plot(x_vec,-1*dxdt_beta3,'color',clr(6,:),'linewidth',2,'linestyle',':');
ylim([-1,10]);

%{
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/special_case_mikes_way','-dpdf','-r400')
%}
%% compute fixed point locations

xP_beta1 = (beta1-sqrt(beta1)*sqrt(beta1-4*N*alpha*mu))/(2*N*alpha); % matches fig
xP_beta2 = (beta2-sqrt(beta2)*sqrt(beta2-4*N*alpha*mu))/(2*N*alpha); % matches fig
xP_beta3 = (beta3-sqrt(beta3)*sqrt(beta3-4*N*alpha*mu))/(2*N*alpha); % DNE, complex, correct

xU_beta1 = (beta1+sqrt(beta1)*sqrt(beta1-4*N*alpha*mu))/(2*N*alpha); % matches fig
xU_beta2 = (beta2+sqrt(beta2)*sqrt(beta2-4*N*alpha*mu))/(2*N*alpha); % fake real #, shouldn't exist
% ???? xU_beta2 is wrong, should not be real, but correct if using parabola approx, there in another intercept
xU_beta3 = (beta3+sqrt(beta3)*sqrt(beta3-4*N*alpha*mu))/(2*N*alpha); % DNE, complex, correct

xW_beta1 = (L+beta1*mu)/(beta1); % matches fig
xW_beta2 = (L+beta2*mu)/(beta2); % fake real #, shouldn't exist
xW_beta3 = (L+beta3*mu)/(beta3); % matches fig

piecewise_cutoff = sqrt(L/(alpha*N)); % matches fig

%% functions

% function dydt = odefcn(t,y,beta,alpha,N,L,FPs)
%     y = reshape(y,N,N);
%     dydt = zeros(N,N);
%     
%     struc_bal = L*tanh(alpha*(y^2)/L);  % values that are too large in absolute value level off
%     
%     dydt = stability(y,beta,FPs) + heaviside(t-4)*struc_bal;  %don't want a square root --- too much
%     
%     dydt = reshape(dydt,[],1);
% end
% 
% function stab = stability(x,beta,FPs)
%     stab = -beta*(x-FPs);
% end