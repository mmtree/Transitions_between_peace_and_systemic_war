sz = 20;

A11 = sprandn(sz/2,sz/2,0.5)+2;
A11 = (A11+A11')/2;
A22 = sprandn(sz/2,sz/2,0.5)+2;
A22 = (A22 + A22')/2;
A12 = sprandn(sz/2,sz/2,0.5) -0.4;
%A21 = randn(sz/2,sz/2);

A = [A11 A12;
    A12' A22];
A = A/sz;

fig = figure('position', [0, 0, 300, 460]);  %500, 850
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

subplot(4,2,1);
imagesc(A);
pbaspect([1 1 1]);
title('A(0)');

[V,D] = eig(A);
eig_vals = diag(D);

subplot(4,2,2); hold on;
scatter(real(eig_vals), imag(eig_vals),'.');
x = -0.3:0.1:2.3; y = -x.*(x-1).*(x-2);
plot(x,0*x,'color',[0.5 0.5 0.5]);
plot(x,y);
plot(x,y-0.3,':','color',[0.5 0.5 0.5]);
plot(x,y+0.3,':','color',[0.5 0.5 0.5]);
ylim([-0.5 0.5]);
pbaspect([1 1 1]);

title('Eigenvalues');

subplot(4,2,3);
imagesc(real(V(:,end)*V(:,end)'));
pbaspect([1 1 1]);
title('1st mode');
caxis([-0.1 0.1])
%colorbar();

subplot(4,2,4);
imagesc(real(V(:,end-1)*V(:,end-1)'));
pbaspect([1 1 1]);
title('2nd mode');
caxis([-0.1 0.1])
%colorbar();

B = zeros(sz);
tspan = 0:0.1:5;
y0_vec = reshape(A,[],1);
[t,y] = ode45(@(t,y) odefcn(t,y,B,sz), tspan, y0_vec);

Afinal = reshape(y(end,:),sz,sz);
subplot(4,2,5);
imagesc(Afinal);
pbaspect([1 1 1]);
caxis([-0.12 0.12])
%colorbar();
title({'A(t) final state','without control'});

eigs_mat = [];
for i = 1:length(tspan)
    tmp_mat = reshape(y(i,:),sz,sz);
    eigs_mat = [eigs_mat eig(tmp_mat)];
end


subplot(4,2,7);
all_time = tspan;
plot(all_time,eigs_mat','.','color','k','markersize',4);
xlabel('time');
title('No control');
ylabel('eigenvalues');

v1 = [ones(1,sz/2) -1*ones(1,sz/2)]'/sqrt(sz);
v2 = ones(sz,1)/sqrt(sz);

B = zeros(sz) - 0.3*v1*v1' + 0.3*v2*v2';
tspan1 = 0:0.1:2;
y0_vec = reshape(A,[],1);
[~,y1] = ode45(@(t,y) odefcn(t,y,B,sz), tspan1, y0_vec);
B = zeros(sz);
tspan = 0:0.1:5;
y0_vec = y1(end,:)';
[t,y] = ode45(@(t,y) odefcn(t,y,B,sz), tspan, y0_vec);
Afinal = reshape(y(end,:),sz,sz);
subplot(4,2,6);
imagesc(Afinal);
pbaspect([1 1 1]);
caxis([-0.12 0.12])
%colorbar();
title({'A(t) final state','with control'});

eigs_mat = [];
for i = 1:length(tspan1)
    tmp_mat = reshape(y1(i,:),sz,sz);
    eigs_mat = [eigs_mat eig(tmp_mat)];
end
for i = 1:length(tspan)
    tmp_mat = reshape(y(i,:),sz,sz);
    eigs_mat = [eigs_mat eig(tmp_mat)];
end

subplot(4,2,8);
all_time = [tspan1 tspan1(end)+tspan];
plot(all_time,eigs_mat','.','color','k','markersize',4);
xlabel('time');
title('Control')

%{
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'figures/community_structure_control','-dpdf','-r800')
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%  functions  %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = odefcn(t,y_vec, B, sz)

A = reshape(y_vec,sz,sz);

dAdt = -A^3 +3*A^2 -2*A +B;

dydt = reshape(dAdt,[],1);
end