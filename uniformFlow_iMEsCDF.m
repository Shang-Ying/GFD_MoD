clear; close all; clc;

xy_res = 51;
H_res = 51;

% Nodes
[xx,yy] = meshgrid(linspace(0,1,xy_res),linspace(0,1,xy_res));
xy = [xx(:), yy(:)];
inps = xy(~(xy(:,1) == 0 | xy(:,1) == 1 | xy(:,2) == 0 | xy(:,2) == 1),:);
Bdps = xy((xy(:,1) == 0 | xy(:,1) == 1 | xy(:,2) == 0 | xy(:,2) == 1),:);
clear xx yy xy
Bdps(:,3) = atan2(Bdps(:,2)-0.5,Bdps(:,1)-0.5);
Bdps = sortrows(Bdps,3);
Bdps(:,3) = [];

Bdps_R = [Bdps;Bdps(1,:)];
Refine_bd = 0.2; % 0 ~ 1
splineXY = (spline(1:length(Bdps_R), Bdps_R', 1:Refine_bd:length(Bdps_R)))';
splineXY(end,:) = [];
Bdps_R = splineXY;

box_xm = 0;
box_xM = 1;
box_ym = 0;
box_yM = 1;
ctps = [0.24,0.5; 0.5,0.5];
b_box = [box_xm box_xM box_ym box_yM];
ninit = 1001;
dotmax = 5e6;
radius = @(p,ctps) 0.020 + 0.05*(min(pdist2(ctps, p)));
tStart = cputime;
ptol = 0.018;
btol = 0.020;
xy = node_drop_2d_ctps_tol(b_box, Bdps_R, ninit, dotmax, radius, ctps, btol, ptol);

CPU_AFN = cputime - tStart;
[in,on] = inpolygon(xy(:,1),xy(:,2),Bdps(:,1),Bdps(:,2));
xy = xy(in,:);
inps_new = xy;

% Type
Bdps(Bdps(:,2) == 0, 3) = 2;
Bdps(Bdps(:,2) == 1, 3) = 2;
Bdps(Bdps(:,1) == 0, 3) = 1;
Bdps(Bdps(:,1) == 1, 3) = 1;

% BC
Bdps(Bdps(:,2) == 0, 4) = 0;
Bdps(Bdps(:,2) == 1, 4) = 0;
Bdps(Bdps(:,1) == 0, 4) = 0.1;
Bdps(Bdps(:,1) == 1, 4) = 1.1;

% Vector
Bdps(Bdps(:,2) == 0, 5) = 0;
Bdps(Bdps(:,2) == 0, 6) = -1;
Bdps(Bdps(:,2) == 1, 5) = 0;
Bdps(Bdps(:,2) == 1, 6) = 1;
Bdps(Bdps(:,1) == 0, 5) = -1;
Bdps(Bdps(:,1) == 0, 6) = 0; 
Bdps(Bdps(:,1) == 1, 5) = 1;
Bdps(Bdps(:,1) == 1, 6) = 0;
Bdps_tv = Bdps(:,3:4);
Bdps_V = [Bdps(:,1:2) Bdps(:,5:6)];



% Reshape
tp = [Bdps(:,1:2); inps_new];
tp_tvV = [Bdps(:,:); [inps_new zeros(size(inps_new,1),size(Bdps,2)-size(inps_new,2))]];

% GFDM coef
order = 2;
ns = 33;
[Wcoef,idx] = fct_GFDM_Coef(order,tp(:,1),tp(:,2),ns);

% SGS
% pos_known = [10,10];
% val_known = [0];
% pos_sim = tp;
% nsim = 1000;
% V = '1 Exp(0.3*3)';
% options.max = 50;
% options.nsim = nsim;
% [Y_sgs] = gstat_krig(pos_known,val_known,pos_sim,V,options);
% Y_sgs = Y_sgs - mean(Y_sgs);

% % Load SGS
% load('Y_5000.mat');

% MCS (parallel)
% parfor rlz = 1:nsim
%     tp_Y = [tp Y_sgs(:,rlz)];
%     u_MC(:,rlz) = fct_Steady_DeterHead_HLG(Wcoef,idx,ns,tp_Y,Bdps_tv,Bdps_V);
% end
% % % Calculate u_mean and u_var
% [u_var_MC,u_mean_MC] = var(u_MC,0,2);


% Moment Eqs
% [Y_mean,Y_cov,Y_var] = fct_SK(pos_known,val_known,pos_est,sill,range,nugget,m);
Y_mean(1:length(tp),1) = 0;
sill = 1.0;
cor_length = 0.05;
range = cor_length * 3;
nugget = 0;
[Y_cov] = fct_UnCovMat(tp,sill,range,nugget);
tp_cov = [tp(:,1:2) Y_mean Y_cov]; % combine with kriging Y covariance
[G_loop,D_loop,N_loop] = fct_2DStencil_GFDM(tp_tvV);
tStart = cputime;
u_mean_ME = fct_Steady_1stMEHead_HLG(Wcoef,idx,ns,tp_cov,tp_tvV,G_loop,D_loop,N_loop);
[CYY,CYu,Cuu] = fct_Steady_2ndMEHead_HLG(Wcoef,idx,ns,tp_cov,tp_tvV,u_mean_ME,G_loop,D_loop,N_loop);
CPU_MEs = cputime - tStart;
u_var_ME(:,1) = diag(Cuu);


scrsz = get(0,'ScreenSize');
resol = 500;


current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
sz = 25;
% c = linspace(1,10,length(x));
scatter(tp(:,1),tp(:,2),sz,'filled')
axis([0 1 0 1]); axis('square');
box on;
xlabel('{\it x}_1', 'FontSize', 28); ylabel('{\it x}_2', 'FontSize', 28);
ax = gca; ax.FontSize = 28;
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')

[Xm,Ym] = ndgrid(linspace(min(tp(:,1)), max(tp(:,1)), resol), linspace(min(tp(:,2)), max(tp(:,2)), resol));
Zm = griddata(tp(:,1), tp(:,2), u_mean_ME, Xm, Ym);
current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
contourf(Xm, Ym, Zm, 20,'LineColor','none','ShowText','off'); grid off;
axis([0 1 0 1]); axis('square');
clim([0.1 1.1]); colorbar;
xlabel('{\it x}_1', 'FontSize', 28); ylabel('{\it x}_2', 'FontSize', 28);
ax = gca; ax.FontSize = 28;
chb = colorbar; ylabel(chb, 'Mean','FontSize',28, 'FontName', 'Times');
tlm = get(chb,'Limits'); Tspc = linspace(tlm(1),tlm(2),11);
set(chb,'Ticks',Tspc); TL = arrayfun(@(x) sprintf('%.1f',x),Tspc,'un',0);
set(chb,'TickLabels',TL)
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')

Zm = griddata(tp(:,1), tp(:,2), u_var_ME, Xm, Ym);
current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
contourf(Xm, Ym, Zm, 20,'LineColor','none','ShowText','off'); grid off;
axis([0 1 0 1]); axis('square');
clim([0 0.03]); colorbar;
xlabel('{\it x}_1', 'FontSize', 28); ylabel('{\it x}_2', 'FontSize', 28);
ax = gca; ax.FontSize = 28;
chb = colorbar; ylabel(chb, 'Variance','FontSize',28, 'FontName', 'Times');
tlm = get(chb,'Limits'); Tspc = linspace(tlm(1),tlm(2),4);
set(chb,'Ticks',Tspc); TL = arrayfun(@(x) sprintf('%.2f',x),Tspc,'un',0);
set(chb,'TickLabels',TL)
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')

current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
sz = 25;
% c = linspace(1,10,length(x));
scatter3(tp(:,1),tp(:,2),u_mean_ME,sz,'filled')
axis([0 1 0 1 0 1.2]); axis('square');
box on;
xlabel('{\it x}_1', 'FontSize', 28); ylabel('{\it x}_2', 'FontSize', 28); zlabel('{\it H}', 'FontSize', 28);
ax = gca; ax.FontSize = 28;
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')



% MoM CDF
% V(x,y)
K_mean = exp(Y_mean);
V = zeros(length(tp),1);
parfor ith =  1:length(tp)
    V(ith,1) = -0.5.*K_mean(ith).*...
        ((Wcoef(3,1:ns,ith)*u_var_ME((idx(ith,1:ns)),1)) + (Wcoef(5,1:ns,ith)*u_var_ME((idx(ith,1:ns)),1)));
end

g(length(tp),1) = 0;
% beta(x,y)
beta = zeros(length(tp),1);
parfor ith =  1:length(tp)
    beta(ith,1) = K_mean(ith).*...
        ((Wcoef(3,1:ns,ith)*u_mean_ME((idx(ith,1:ns)),1)) + (Wcoef(5,1:ns,ith)*u_mean_ME((idx(ith,1:ns)),1)))+...
        g(ith,1); % no sink/source
end

% alpha
alpha = zeros(length(tp),1);
parfor ith =  1:length(tp)
    alpha(ith,1) = (K_mean(ith).*...
        (((Wcoef(1,1:ns,ith)*u_mean_ME((idx(ith,1:ns)),1))^2 + (Wcoef(2,1:ns,ith)*u_mean_ME((idx(ith,1:ns)),1))^2)) -...
        V(ith,1)) ./ u_var_ME(ith,1); 
end

% % Reduced tpH nodes (mean head as control points)
% b_box_3D = [0,1,0,1,0.1,1.1];
% ninit_3D = [501,501];
% dotmax_3D = 5e5;
% radius_3D = 0.04;
% % ctps = [tp(:,1) tp(:,2) u_mean_ME];
% N = 31;
% z_even = linspace(0.1, 1.1, H_res); % N is the number of points you want along the z-axis
% ctps_3D = [repmat(0.25, size(z_even))' repmat(0.5, size(z_even))' z_even';...
%             repmat(0.5, size(z_even))' repmat(0.5, size(z_even))' z_even'];
% 
% xyz = node_drop_3d_ctps_RH(b_box_3D, ninit_3D, dotmax_3D, radius_3D, ctps_3D);

% x = linspace(0, 1, N);
% y = linspace(0, 1, N);
% z = linspace(0.1, 1.1, N);  % Adjusted z range from 0.1 to 1.1
% [X, Y] = meshgrid(x, y);
% bottomFace = [X(:), Y(:), 0.1*ones(numel(X), 1)];  % z = 0.1
% topFace = [X(:), Y(:), 1.1*ones(numel(X), 1)];  % z = 1.1
% [X, Z] = meshgrid(x, z);
% frontFace = [X(:), zeros(numel(X), 1), Z(:)];  % y = 0
% backFace = [X(:), ones(numel(X), 1), Z(:)];   % y = 1
% [Y, Z] = meshgrid(y, z);
% leftFace = [zeros(numel(Y), 1), Y(:), Z(:)];  % x = 0
% rightFace = [ones(numel(Y), 1), Y(:), Z(:)];   % x = 1
% Bd_tpH = unique([bottomFace; topFace; frontFace; backFace; leftFace; rightFace], 'rows');
% dist2Bd = pdist2(Bd_tpH, xyz);
% [mdis, mindex] = min(dist2Bd);
% delete = find(mdis<(radius/3));
% xyz(delete,:) = [];
% delete = find((xyz(:,3)-0.1)<(radius/3));
% xyz(delete,:) = [];
% delete = find((1.1-xyz(:,3))<(radius/3));
% xyz(delete,:) = [];
% delete = find((xyz(:,2)-0)<(radius/3));
% xyz(delete,:) = [];
% delete = find((1-xyz(:,2))<(radius/3));
% xyz(delete,:) = [];
% delete = find((xyz(:,1)-0)<(radius/3));
% xyz(delete,:) = [];
% delete = find((1-xyz(:,1))<(radius/3));
% xyz(delete,:) = [];
% xyz_Bd = [xyz; Bd_tpH];

[~,HH] = meshgrid(linspace(0,1,size(tp,1)),linspace(0.1,1.1,H_res));
tpH = [repmat(tp,H_res,1) HH(:)];
alpha = [repmat(alpha,H_res,1)];
beta = [repmat(beta,H_res,1)];
g = [repmat(g,H_res,1)];


% current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
% scatter3(xyz(:,1), xyz(:,2), xyz(:,3), sz, xyz(:,3), 'filled');
% axis([0 1 0 1 0 1.2]); axis('square');
% box on;
% xlabel('{\it x}_1', 'FontSize', 28); ylabel('{\it x}_2', 'FontSize', 28); zlabel('{\it H}', 'FontSize', 28);
% ax = gca; ax.FontSize = 28;
% set(findall(gcf,'type','text'), 'FontName', 'Times');
% set(gca, 'FontName', 'Times')

current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
scatter3(tpH(:,1), tpH(:,2), tpH(:,3), sz, tpH(:,3), 'filled');
axis([0 1 0 1 0 1.2]); axis('square');
box on;
xlabel('{\it x}_1', 'FontSize', 28); ylabel('{\it x}_2', 'FontSize', 28); zlabel('{\it H}', 'FontSize', 28);
ax = gca; ax.FontSize = 28;
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')

% tpH nodes
% tpH = xyz_Bd;
alpha_tpH = zeros(size(tpH,1),1);
beta_tpH = zeros(size(tpH,1),1);
g_tpH = zeros(size(tpH,1),1);
Y_mean_CDF = zeros(size(tpH,1),1);
u_mean_CDF = zeros(size(tpH,1),1);
for i = 1:size(tp, 1)
    st = find(ismember(tpH(:,1:2), tp(i, :), 'rows'));
    alpha_tpH(st) = alpha(i);
    beta_tpH(st) = beta(i);
    g_tpH(st) = g(i);
    Y_mean_CDF(st) = Y_mean(i);
    u_mean_CDF(st) = u_mean_ME(i);
end
for i = 1:size(Bdps, 1)
    st = find(ismember(tpH(:,1:2), Bdps(i, 1:2), 'rows'));
    tpH(st,4:5) = Bdps_tv(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% [~,~,HH] = meshgrid(linspace(0,1,xy_res),linspace(0,1,xy_res),linspace(0.1,1.1,H_res));
% tpH = [repmat(tp,H_res,1) HH(:)];
% alpha = [repmat(alpha,H_res,1)];
% beta = [repmat(beta,H_res,1)];
% g = [repmat(g,H_res,1)];
% Define Fh boundary conditions
% % Fh(0.1; x,y) = 0
% % Fh(1.1; x,y) = 1
% % Fh(H; x,y) = F_phi(H; x,y)  (x,y) Gamma_D
% % Grad(Fh(H; x,y).n(x,y)) = 0  (x,y) Gamma_N
% tpH(:,4:5) = repmat([Bdps_tv; zeros(size(inps,1),2)],H_res,1);
% % tpH(tpH(:,1) == 0 & tpH(:,2)>=0,5) = 1; % Value of BCs
% % tpH(tpH(:,1) == 1 & tpH(:,2)>=1,5) = 0; % Value of BCs
tpH((tpH(:,2) == 0),4) = 2; % 2nd BCs
tpH((tpH(:,2) == 0),5) = 0; % Value of BCs
tpH((tpH(:,2) == 1),4) = 2; % 2nd BCs
tpH((tpH(:,2) == 1),5) = 0; % Value of BCs
tpH((tpH(:,1) == 0),4) = 1; % 1nd BCs
tpH((tpH(:,1) == 0),5) = 1; % Value of BCs
tpH((tpH(:,1) == 1),4) = 1; % 1nd BCs
tpH((tpH(:,1) == 1),5) = 0; % Value of BCs
tpH((tpH(:,3) == 0.1),4) = 1; % 1st BCs
tpH((tpH(:,3) == 0.1),5) = 0; % Value of BCs
tpH((tpH(:,3) == 1.1),4) = 1; % 1st BCs
tpH((tpH(:,3) == 1.1),5) = 1; % Value of BCs
tpH(tpH(:,1) == 0,6:8) = repmat([-1 0 0],sum(tpH(:,1) == 0),1); % normal vector
tpH(tpH(:,1) == 1,6:8) = repmat([1 0 0],sum(tpH(:,1) == 1),1); % normal vector
tpH((tpH(:,3) == 0.1),6:8) = repmat([0 0 -1],sum(tpH(:,3) == 0.1),1); % normal vector
tpH((tpH(:,3) == 1.1),6:8) = repmat([0 0 1],sum(tpH(:,3) == 1.1),1); % normal vector
tpH(tpH(:,2) == 0,6:8) = repmat([0 -1 0],sum(tpH(:,2) == 0),1); % normal vector
tpH(tpH(:,2) == 1,6:8) = repmat([0 1 0],sum(tpH(:,2) == 1),1); % normal vector

current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
scatter3(tpH(:,1), tpH(:,2), tpH(:,3), sz, tpH(:,3), 'filled');
axis([0 1 0 1 0 1.2]); axis('square');
box on;
xlabel('{\it x}_1', 'FontSize', 28); ylabel('{\it x}_2', 'FontSize', 28); zlabel('{\it H}', 'FontSize', 28);
ax = gca; ax.FontSize = 28;
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')

% GFDM coef 3D
ns = 53;
[Wcoef_3D,idx] = fct_GFDM_Coef_3D(order,tpH(:,1),tpH(:,2),tpH(:,3),ns);
tStart = cputime;
tic
% GE: K(4+7+9) + [alpha(x,y)(H-h(x,y))+beta(x,y)+g(x,y)](3) = 0
Fh = fct_Steady_CdfHead(Wcoef_3D,idx,ns,tpH,alpha_tpH,beta_tpH,g_tpH,Y_mean_CDF,u_mean_CDF);
CPU_CDF = cputime - tStart;
CPU_CDF_toc = toc;
% % P[h(x) > H = 0.8] = 1 − Fh(H = 0.8; x)
% H = 0.8;
% n_start = int16(1 + ((H-0.1)/((1.1-0.1)/(H_res-1))) * xy_res*xy_res);
% n_end = int16(xy_res*xy_res + ((H-0.1)/((1.1-0.1)/(H_res-1))) * xy_res*xy_res);
% P_H = [tpH(n_start:n_end,1) tpH(n_start:n_end,2) Fh(n_start:n_end)];
% 
% Zm = griddata(P_H(:,1), P_H(:,2), 1-P_H(:,3), Xm, Ym);
% current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
% contourf(Xm, Ym, Zm, 20,'LineColor','none','ShowText','off'); grid off;
% axis([0 1 0 1]); axis('square');
% clim([0 1]); colorbar;
% xlabel('{\it x}_1', 'FontSize', 28); ylabel('{\it x}_2', 'FontSize', 28);
% ax = gca; ax.FontSize = 28;
% chb = colorbar; ylabel(chb, 'P [{\it h} > 0.8] ','FontSize',28, 'FontName', 'Times');
% tlm = get(chb,'Limits'); Tspc = linspace(tlm(1),tlm(2),6);
% set(chb,'Ticks',Tspc); TL = arrayfun(@(x) sprintf('%.1f',x),Tspc,'un',0);
% set(chb,'TickLabels',TL)
% set(findall(gcf,'type','text'), 'FontName', 'Times');
% set(gca, 'FontName', 'Times')
% 
% % P[h(x) > H = 0.6] = 1 − Fh(H = 0.6; x)
% H = 0.6;
% n_start = int16(1 + ((H-0.1)/((1.1-0.1)/(H_res-1))) * xy_res*xy_res);
% n_end = int16(xy_res*xy_res + ((H-0.1)/((1.1-0.1)/(H_res-1))) * xy_res*xy_res);
% P_H = [tpH(n_start:n_end,1) tpH(n_start:n_end,2) Fh(n_start:n_end)];
% 
% Zm = griddata(P_H(:,1), P_H(:,2), 1-P_H(:,3), Xm, Ym);
% current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
% contourf(Xm, Ym, Zm, 20,'LineColor','none','ShowText','off'); grid off;
% axis([0 1 0 1]); axis('square');
% clim([0 1]); colorbar;
% xlabel('{\it x}_1', 'FontSize', 28); ylabel('{\it x}_2', 'FontSize', 28);
% ax = gca; ax.FontSize = 28;
% chb = colorbar; ylabel(chb, 'P [{\it h} > 0.6] ','FontSize',28, 'FontName', 'Times');
% tlm = get(chb,'Limits'); Tspc = linspace(tlm(1),tlm(2),6);
% set(chb,'Ticks',Tspc); TL = arrayfun(@(x) sprintf('%.1f',x),Tspc,'un',0);
% set(chb,'TickLabels',TL)
% set(findall(gcf,'type','text'), 'FontName', 'Times');
% set(gca, 'FontName', 'Times')
% 
% % P[h(x) > H = 0.4] = 1 − Fh(H = 0.4; x)
% H = 0.4;
% n_start = int16(1 + ((H-0.1)/((1.1-0.1)/(H_res-1))) * xy_res*xy_res);
% n_end = int16(xy_res*xy_res + ((H-0.1)/((1.1-0.1)/(H_res-1))) * xy_res*xy_res);
% P_H = [tpH(n_start:n_end,1) tpH(n_start:n_end,2) Fh(n_start:n_end)];
% 
% Zm = griddata(P_H(:,1), P_H(:,2), 1-P_H(:,3), Xm, Ym);
% current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
% contourf(Xm, Ym, Zm, 20,'LineColor','none','ShowText','off'); grid off;
% axis([0 1 0 1]); axis('square');
% clim([0 1]); colorbar;
% xlabel('{\it x}_1', 'FontSize', 28); ylabel('{\it x}_2', 'FontSize', 28);
% ax = gca; ax.FontSize = 28;
% chb = colorbar; ylabel(chb, 'P [{\it h} > 0.4] ','FontSize',28, 'FontName', 'Times');
% tlm = get(chb,'Limits'); Tspc = linspace(tlm(1),tlm(2),6);
% set(chb,'Ticks',Tspc); TL = arrayfun(@(x) sprintf('%.1f',x),Tspc,'un',0);
% set(chb,'TickLabels',TL)
% set(findall(gcf,'type','text'), 'FontName', 'Times');
% set(gca, 'FontName', 'Times')
% 
% figure
H_pick = tpH(:,3);
parfor refer = 1:size(ctps,1)
    pick = find(tpH(:,1) == ctps(refer,1) & tpH(:,2) == ctps(refer,2));
    cdf{refer} = sortrows([H_pick(pick) Fh(pick,1)]);
end

load("FhH_2601.mat")

current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
uq = plot(cdf{1}(:,1), cdf{1}(:,2),":r",'LineWidth',4);
hold on
mcs = plot(linspace(0.1,1.1,H_res), FhH_2601(764,:),"-k",'LineWidth',2);
axis([0.1 1.1 0 1]); axis('square');
lgd = legend([uq mcs],{'iMEs-CDF','MCS'},...
    'Location','southeast');
xlabel('{\it H}', 'FontSize', 28); ylabel('{\it F_h}({\it H}; {\it x}_1= 0.24,{\it x}_2= 0.5)', 'FontSize', 28);
ax = gca; ax.FontSize = 28;
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')

mean(abs(cdf{1}(:,2)-FhH_2601(764,:)'),1)


current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
uq = plot(cdf{2}(:,1), cdf{2}(:,2),":r",'LineWidth',4);
hold on
mcs = plot(linspace(0.1,1.1,H_res), FhH_2601(1401,:),"-k",'LineWidth',2);
axis([0.1 1.1 0 1]); axis('square');
lgd = legend([uq mcs],{'iMEs-CDF','MCS'},...
    'Location','southeast');
xlabel('{\it H}', 'FontSize', 28); ylabel('{\it F_h}({\it H}; {\it x}_1= 0.5,{\it x}_2= 0.5)', 'FontSize', 28);
ax = gca; ax.FontSize = 28;
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')

mean(abs(cdf{2}(:,2)-FhH_2601(1401,:)'),1)