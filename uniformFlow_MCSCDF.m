clear; close all; clc;

xy_res = 51;
H_res = 51;

% Nodes
[xx,yy] = meshgrid(linspace(0,1,xy_res),linspace(0,1,xy_res));
xy = [xx(:), yy(:)];
inps = xy(~(xy(:,1) == 0 | xy(:,1) == 1 | xy(:,2) == 0 | xy(:,2) == 1),:);
Bdps = xy((xy(:,1) == 0 | xy(:,1) == 1 | xy(:,2) == 0 | xy(:,2) == 1),:);
clear xx yy
Bdps(:,3) = atan2(Bdps(:,2)-0.5,Bdps(:,1)-0.5);
Bdps = sortrows(Bdps,3);
Bdps(:,3) = [];

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
tp = [Bdps(:,1:2); inps];
tp_tvV = [Bdps(:,:); [inps zeros(size(inps,1),size(Bdps,2)-size(inps,2))]];

% % GFDM coef
% order = 2;
% ns = 13;
% [Wcoef,idx] = fct_GFDM_Coef(order,tp(:,1),tp(:,2),ns);

% SGS S-GeMS
setenv('SEMIVAR_DEF','SGeMS')
nx = xy_res; % number of cell along x
ny = xy_res; % number of cell along y
S = sgems_get_par('sgsim');
S.dim.nx = nx;
S.dim.ny = ny;
S.dim.nz = 1;
% grid cell size
S.dim.dx = 1/(nx-1);
S.dim.dy = 1/(ny-1);
S.dim.dz = 1;
% grid origin
S.dim.x0 = 0;
S.dim.y0 = 0;
S.dim.z0 = 0;
S.XML.parameters.Variogram = sgems_variogram_xml('1 Exp(0.9)');
m = 2000;
S.XML.parameters.Nb_Realizations.value = m;
tStart = cputime;
S = sgems_grid(S);
CPU_SGS = cputime - tStart;
Rest = squeeze(S.D);
SGS_xy = reshape(Rest,ny*nx,m);
indexing_tp = zeros(size(tp, 1), 1);
for i = 1:size(tp, 1)
    indexing_tp(i) = find(ismember(xy, tp(i, :), 'rows'));
end
SGS_reordered = SGS_xy(indexing_tp, :);
clear("SGS_xy","S","Rest")

[G_loop,D_loop,Nxm_loop,NxM_loop,Nym_loop,NyM_loop,n_right,n_left,n_upper,n_lower] = fct_2DStencil_FDM(tp_tvV);
% MCS (parallel)
tStart = cputime;
for rlz = 1:m
    tp_Y = [tp SGS_reordered(:,rlz)];
    u_MC(:,rlz) = fct_Steady_DeterHead_HLG_FDM(tp_Y,tp_tvV,G_loop,D_loop,Nxm_loop,NxM_loop,Nym_loop,NyM_loop,n_right,n_left,n_upper,n_lower);
end
CPU_MCS = cputime - tStart;
% % Calculate u_mean and u_var
[u_var_MC,u_mean_MC] = var(u_MC,0,2);

% % Moment Eqs
% % [Y_mean,Y_cov,Y_var] = fct_SK(pos_known,val_known,pos_est,sill,range,nugget,m);

scrsz = get(0,'ScreenSize');
resol = 500;

[Xm,Ym] = ndgrid(linspace(min(tp(:,1)), max(tp(:,1)), resol), linspace(min(tp(:,2)), max(tp(:,2)), resol));
Zm = griddata(tp(:,1), tp(:,2), mean(SGS_reordered,2), Xm, Ym);
current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
contourf(Xm, Ym, Zm, 30,'LineColor','none','ShowText','off'); grid off;
axis([0 1 0 1]); axis('square');
clim([-5 5]); colorbar;
xlabel('{\it x}_1', 'FontSize', 28); ylabel('{\it x}_2', 'FontSize', 28);
ax = gca; ax.FontSize = 28;
colormap(othercolor('YlGnBu9'));
chb = colorbar; ylabel(chb, 'Log hydraulic conductivity','FontSize',28, 'FontName', 'Times');
tlm = get(chb,'Limits'); Tspc = linspace(tlm(1),tlm(2),11);
set(chb,'Ticks',Tspc); TL = arrayfun(@(x) sprintf('%.0f',x),Tspc,'un',0);
set(chb,'TickLabels',TL)
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')

[Xm,Ym] = ndgrid(linspace(min(tp(:,1)), max(tp(:,1)), resol), linspace(min(tp(:,2)), max(tp(:,2)), resol));
Zm = griddata(tp(:,1), tp(:,2), u_mean_MC, Xm, Ym);
current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
contourf(Xm, Ym, Zm, 50,'LineColor','none','ShowText','off'); grid off;
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

Zm = griddata(tp(:,1), tp(:,2), u_var_MC, Xm, Ym);
current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
contourf(Xm, Ym, Zm, 50,'LineColor','none','ShowText','off'); grid off;
axis([0 1 0 1]); axis('square');
clim([0 0.03]);
colorbar;
xlabel('{\it x}_1', 'FontSize', 28); ylabel('{\it x}_2', 'FontSize', 28);
ax = gca; ax.FontSize = 28;
chb = colorbar; ylabel(chb, 'Variance','FontSize',28, 'FontName', 'Times');
tlm = get(chb,'Limits'); Tspc = linspace(tlm(1),tlm(2),4);
set(chb,'Ticks',Tspc); TL = arrayfun(@(x) sprintf('%.2f',x),Tspc,'un',0);
set(chb,'TickLabels',TL)
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')



% MoM CDF
g(size(tp,1),1) = 0;
tp_Ym = [tp mean(SGS_reordered,2)];
[G_loop,D_loop,Nxm_loop,NxM_loop,Nym_loop,NyM_loop,n_right,n_left,n_upper,n_lower] = fct_CLS2DStencil_FDM(tp_tvV);
tStart = cputime;
[alpha,beta] = fct_Steady_Closure_FDM(tp_Ym,tp_tvV,u_mean_MC,u_var_MC,g,G_loop,D_loop,Nxm_loop,NxM_loop,Nym_loop,NyM_loop,n_right,n_left,n_upper,n_lower);
CPU_CLS = cputime - tStart;

% tpH nodes
[~,~,HH] = meshgrid(linspace(0,1,xy_res),linspace(0,1,xy_res),linspace(0.1,1.1,H_res));
tpABgYu = [tp alpha beta g tp_Ym(:,3) u_mean_MC];
tpABgYu_sort = sortrows(sortrows(tpABgYu,2),1);
% tpH = [repmat(tp,H_res,1) HH(:)];
tpH = [repmat(tpABgYu_sort(:,1:2),H_res,1) HH(:)];
alpha_CDF = [repmat(tpABgYu_sort(:,3),H_res,1)];
beta_CDF = [repmat(tpABgYu_sort(:,4),H_res,1)];
g_CDF = [repmat(tpABgYu_sort(:,5),H_res,1)];
Y_mean_CDF = repmat(tpABgYu_sort(:,6),H_res,1);
u_mean_CDF = repmat(tpABgYu_sort(:,7),H_res,1);
% Define Fh boundary conditions
% Fh(0.1; x,y) = 0
% Fh(1.1; x,y) = 1
% Fh(H; x,y) = F_phi(H; x,y)  (x,y) Gamma_D
% Grad(Fh(H; x,y).n(x,y)) = 0  (x,y) Gamma_N
% tpH(:,4:5) = repmat([Bdps_tv; zeros(size(inps,1),2)],H_res,1);
% tpH(tpH(:,1) == 0 & tpH(:,2)>=0,5) = 1; % Value of BCs
% tpH(tpH(:,1) == 1 & tpH(:,2)>=1,5) = 0; % Value of BCs
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
sz = 10;
scatter3(tpH(:,1), tpH(:,2), tpH(:,3), sz, 'filled','MarkerFaceColor','k','MarkerEdgeColor','k',...
    'MarkerFaceAlpha',.25,'MarkerEdgeAlpha',.05);
axis([0 1 0 1 0 1.2]); axis('square');
box on;
xlabel('{\it x}_1', 'FontSize', 28); ylabel('{\it x}_2', 'FontSize', 28); zlabel('{\it H}', 'FontSize', 28);
ax = gca; ax.FontSize = 28;
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')



% CDF Governing Equation: K(4+7+9) + [alpha(x,y)(H-h(x,y))+beta(x,y)+g(x,y)](3) = 0
dH = (1.1-0.1)/(H_res-1);
[Sten.G_loop,Sten.D_loop,Sten.Nxm_loop,Sten.NxM_loop,Sten.Nym_loop,Sten.NyM_loop,Sten.Nzm_loop,Sten.NzM_loop,Sten.n_right,Sten.n_left,Sten.n_front,Sten.n_back,Sten.n_upper,Sten.n_lower] = fct_CDF3DStencil_FDM(tpH,xy_res,H_res);
tStart = cputime;
Fh = fct_Steady_CdfHead_FDM(tpH,alpha_CDF,beta_CDF,g_CDF,Y_mean_CDF,u_mean_CDF,dH,Sten,xy_res);
CPU_CDF = cputime - tStart;
% % figure

current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
sz = 5;
scatter3(tpH(:,1), tpH(:,2), tpH(:,3), sz, Fh, 'filled','MarkerEdgeColor','k',...
    'MarkerFaceAlpha',1.0,'MarkerEdgeAlpha',.05);
axis([0 1 0 1 0 1.2]); axis('square');
box on;
clim([0 1]);
colorbar;
fullColormap = othercolor('Spectral7');
nColors = size(fullColormap, 1);       % Total number of colors in the colormap
upperHalf = fullColormap(ceil(nColors / 2):end, :); % Extract the upper half
colormap(flipud(upperHalf));
chb = colorbar; ylabel(chb, '{\it F}_{\ith}','FontSize',28, 'FontName', 'Times');
tlm = get(chb,'Limits'); Tspc = linspace(0,1,6);
set(chb,'Ticks',Tspc); TL = arrayfun(@(x) sprintf('%.1f',x),Tspc,'un',0);
set(chb,'TickLabels',TL)
xlabel('{\it x}_1', 'FontSize', 28); ylabel('{\it x}_2', 'FontSize', 28); zlabel('{\it H}', 'FontSize', 28);
ax = gca; ax.FontSize = 28;
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')

parfor refer = 1:size(tp,1)
    pick = find(tpH(:,1) == tp(refer,1) & tpH(:,2) == tp(refer,2));
    cdf(refer,:) = Fh(pick);
end

% 
% P[h(x) > H = 0.8] = 1 − Fh(H = 0.8; x)
H = 0.8;
Zm = griddata(tpH(abs(tpH(:,3) - H) < eps,1), tpH(abs(tpH(:,3) - H) < eps,2), 1-Fh(abs(tpH(:,3) - H) < eps,1), Xm, Ym);
current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
contourf(Xm, Ym, Zm, 20,'LineColor','none','ShowText','off'); grid off;
axis([0 1 0 1]); axis('square');
xlabel('{\it x}_1', 'FontSize', 28); ylabel('{\it x}_2', 'FontSize', 28);
ax = gca; ax.FontSize = 28;
clim([0 1]);
colorbar;
fullColormap = othercolor('Spectral7');
nColors = size(fullColormap, 1);       % Total number of colors in the colormap
upperHalf = fullColormap(ceil(nColors / 2):end, :); % Extract the upper half
colormap(flipud(upperHalf));
chb = colorbar; ylabel(chb, 'P [{\it h} > 0.8] ','FontSize',28, 'FontName', 'Times');
tlm = get(chb,'Limits'); Tspc = linspace(0,1,6);
set(chb,'Ticks',Tspc); TL = arrayfun(@(x) sprintf('%.1f',x),Tspc,'un',0);
set(chb,'TickLabels',TL)
tlm = get(chb,'Limits'); Tspc = linspace(tlm(1),tlm(2),6);
set(chb,'Ticks',Tspc); TL = arrayfun(@(x) sprintf('%.1f',x),Tspc,'un',0);
set(chb,'TickLabels',TL)
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')
% 
% % P[h(x) > H = 0.6] = 1 − Fh(H = 0.6; x)
H = 0.6;
Zm = griddata(tpH(abs(tpH(:,3) - H) < eps,1), tpH(abs(tpH(:,3) - H) < eps,2), 1-Fh(abs(tpH(:,3) - H) < eps,1), Xm, Ym);
current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
contourf(Xm, Ym, Zm, 20,'LineColor','none','ShowText','off'); grid off;
axis([0 1 0 1]); axis('square');
xlabel('{\it x}_1', 'FontSize', 28); ylabel('{\it x}_2', 'FontSize', 28);
ax = gca; ax.FontSize = 28;
clim([0 1]);
colorbar;
fullColormap = othercolor('Spectral7');
nColors = size(fullColormap, 1);       % Total number of colors in the colormap
upperHalf = fullColormap(ceil(nColors / 2):end, :); % Extract the upper half
colormap(flipud(upperHalf));
chb = colorbar; ylabel(chb, 'P [{\it h} > 0.6] ','FontSize',28, 'FontName', 'Times');
tlm = get(chb,'Limits'); Tspc = linspace(0,1,6);
set(chb,'Ticks',Tspc); TL = arrayfun(@(x) sprintf('%.1f',x),Tspc,'un',0);
set(chb,'TickLabels',TL)
tlm = get(chb,'Limits'); Tspc = linspace(tlm(1),tlm(2),6);
set(chb,'Ticks',Tspc); TL = arrayfun(@(x) sprintf('%.1f',x),Tspc,'un',0);
set(chb,'TickLabels',TL)
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')
% 
% % P[h(x) > H = 0.4] = 1 − Fh(H = 0.4; x)
H = 0.4;
Zm = griddata(tpH(abs(tpH(:,3) - H) < eps,1), tpH(abs(tpH(:,3) - H) < eps,2), 1-Fh(abs(tpH(:,3) - H) < eps,1), Xm, Ym);
current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
contourf(Xm, Ym, Zm, 20,'LineColor','none','ShowText','off'); grid off;
axis([0 1 0 1]); axis('square');
xlabel('{\it x}_1', 'FontSize', 28); ylabel('{\it x}_2', 'FontSize', 28);
ax = gca; ax.FontSize = 28;
clim([0 1]);
colorbar;
fullColormap = othercolor('Spectral7');
nColors = size(fullColormap, 1);       % Total number of colors in the colormap
upperHalf = fullColormap(ceil(nColors / 2):end, :); % Extract the upper half
colormap(flipud(upperHalf));
chb = colorbar; ylabel(chb, 'P [{\it h} > 0.4] ','FontSize',28, 'FontName', 'Times');
tlm = get(chb,'Limits'); Tspc = linspace(0,1,6);
set(chb,'Ticks',Tspc); TL = arrayfun(@(x) sprintf('%.1f',x),Tspc,'un',0);
set(chb,'TickLabels',TL)
tlm = get(chb,'Limits'); Tspc = linspace(tlm(1),tlm(2),6);
set(chb,'Ticks',Tspc); TL = arrayfun(@(x) sprintf('%.1f',x),Tspc,'un',0);
set(chb,'TickLabels',TL)
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')

load("FhH_2601.mat")
D_map = mean(abs(cdf-FhH_2601),2);
D_ave = mean(D_map);

Zm = griddata(tp(:,1), tp(:,2), D_map, Xm, Ym);
current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
contourf(Xm, Ym, Zm, 50,'LineColor','none','ShowText','off'); grid off;
axis([0 1 0 1]); axis('square');
xlabel('{\it x}_1', 'FontSize', 28); ylabel('{\it x}_2', 'FontSize', 28);
ax = gca; ax.FontSize = 28;
clim([0 0.05]);
colorbar;
fullColormap = othercolor('YlOrBr9');
colormap(fullColormap);
chb = colorbar; ylabel(chb, 'CDF Discrepancy','FontSize',28, 'FontName', 'Times');
set(chb,'Ticks',Tspc); TL = arrayfun(@(x) sprintf('%.2f',x),Tspc,'un',0);
set(chb,'TickLabels',TL)
tlm = get(chb,'Limits'); Tspc = linspace(tlm(1),tlm(2),6);
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')

current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
x1 = 0.24;
x2 = 0.5;
pick = tp(:,1) == x1 & tp(:,2) == x2;
uq = plot(linspace(0.1,1.1,H_res), cdf(pick,:),":r",'LineWidth',4);
hold on
mcs = plot(linspace(0.1,1.1,H_res), FhH_2601(pick,:),"-k",'LineWidth',2);
axis([0.1 1.1 0 1]); axis('square');
lgd = legend([uq mcs],{'MCS-CDF','MCS'},...
    'Location','southeast');
xlabel('{\it H}', 'FontSize', 28); ylabel('{\it F_h}({\it H}; {\it x}_1= 0.24,{\it x}_2= 0.5)', 'FontSize', 28);
ax = gca; ax.FontSize = 28;
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')

current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
x1 = 0.5;
x2 = 0.5;
pick = tp(:,1) == x1 & tp(:,2) == x2;
uq = plot(linspace(0.1,1.1,H_res), cdf(pick,:),":r",'LineWidth',4);
hold on
mcs = plot(linspace(0.1,1.1,H_res), FhH_2601(pick,:),"-k",'LineWidth',2);
axis([0.1 1.1 0 1]); axis('square');
lgd = legend([uq mcs],{'MCS-CDF','MCS'},...
    'Location','southeast');
xlabel('{\it H}', 'FontSize', 28); ylabel('{\it F_h}({\it H}; {\it x}_1= 0.5,{\it x}_2= 0.5)', 'FontSize', 28);
ax = gca; ax.FontSize = 28;
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')
