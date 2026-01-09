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

% GFDM coef
order = 2;
ns = 33;
[Wcoef,idx] = fct_GFDM_Coef(order,tp(:,1),tp(:,2),ns);

% Moment Eqs
Y_mean(1:length(tp),1) = 0;
sill = 1.0;
cor_length = 0.05;
range = cor_length * 3;
nugget = 0;
[Y_cov] = fct_UnCovMat(tp,sill,range,nugget);
tp_cov = [tp(:,1:2) Y_mean Y_cov]; % combine with kriging Y covariance
% u_mean_ME = fct_Steady_1stMEHead_HLG_FDM(tp_cov,tp_tvV);
[G_loop,D_loop,N_loop] = fct_2DStencil_GFDM(tp_tvV);
tStart = cputime;
tic
u_mean_ME = fct_Steady_1stMEHead_HLG(Wcoef,idx,ns,tp_cov,tp_tvV,G_loop,D_loop,N_loop);
[CYY,CYu,Cuu] = fct_Steady_2ndMEHead_HLG(Wcoef,idx,ns,tp_cov,tp_tvV,u_mean_ME,G_loop,D_loop,N_loop);
CPU_ME = cputime - tStart;
CPU_ME_toc = toc;
u_var_ME(:,1) = full(diag(Cuu));


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
g(length(tp),1) = 0;
[alpha,beta] = fct_Steady_Closure(Wcoef,idx,ns,Y_mean,u_mean_ME,u_var_ME,g);

% % Reduced tpH nodes (mean head as control points)
% b_box = [0,1,0,1,0.1,1.1];
% ninit = [25,25];
% dotmax = 5e5;
% radius = 0.03;
% ctps = [tp(:,1) tp(:,2) u_mean_ME];
% xyz = node_drop_3d_ctps_RH(b_box, ninit, dotmax, radius, ctps);
% current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
% scatter3(xyz(:,1), xyz(:,2), xyz(:,3), sz, xyz(:,3), 'filled');
% axis([0 1 0 1 0 1.2]); axis('square');
% box on;
% xlabel('{\it x}_1', 'FontSize', 28); ylabel('{\it x}_2', 'FontSize', 28); zlabel('{\it H}', 'FontSize', 28);
% ax = gca; ax.FontSize = 28;
% set(findall(gcf,'type','text'), 'FontName', 'Times');
% set(gca, 'FontName', 'Times')

% tpH nodes
[~,~,HH] = meshgrid(linspace(0,1,xy_res),linspace(0,1,xy_res),linspace(0.1,1.1,H_res));
tpH = [repmat(tp,H_res,1) HH(:)];
alpha = [repmat(alpha,H_res,1)];
beta = [repmat(beta,H_res,1)];
g = [repmat(g,H_res,1)];
% Define Fh boundary conditions
% Fh(0.1; x,y) = 0
% Fh(1.1; x,y) = 1
% Fh(H; x,y) = F_phi(H; x,y)  (x,y) Gamma_D
% Grad(Fh(H; x,y).n(x,y)) = 0  (x,y) Gamma_N
tpH(:,4:5) = repmat([Bdps_tv; zeros(size(inps,1),2)],H_res,1);
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
scatter3(tpH(:,1), tpH(:,2), tpH(:,3), sz, tpH(:,3), 'filled');
axis([0 1 0 1 0 1.2]); axis('square');
box on;
xlabel('{\it x}_1', 'FontSize', 28); ylabel('{\it x}_2', 'FontSize', 28); zlabel('{\it H}', 'FontSize', 28);
ax = gca; ax.FontSize = 28;
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')

Y_mean_CDF = repmat(Y_mean,H_res,1);
u_mean_CDF = repmat(u_mean_ME,H_res,1);

% GFDM coef 3D
ns = 53;
[Wcoef_3D,idx] = fct_GFDM_Coef_3D(order,tpH(:,1),tpH(:,2),tpH(:,3),ns);
% Identify node types


% GE: K(4+7+9) + [alpha(x,y)(H-h(x,y))+beta(x,y)+g(x,y)](3) = 0
tStart = cputime;
tic
profile on
Fh = fct_Steady_CdfHead(Wcoef_3D,idx,ns,tpH,alpha,beta,g,Y_mean_CDF,u_mean_CDF);
profile viewer
CPU_CDF = cputime - tStart;
CPU_CDF_toc = toc;
% figure
parfor refer = 1:size(tp,1)
    pick = find(tpH(:,1) == tp(refer,1) & tpH(:,2) == tp(refer,2));
    cdf(refer,:) = Fh(pick);
end
% plot(cdf')

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
tlm = get(chb,'Limits'); Tspc = linspace(tlm(1),tlm(2),6);
set(chb,'Ticks',Tspc); TL = arrayfun(@(x) sprintf('%.2f',x),Tspc,'un',0);
set(chb,'TickLabels',TL)
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
lgd = legend([uq mcs],{'mMEs-CDF','MCS'},...
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
lgd = legend([uq mcs],{'mMEs-CDF','MCS'},...
    'Location','southeast');
xlabel('{\it H}', 'FontSize', 28); ylabel('{\it F_h}({\it H}; {\it x}_1= 0.5,{\it x}_2= 0.5)', 'FontSize', 28);
ax = gca; ax.FontSize = 28;
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')