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
m = 10000;
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
tic
parfor rlz = 1:m
    tp_Y = [tp SGS_reordered(:,rlz)];
    u_MC(:,rlz) = fct_Steady_DeterHead_HLG_FDM(tp_Y,tp_tvV,G_loop,D_loop,Nxm_loop,NxM_loop,Nym_loop,NyM_loop,n_right,n_left,n_upper,n_lower);
end
CPU_MCS = cputime - tStart;
CPU_MCS_toc = toc;
% % Calculate u_mean and u_var
[u_var_MC,u_mean_MC] = var(u_MC,0,2);
[u_var_MC10,u_mean_MC10] = var(u_MC(:,1:10),0,2);
[u_var_MC20,u_mean_MC20] = var(u_MC(:,1:20),0,2);
[u_var_MC50,u_mean_MC50] = var(u_MC(:,1:50),0,2);
[u_var_MC100,u_mean_MC100] = var(u_MC(:,1:100),0,2);
[u_var_MC200,u_mean_MC200] = var(u_MC(:,1:200),0,2);
[u_var_MC500,u_mean_MC500] = var(u_MC(:,1:500),0,2);
[u_var_MC1000,u_mean_MC1000] = var(u_MC(:,1:1000),0,2);
[u_var_MC2000,u_mean_MC2000] = var(u_MC(:,1:2000),0,2);
% [u_var_MC5000,u_mean_MC5000] = var(u_MC(:,1:5000),0,2);
% [u_var_MC10000,u_mean_MC10000] = var(u_MC(:,1:10000),0,2);
% [u_var_MC15000,u_mean_MC15000] = var(u_MC(:,1:15000),0,2);
% [u_var_MC20000,u_mean_MC20000] = var(u_MC(:,1:20000),0,2);

% % Moment Eqs
% % [Y_mean,Y_cov,Y_var] = fct_SK(pos_known,val_known,pos_est,sill,range,nugget,m);
Y_mean(1:length(tp),1) = 0; 

scrsz = get(0,'ScreenSize');
resol = 500;

[Xm,Ym] = ndgrid(linspace(min(tp(:,1)), max(tp(:,1)), resol), linspace(min(tp(:,2)), max(tp(:,2)), resol));
Zm = griddata(tp(:,1), tp(:,2), SGS_reordered(:,35), Xm, Ym);
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
Zm = griddata(tp(:,1), tp(:,2), u_MC(:,28), Xm, Ym);
current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
contourf(Xm, Ym, Zm, 20,'LineColor','none','ShowText','off'); grid off;
axis([0 1 0 1]); axis('square');
clim([0.1 1.1]); colorbar;
xlabel('{\it x}_1', 'FontSize', 28); ylabel('{\it x}_2', 'FontSize', 28);
ax = gca; ax.FontSize = 28;
chb = colorbar; ylabel(chb, 'Hydraulic head','FontSize',28, 'FontName', 'Times');
tlm = get(chb,'Limits'); Tspc = linspace(tlm(1),tlm(2),11);
set(chb,'Ticks',Tspc); TL = arrayfun(@(x) sprintf('%.1f',x),Tspc,'un',0);
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

current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
plt_x = tp((tp(:,2) == 0.5),1);
plt_x = [plt_x(2:end); plt_x(1)];
plt_uMC = u_MC((tp(:,2) == 0.5),:);
plt_uMC = [plt_uMC(2:end,:); plt_uMC(1,:)];
plot(plt_x,plt_uMC(:,101:200))
xlabel('{\it x}_1', 'FontSize', 28); ylabel('Hydraulic head', 'FontSize', 28);
text(0.02, 1.03, '{\it x}_2 = 0.5 (MCS #1-100)','FontSize',24)
yticks(0.1:0.2:1.1);
axis([0 1 0.1 1.1]); axis('square');
ax = gca; ax.FontSize = 28;
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')

current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
numPlots = 5;
ramp = linspace(0.75, 0, numPlots);
listOfGrayColors = [ramp; ramp; ramp]';
plt_x = tp((tp(:,2) == 0.5),1);
plt_x = [plt_x(2:end); plt_x(1)];
plt_umMC = u_mean_MC10((tp(:,2) == 0.5),:);
plt_umMC = [plt_umMC(2:end,:); plt_umMC(1,:)];
Nmc10 = plot(plt_x,plt_umMC,'LineWidth',1.6, 'Color', listOfGrayColors(1, :));
hold on
plt_umMC = u_mean_MC20((tp(:,2) == 0.5),:);
plt_umMC = [plt_umMC(2:end,:); plt_umMC(1,:)];
Nmc20 = plot(plt_x,plt_umMC,'LineWidth',1.6, 'Color', listOfGrayColors(2, :));
plt_umMC = u_mean_MC50((tp(:,2) == 0.5),:);
plt_umMC = [plt_umMC(2:end,:); plt_umMC(1,:)];
Nmc50 = plot(plt_x,plt_umMC,'LineWidth',1.6, 'Color', listOfGrayColors(3, :));
plt_umMC = u_mean_MC100((tp(:,2) == 0.5),:);
plt_umMC = [plt_umMC(2:end,:); plt_umMC(1,:)];
Nmc100 = plot(plt_x,plt_umMC,'LineWidth',1.6, 'Color', listOfGrayColors(4, :));
plt_umMC = u_mean_MC200((tp(:,2) == 0.5),:);
plt_umMC = [plt_umMC(2:end,:); plt_umMC(1,:)];
Nmc200 = plot(plt_x,plt_umMC,'LineWidth',1.6, 'Color', listOfGrayColors(5, :));
lgd = legend([Nmc10 Nmc20 Nmc50 Nmc100 Nmc200],{'{\it N}_{\it MC} = 10','{\it N}_{\it MC} = 20','{\it N}_{\it MC} = 50','{\it N}_{\it MC} = 100','{\it N}_{\it MC} = 200'},...
    'Location','southeast');
xlabel('{\it x}_1', 'FontSize', 28); ylabel('Mean of hydraulic head', 'FontSize', 28);
text(0.02, 1.03, '{\it x}_2 = 0.5','FontSize',24)
yticks(0.1:0.2:1.1);
axis([0 1 0.1 1.1]); axis('square');
ax = gca; ax.FontSize = 28;
fontsize(lgd, 22, 'points')
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')

current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
numPlots = 5;
ramp = linspace(0.75, 0, numPlots);
listOfGrayColors = [ramp; ramp; ramp]';
plt_x = tp((tp(:,2) == 0.5),1);
plt_x = [plt_x(2:end); plt_x(1)];
plt_uvMC = u_var_MC100((tp(:,2) == 0.5),:);
plt_uvMC = [plt_uvMC(2:end,:); plt_uvMC(1,:)];
Nvc100 = plot(plt_x,plt_uvMC,'LineWidth',1.6, 'Color', listOfGrayColors(1, :));
hold on
plt_uvMC = u_var_MC200((tp(:,2) == 0.5),:);
plt_uvMC = [plt_uvMC(2:end,:); plt_uvMC(1,:)];
Nvc200 = plot(plt_x,plt_uvMC,'LineWidth',1.6, 'Color', listOfGrayColors(2, :));
plt_uvMC = u_var_MC500((tp(:,2) == 0.5),:);
plt_uvMC = [plt_uvMC(2:end,:); plt_uvMC(1,:)];
Nvc500 = plot(plt_x,plt_uvMC,'LineWidth',1.6, 'Color', listOfGrayColors(3, :));
plt_uvMC = u_var_MC1000((tp(:,2) == 0.5),:);
plt_uvMC = [plt_uvMC(2:end,:); plt_uvMC(1,:)];
Nvc1000 = plot(plt_x,plt_uvMC,'LineWidth',1.6, 'Color', listOfGrayColors(4, :));
plt_uvMC = u_var_MC2000((tp(:,2) == 0.5),:);
plt_uvMC = [plt_uvMC(2:end,:); plt_uvMC(1,:)];
Nvc2000 = plot(plt_x,plt_uvMC,'LineWidth',1.6, 'Color', listOfGrayColors(5, :));
lgd = legend([Nvc100 Nvc200 Nvc500 Nvc1000 Nvc2000],{'{\it N}_{\it MC} = 100','{\it N}_{\it MC} = 200','{\it N}_{\it MC} = 500','{\it N}_{\it MC} = 1000','{\it N}_{\it MC} = 2000'},...
    'Location','northeast');
xlabel('{\it x}_1', 'FontSize', 28); ylabel('Variance of hydraulic head', 'FontSize', 28);
text(0.02, 0.038, '{\it x}_2 = 0.5','FontSize',24)
yticks(0:0.01:0.04);
axis([0 1 0 0.04]); axis('square');
ax = gca; ax.FontSize = 28;
fontsize(lgd, 22, 'points')
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')


current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
plt_x = tp((tp(:,2) == 0),1);
% plt_x = [plt_x(2:end); plt_x(1)];
plt_uMC = u_MC((tp(:,2) == 0),:);
% plt_uMC = [plt_uMC(2:end,:); plt_uMC(1,:)];
plot(plt_x,plt_uMC(:,101:200))
xlabel('{\it x}_1', 'FontSize', 28); ylabel('Hydraulic head', 'FontSize', 28);
text(0.02, 1.03, '{\it x}_2 = 0 (MCS #1-100)','FontSize',24)
yticks(0.1:0.2:1.1);
axis([0 1 0.1 1.1]); axis('square');
ax = gca; ax.FontSize = 28;
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')

current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
numPlots = 5;
ramp = linspace(0.75, 0, numPlots);
listOfGrayColors = [ramp; ramp; ramp]';
plt_x = tp((tp(:,2) == 0),1);
% plt_x = [plt_x(2:end); plt_x(1)];
plt_umMC = u_mean_MC10((tp(:,2) == 0),:);
% plt_umMC = [plt_umMC(2:end,:); plt_umMC(1,:)];
Nmc10 = plot(plt_x,plt_umMC,'LineWidth',1.6, 'Color', listOfGrayColors(1, :));
hold on
plt_umMC = u_mean_MC20((tp(:,2) == 0),:);
% plt_umMC = [plt_umMC(2:end,:); plt_umMC(1,:)];
Nmc20 = plot(plt_x,plt_umMC,'LineWidth',1.6, 'Color', listOfGrayColors(2, :));
plt_umMC = u_mean_MC50((tp(:,2) == 0),:);
% plt_umMC = [plt_umMC(2:end,:); plt_umMC(1,:)];
Nmc50 = plot(plt_x,plt_umMC,'LineWidth',1.6, 'Color', listOfGrayColors(3, :));
plt_umMC = u_mean_MC100((tp(:,2) == 0),:);
% plt_umMC = [plt_umMC(2:end,:); plt_umMC(1,:)];
Nmc100 = plot(plt_x,plt_umMC,'LineWidth',1.6, 'Color', listOfGrayColors(4, :));
plt_umMC = u_mean_MC200((tp(:,2) == 0),:);
% plt_umMC = [plt_umMC(2:end,:); plt_umMC(1,:)];
Nmc200 = plot(plt_x,plt_umMC,'LineWidth',1.6, 'Color', listOfGrayColors(5, :));
lgd = legend([Nmc10 Nmc20 Nmc50 Nmc100 Nmc200],{'{\it N}_{\it MC} = 10','{\it N}_{\it MC} = 20','{\it N}_{\it MC} = 50','{\it N}_{\it MC} = 100','{\it N}_{\it MC} = 200'},...
    'Location','southeast');
xlabel('{\it x}_1', 'FontSize', 28); ylabel('Mean of hydraulic head', 'FontSize', 28);
text(0.02, 1.03, '{\it x}_2 = 0','FontSize',24)
yticks(0.1:0.2:1.1);
axis([0 1 0.1 1.1]); axis('square');
ax = gca; ax.FontSize = 28;
fontsize(lgd, 22, 'points')
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')

current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
numPlots = 5;
ramp = linspace(0.75, 0, numPlots);
listOfGrayColors = [ramp; ramp; ramp]';
plt_x = tp((tp(:,2) == 0.5),1);
plt_x = [plt_x(2:end); plt_x(1)];
plt_uvMC = u_var_MC100((tp(:,2) == 0),:);
% plt_uvMC = [plt_uvMC(2:end,:); plt_uvMC(1,:)];
Nvc100 = plot(plt_x,plt_uvMC,'LineWidth',1.6, 'Color', listOfGrayColors(1, :));
hold on
plt_uvMC = u_var_MC200((tp(:,2) == 0),:);
% plt_uvMC = [plt_uvMC(2:end,:); plt_uvMC(1,:)];
Nvc200 = plot(plt_x,plt_uvMC,'LineWidth',1.6, 'Color', listOfGrayColors(2, :));
plt_uvMC = u_var_MC500((tp(:,2) == 0),:);
% plt_uvMC = [plt_uvMC(2:end,:); plt_uvMC(1,:)];
Nvc500 = plot(plt_x,plt_uvMC,'LineWidth',1.6, 'Color', listOfGrayColors(3, :));
plt_uvMC = u_var_MC1000((tp(:,2) == 0),:);
% plt_uvMC = [plt_uvMC(2:end,:); plt_uvMC(1,:)];
Nvc1000 = plot(plt_x,plt_uvMC,'LineWidth',1.6, 'Color', listOfGrayColors(4, :));
plt_uvMC = u_var_MC2000((tp(:,2) == 0),:);
% plt_uvMC = [plt_uvMC(2:end,:); plt_uvMC(1,:)];
Nvc2000 = plot(plt_x,plt_uvMC,'LineWidth',1.6, 'Color', listOfGrayColors(5, :));
lgd = legend([Nvc100 Nvc200 Nvc500 Nvc1000 Nvc2000],{'{\it N}_{\it MC} = 100','{\it N}_{\it MC} = 200','{\it N}_{\it MC} = 500','{\it N}_{\it MC} = 1000','{\it N}_{\it MC} = 2000'},...
    'Location','northeast');
xlabel('{\it x}_1', 'FontSize', 28); ylabel('Variance of hydraulic head', 'FontSize', 28);
text(0.02, 0.038, '{\it x}_2 = 0','FontSize',24)
yticks(0:0.01:0.04);
axis([0 1 0 0.04]); axis('square');
ax = gca; ax.FontSize = 28;
fontsize(lgd, 22, 'points')
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')

uvMC = zeros(size(u_MC,1),size(u_MC,2));
parfor iMC = 1:size(u_MC,2)
    uvMC(:,iMC) = var(u_MC(:,1:iMC),0,2);
end

current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
Cvgc_1 = uvMC((tp(:,1) == 0.5 & tp(:,2) == 0.5),:);
plot(2:size(Cvgc_1,2),Cvgc_1(2:end),'LineWidth',1.6)
xlabel('{\it N}_{\it MC}', 'FontSize', 28); ylabel('Variance of hydraulic head', 'FontSize', 28);
text(6600, 0.0314, '{\it x}_1 = 0.5, {\it x}_2 = 0.5','FontSize',24)
ylim([0.022 0.032]); axis('square');
ytickformat('%.3f')
ax = gca; ax.FontSize = 28;
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')

current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
Cvgc_2 = uvMC((tp(:,1) == 0.5 & tp(:,2) == 0),:);
plot(2:size(Cvgc_2,2),Cvgc_2(2:end),'LineWidth',1.6)
xlabel('{\it N}_{\it MC}', 'FontSize', 28); ylabel('Variance of hydraulic head', 'FontSize', 28);
text(6680, 0.0314, '{\it x}_1 = 0.5, {\it x}_2 = 0','FontSize',24)
ylim([0.022 0.032]); axis('square');
ytickformat('%.3f')
ax = gca; ax.FontSize = 28;
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')

current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
numPlots = 4;
ramp = linspace(0.75, 0, numPlots);
listOfGrayColors = [ramp; ramp; ramp]';
[f,x] = ecdf(u_MC((tp(:,1) == 0.24 & tp(:,2) == 0.5),1:200));
ncdf200 = plot(x,f,'LineWidth',1.6, 'Color', listOfGrayColors(1, :));
hold on
[f,x] = ecdf(u_MC((tp(:,1) == 0.24 & tp(:,2) == 0.5),1:500));
ncdf500 = plot(x,f,'LineWidth',1.6, 'Color', listOfGrayColors(2, :));
[f,x] = ecdf(u_MC((tp(:,1) == 0.24 & tp(:,2) == 0.5),1:1000));
ncdf1000 = plot(x,f,'LineWidth',1.6, 'Color', listOfGrayColors(3, :));
[f,x] = ecdf(u_MC((tp(:,1) == 0.24 & tp(:,2) == 0.5),1:5000));
ncdf5000 = plot(x,f,'LineWidth',1.6, 'Color', listOfGrayColors(4, :));
hold off
lgd = legend([ncdf200 ncdf500 ncdf1000 ncdf5000],{'{\it N}_{\it MC} = 200','{\it N}_{\it MC} = 500','{\it N}_{\it MC} = 1000','{\it N}_{\it MC} = 5000'},...
    'Location','southeast');
fontsize(lgd, 22, 'points')
text(0.115, 0.95, '{\it x}_1 = 0.24, {\it x}_2 = 0.5','FontSize',24)
xticks(0.1:0.2:1.1);
xlabel('{\it H}')
ylabel('{\it F}_{\ith}({\itH} ; · )')
xlim([0.1 1.1])
ylim([0 1])
ax = gca; ax.FontSize = 28;
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')


current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
numPlots = 4;
ramp = linspace(0.75, 0, numPlots);
listOfGrayColors = [ramp; ramp; ramp]';
[f,x] = ecdf(u_MC((tp(:,1) == 0.5 & tp(:,2) == 0.5),1:200));
ncdf200 = plot(x,f,'LineWidth',1.6, 'Color', listOfGrayColors(1, :));
hold on
[f,x] = ecdf(u_MC((tp(:,1) == 0.5 & tp(:,2) == 0.5),1:500));
ncdf500 = plot(x,f,'LineWidth',1.6, 'Color', listOfGrayColors(2, :));
[f,x] = ecdf(u_MC((tp(:,1) == 0.5 & tp(:,2) == 0.5),1:1000));
ncdf1000 = plot(x,f,'LineWidth',1.6, 'Color', listOfGrayColors(3, :));
[f,x] = ecdf(u_MC((tp(:,1) == 0.5 & tp(:,2) == 0.5),1:5000));
ncdf5000 = plot(x,f,'LineWidth',1.6, 'Color', listOfGrayColors(4, :));
hold off
lgd = legend([ncdf200 ncdf500 ncdf1000 ncdf5000],{'{\it N}_{\it MC} = 200','{\it N}_{\it MC} = 500','{\it N}_{\it MC} = 1000','{\it N}_{\it MC} = 5000'},...
    'Location','southeast');
fontsize(lgd, 22, 'points')
text(0.115, 0.95, '{\it x}_1 = 0.5, {\it x}_2 = 0.5','FontSize',24)
xticks(0.1:0.2:1.1);
xlabel('{\it H}')
ylabel('{\it F}_{\ith}({\itH} ; · )')
xlim([0.1 1.1])
ylim([0 1])
ax = gca; ax.FontSize = 28;
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')

% CDF x1 = 0.3 0.4 0.5 0.6 0.7
current = figure('OuterPosition',[0 0 (scrsz(4)+500) scrsz(4)-100]);
numPlots = 5;
ramp = linspace(0.75, 0, numPlots);
listOfGrayColors = [ramp; ramp; ramp]';
[f,x] = ecdf(u_MC((tp(:,1) == 0.3 & tp(:,2) == 0.5),1:10000));
ncdf0d3 = plot(x,f,'LineWidth',1.6, 'Color', listOfGrayColors(1, :));
hold on
[f,x] = ecdf(u_MC((tp(:,1) == 0.4 & tp(:,2) == 0.5),1:10000));
ncdf0d4 = plot(x,f,'LineWidth',1.6, 'Color', listOfGrayColors(2, :));
[f,x] = ecdf(u_MC((tp(:,1) == 0.5 & tp(:,2) == 0.5),1:10000));
ncdf0d5 = plot(x,f,'LineWidth',1.6, 'Color', listOfGrayColors(3, :));
[f,x] = ecdf(u_MC((tp(:,1) == 0.6 & tp(:,2) == 0.5),1:10000));
ncdf0d6 = plot(x,f,'LineWidth',1.6, 'Color', listOfGrayColors(4, :));
[f,x] = ecdf(u_MC((tp(:,1) == 0.7 & tp(:,2) == 0.5),1:10000));
ncdf0d7 = plot(x,f,'LineWidth',1.6, 'Color', listOfGrayColors(5, :));
hold off
lgd = legend([ncdf0d3 ncdf0d4 ncdf0d5 ncdf0d6 ncdf0d7],{'{\it x}_1 = 0.3','{\it x}_1 = 0.4','{\it x}_1 = 0.5','{\it x}_1 = 0.6','{\it x}_1 = 0.7'},...
    'Location','southeast');
fontsize(lgd, 22, 'points')
text(0.115, 0.95, '{\it x}_2 = 0.5; {\it N}_{\it MC} = 10000','FontSize',24)
xticks(0.1:0.2:1.1);
xlabel('{\it H}')
ylabel('{\it F}_{\ith}({\itH} ; · )')
xlim([0.1 1.1])
ylim([0 1])
ax = gca; ax.FontSize = 28;
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')



% Compute Fh values for evenly spaced H
[f,x] = ecdf(u_MC((tp(:,1) == 0.24 & tp(:,2) == 0.5),1:10000));
x_even = linspace(0.1, 1.1, H_res);
f_even = zeros(size(x_even));
for i = 1:length(x_even)
    idx = find(x >= x_even(i), 1); % Find the first index where x >= x_even(i)
    if isempty(idx)
        f_even(i) = 1; % Assign 1 if no x value is found (beyond the range of x)
    else
        f_even(i) = f(idx); % Assign the corresponding f value
    end
end
FhH_ctps1 = f_even;

[f,x] = ecdf(u_MC((tp(:,1) == 0.5 & tp(:,2) == 0.5),1:10000));
x_even = linspace(0.1, 1.1, H_res);
f_even = zeros(size(x_even));
for i = 1:length(x_even)
    idx = find(x >= x_even(i), 1); % Find the first index where x >= x_even(i)
    if isempty(idx)
        f_even(i) = 1; % Assign 1 if no x value is found (beyond the range of x)
    else
        f_even(i) = f(idx); % Assign the corresponding f value
    end
end
FhH_ctps2 = f_even;

parfor ndi = 1:size(tp,1)
    [f,x] = ecdf(u_MC(ndi,1:10000));
x_even = linspace(0.1, 1.1, H_res);
f_even = zeros(size(x_even));
for i = 1:length(x_even)
    idx = find(x >= x_even(i), 1); % Find the first index where x >= x_even(i)
    if isempty(idx)
        f_even(i) = 1; % Assign 1 if no x value is found (beyond the range of x)
    else
        f_even(i) = f(idx); % Assign the corresponding f value
    end
end
FhH_2601(ndi,:) = f_even;
end

% % MoM CDF
% [alpha,beta] = fct_Steady_Closure_FDM(tpY_cov,tp_tvV,u_mean_MC,u_var_MC,g);
% 
% % tpH nodes
% [~,~,HH] = meshgrid(linspace(0,1,xy_res),linspace(0,1,xy_res),linspace(0.1,1.1,H_res));
% tpH = [repmat(tp,H_res,1) HH(:)];
% alpha_CDF = [repmat(alpha,H_res,1)];
% beta_CDF = [repmat(beta,H_res,1)];
% g_CDF = [repmat(g,H_res,1)];
% % Define Fh boundary conditions
% % Fh(0.1; x,y) = 0
% % Fh(1.1; x,y) = 1
% % Fh(H; x,y) = F_phi(H; x,y)  (x,y) Gamma_D
% % Grad(Fh(H; x,y).n(x,y)) = 0  (x,y) Gamma_N
% tpH(:,4:5) = repmat([Bdps_tv; zeros(size(inps,1),2)],H_res,1);
% % tpH(tpH(:,1) == 0 & tpH(:,2)>=0,5) = 1; % Value of BCs
% % tpH(tpH(:,1) == 1 & tpH(:,2)>=1,5) = 0; % Value of BCs
% tpH((tpH(:,2) == 0),4) = 2; % 2nd BCs
% tpH((tpH(:,2) == 0),5) = 0; % Value of BCs
% tpH((tpH(:,2) == 1),4) = 2; % 2nd BCs
% tpH((tpH(:,2) == 1),5) = 0; % Value of BCs
% tpH((tpH(:,1) == 0),4) = 1; % 1nd BCs
% tpH((tpH(:,1) == 0),5) = 1; % Value of BCs
% tpH((tpH(:,1) == 1),4) = 1; % 1nd BCs
% tpH((tpH(:,1) == 1),5) = 0; % Value of BCs
% tpH((tpH(:,3) == 0.1),4) = 1; % 1st BCs
% tpH((tpH(:,3) == 0.1),5) = 0; % Value of BCs
% tpH((tpH(:,3) == 1.1),4) = 1; % 1st BCs
% tpH((tpH(:,3) == 1.1),5) = 1; % Value of BCs
% tpH(tpH(:,1) == 0,6:8) = repmat([-1 0 0],sum(tpH(:,1) == 0),1); % normal vector
% tpH(tpH(:,1) == 1,6:8) = repmat([1 0 0],sum(tpH(:,1) == 1),1); % normal vector
% tpH((tpH(:,3) == 0.1),6:8) = repmat([0 0 -1],sum(tpH(:,3) == 0.1),1); % normal vector
% tpH((tpH(:,3) == 1.1),6:8) = repmat([0 0 1],sum(tpH(:,3) == 1.1),1); % normal vector
% tpH(tpH(:,2) == 0,6:8) = repmat([0 -1 0],sum(tpH(:,2) == 0),1); % normal vector
% tpH(tpH(:,2) == 1,6:8) = repmat([0 1 0],sum(tpH(:,2) == 1),1); % normal vector
% 
% current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
% scatter3(tpH(:,1), tpH(:,2), tpH(:,3), sz, tpH(:,3), 'filled');
% axis([0 1 0 1 0 1.2]); axis('square');
% box on;
% xlabel('{\it x}_1', 'FontSize', 28); ylabel('{\it x}_2', 'FontSize', 28); zlabel('{\it H}', 'FontSize', 28);
% ax = gca; ax.FontSize = 28;
% set(findall(gcf,'type','text'), 'FontName', 'Times');
% set(gca, 'FontName', 'Times')
% 
% Y_mean_CDF = repmat(Y_mean,H_res,1);
% u_mean_CDF = repmat(u_mean_ME,H_res,1);
% 
% % GFDM coef 3D
% ns = 33;
% order = 2;
% [Wcoef_3D,idx] = fct_GFDM_Coef_3D(order,tpH(:,1),tpH(:,2),tpH(:,3),ns);
% 
% % GE: K(4+7+9) + [alpha(x,y)(H-h(x,y))+beta(x,y)+g(x,y)](3) = 0
% Fh = fct_Steady_CdfHead(Wcoef_3D,idx,ns,tpH,alpha_CDF,beta_CDF,g_CDF,Y_mean_CDF,u_mean_CDF);
% 
% % figure
% parfor refer = 1:size(tp,1)
%     pick = find(tpH(:,1) == tp(refer,1) & tpH(:,2) == tp(refer,2));
%     cdf(refer,:) = Fh(pick);
% end
% % plot(cdf')
% 
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
% current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
% x1 = 0.25;
% x2 = 0.5;
% pick = tp(:,1) == x1 & tp(:,2) == x2;
% plot(linspace(0.1,1.1,H_res), cdf(pick,:),":r",'LineWidth',4)
% xlabel('{\it H}', 'FontSize', 28); ylabel('{\it F_h}({\it H}; {\it x}_1= 0.25,{\it x}_2= 0.5)', 'FontSize', 28);
% ax = gca; ax.FontSize = 28;
% set(findall(gcf,'type','text'), 'FontName', 'Times');
% set(gca, 'FontName', 'Times')
% 
% current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
% x1 = 0.5;
% x2 = 0.5;
% pick = tp(:,1) == x1 & tp(:,2) == x2;
% plot(linspace(0.1,1.1,H_res), cdf(pick,:),":r",'LineWidth',4)
% xlabel('{\it H}', 'FontSize', 28); ylabel('{\it F_h}({\it H}; {\it x}_1= 0.5,{\it x}_2= 0.5)', 'FontSize', 28);
% ax = gca; ax.FontSize = 28;
% set(findall(gcf,'type','text'), 'FontName', 'Times');
% set(gca, 'FontName', 'Times')
