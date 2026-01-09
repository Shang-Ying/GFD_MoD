clear; close all; clc;

% Boundary Nodes
BCNodes = readtable("FlowDomainVertices.csv");

% Monthly Average Head Observation in May 2021 during the severe drought
Heads = readtable("Parameters.xlsx");
% Monthly Average Head Observation in May 2021 during the severe drought
RiverBed = readtable("RiverBedHeight.xlsx");

% Geostatistics of Transmissivity
pos_T = [Heads.hX Heads.hY];
val_T = Heads.T_m2_s;
pos_T = pos_T(~isnan(val_T));
val_T = val_T(~isnan(val_T));
[gamma, sepa] = semivar_exp(pos_T,val_T); % data points do not provide enough structure

% Save a text file for loading using make_domain
writematrix([BCNodes.X BCNodes.Y],'NBasin.txt')

% Discrete boundary smoothing
[b] = make_domain('NBasin.txt');
hbdy   = 400; % density
[bdy] = bsmooth(b.xy, hbdy);
b_grad = [-1*gradient(bdy(:,2)) gradient(bdy(:,1))];
b_nv = zeros(size(bdy,1),1);
for i = 1:size(bdy,1)
    b_nv(i,1:2) = b_grad(i,:)./norm(b_grad(i,:));
end

% % Examine normal vectors
% scrsz = get(0,'ScreenSize');
% current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
% Lg_BDY = plot([bdy(:,1);bdy(1,1)],[bdy(:,2);bdy(1,2)],'-','LineWidth',3); hold on;
% Lg_NV = quiver(bdy(:,1),bdy(:,2),b_nv(:,1),b_nv(:,2),'r','LineWidth',1.5,'AutoScaleFactor',0.3);
% legend([Lg_BDY Lg_NV],{'Domain boundary','Normal vector'},'Location','southeast','FontSize',18);
% % axis([0 5000 0 5000]);
% axis('square');
% xlabel('{\it x} (m)', 'FontSize', 28); ylabel('{\it y} (m)', 'FontSize', 28);
% ax = gca; ax.FontSize = 28;
% set(findall(gcf,'type','text'), 'FontName', 'Times');
% set(gca, 'FontName', 'Times')
% hold off;

% (Universal) Kriging for Constant Head
pos_est = bdy;
xyh_obs = [Heads.hX Heads.hY Heads.May2021Drought; RiverBed.E RiverBed.N RiverBed.Height];
pos_known = xyh_obs(:,1:2);
val_known = xyh_obs(:,3);
V = '1 Sph(800000)';
% options.polytrend = [2 1];
[d_est,d_var] = krig(pos_known,val_known,pos_est,V);

% Control Points
Wells = readtable("Wells.csv");

% Define prescribed head boundary indices
S = shaperead('BC_Diri\BC_Diri_MATLAB.shp');
in_dajia = inpolygon(pos_est(:,1), pos_est(:,2), S(1).X, S(1).Y); % Use inpolygon to check which points are inside current polygon
in_dali = inpolygon(pos_est(:,1), pos_est(:,2), S(2).X, S(2).Y); % Use inpolygon to check which points are inside current polygon
in_any_polygon = in_dajia | in_dali; % Update master list of points inside any polygon
nodes_dajia = find(in_dajia);
nodes_dali = find(in_dali);
nodes_Diri = find(in_any_polygon); % Get indices of nodes inside polygons

% Dimensionless spatial coordinates
x = bdy(:,1);
y = bdy(:,2);
y_range = (max(y) - min(y));
x_min = min(x);
y_min = min(y);
x_dimless = (x - x_min) / y_range;  % Shift so min is 0, then scale
y_dimless = (y - y_min) / y_range;  % Shift so min is 0, then scale
bdy_dimless = [x_dimless, y_dimless];
x_ctps = (Wells.E - x_min) / y_range;
y_ctps = (Wells.N - y_min) / y_range; 
ctps_dimless = [x_ctps, y_ctps];

% Dimensionless hydraulic head
u_range = max(d_est) - min(d_est);
d_est_dimless = (d_est - min(d_est)) / u_range;

% Boundary nodes types and values
BC_type = 2*ones(size(bdy,1),1);
BC_type(nodes_Diri) = 1;
BC_value = zeros(size(bdy,1),1);
BC_value(nodes_Diri) = d_est_dimless(nodes_Diri);

% Advancing Front Node Generation
box_xm = 0;
box_xM = 1;
box_ym = 0;
box_yM = 1;
b_box = [box_xm box_xM box_ym box_yM];
ninit = 1001;
dotmax = 5e6;
ptol = 0.020;
btol = 0.010;
radius = @(p,ctps_dimless,bdy_dimless) min((ptol/3 + 0.15*(min(pdist2(ctps_dimless, p)))),(btol/4 + 0.25*(min(pdist2(bdy_dimless, p)))));
xy = node_drop_2d_ctps_tol(b_box, bdy_dimless, ninit, dotmax, radius, ctps_dimless, btol/1.8, ptol/3);
[in,on] = inpolygon(xy(:,1),xy(:,2),bdy_dimless(:,1),bdy_dimless(:,2));
xy = xy(in,:);

% Generalized Finite Difference Coefficient
tp = [bdy_dimless; xy];
order = 2;
ns = 15;
[Wcoef,idx] = fct_GFDM_Coef(order,tp(:,1),tp(:,2),ns);

% Moment Equations
Y_mean(1:length(tp),1) = 0;
sill = 1.0;
cor_length = 0.02;
range = cor_length * 3;
nugget = 0;
[Y_cov] = fct_UnCovMat(tp,sill,range,nugget);
tp_cov = [tp(:,1:2) Y_mean Y_cov]; % combine with kriging Y covariance
Bdps = [bdy_dimless BC_type BC_value b_nv];
tp_tvV = [Bdps(:,:); [xy zeros(size(xy,1),size(Bdps,2)-size(xy,2))]];
[G_loop,D_loop,N_loop] = fct_2DStencil_GFDM(tp_tvV);
u_mean_ME = fct_Steady_1stMEHead_HLG(Wcoef,idx,ns,tp_cov,tp_tvV,G_loop,D_loop,N_loop);
[CYY,CYu,Cuu] = fct_Steady_2ndMEHead_HLG(Wcoef,idx,ns,tp_cov,tp_tvV,u_mean_ME,G_loop,D_loop,N_loop);
u_var_ME(:,1) = diag(Cuu);

% Closure Variables
K_mean = exp(Y_mean);
V = zeros(length(tp),1);
g(length(tp),1) = 0;
beta = zeros(length(tp),1);
alpha = zeros(length(tp),1);
parfor ith =  1:length(tp)
    V(ith,1) = -0.5.*K_mean(ith).*...
        ((Wcoef(3,1:ns,ith)*u_var_ME((idx(ith,1:ns)),1)) + ...
        (Wcoef(5,1:ns,ith)*u_var_ME((idx(ith,1:ns)),1)));
    beta(ith,1) = K_mean(ith).*...
        ((Wcoef(3,1:ns,ith)*u_mean_ME((idx(ith,1:ns)),1)) + ...
        (Wcoef(5,1:ns,ith)*u_mean_ME((idx(ith,1:ns)),1))) +...
        g(ith,1);
    alpha(ith,1) = (K_mean(ith).*...
        (((Wcoef(1,1:ns,ith)*u_mean_ME((idx(ith,1:ns)),1))^2 + ...
        (Wcoef(2,1:ns,ith)*u_mean_ME((idx(ith,1:ns)),1))^2)) -...
        V(ith,1)) ./ u_var_ME(ith,1);
end

% Parameters for solving Fh
H_res = 51;
[~,HH] = meshgrid(linspace(0,1,size(tp,1)),linspace(0,1,H_res));
HH = HH';
tp_Fh = [repmat(tp,H_res,1) HH(:)];
Bot_logic = tp_Fh(:,3) == 0;
Top_logic = tp_Fh(:,3) == 1;
tp_type = tp_tvV(:,3);
tp_Fh_type = [repmat(tp_type,H_res,1)];
tp_Fh_type(Bot_logic) = 1;
tp_Fh_type(Top_logic) = 1;
tp_value = tp_tvV(:,4);
tp_DiriHead = [repmat(tp_value,H_res,1)];
% tp_value(in_dajia) = 0; % wrong
% tp_value(in_dali) = 1; % wrong
tp_in_dajia = zeros(size(tp,1),1);
tp_in_dajia(1:size(in_dajia,1)) = in_dajia;
tp_Fh_in_dajia = logical([repmat(tp_in_dajia,H_res,1)]);
tp_in_dali = zeros(size(tp,1),1);
tp_in_dali(1:size(in_dali,1)) = in_dali;
tp_Fh_in_dali = logical([repmat(tp_in_dali,H_res,1)]);
tp_Fh_value = zeros(size(tp_Fh,1),1);
tp_Fh_value(Bot_logic) = 0;
tp_Fh_value(Top_logic) = 1;
tp_Fh_value(tp_Fh_in_dajia) = tp_Fh(tp_Fh_in_dajia,3) >= tp_DiriHead(tp_Fh_in_dajia,1);
tp_Fh_value(tp_Fh_in_dali) = tp_Fh(tp_Fh_in_dali,3) >= tp_DiriHead(tp_Fh_in_dali,1);
tp_Fh_nv1 = [repmat(tp_tvV(:,5),H_res,1)];
tp_Fh_nv1(Bot_logic) = 0;
tp_Fh_nv1(Top_logic) = 0;
tp_Fh_nv2 = [repmat(tp_tvV(:,6),H_res,1)];
tp_Fh_nv2(Bot_logic) = 0;
tp_Fh_nv2(Top_logic) = 0;
tp_Fh_nv3 = zeros(size(tp_Fh,1),1);
tp_Fh_nv3(Bot_logic) = -1;
tp_Fh_nv3(Top_logic) = 1;
beta_Fh = [repmat(beta,H_res,1)];
alpha_Fh = [repmat(alpha,H_res,1)];
g_Fh = [repmat(g,H_res,1)];
Y_mean_Fh = [repmat(Y_mean,H_res,1)];
u_mean_Fh = [repmat(u_mean_ME,H_res,1)];
tpH = [tp_Fh tp_Fh_type tp_Fh_value tp_Fh_nv1 tp_Fh_nv2 tp_Fh_nv3];

% GFDM 3D Coefficients, Define ns_Fh
ns_Fh = 27;
[Wcoef_3D,idx_3D] = fct_GFDM_Coef_3D(order,tpH(:,1),tpH(:,2),tpH(:,3),ns_Fh);

% CdfHead_Taichung
Fh = fct_Steady_CdfHead(Wcoef_3D,idx_3D,ns_Fh,tpH,alpha_Fh,beta_Fh,g_Fh,Y_mean_Fh,u_mean_Fh);

%% Computational Node Set
scrsz = get(0,'ScreenSize');
current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
Lg_ND = plot(tp(:,1)*y_range+x_min, tp(:,2)*y_range+y_min,'.','color',[0.8 0.8 0.8],'MarkerSize', 12); axis('square'); hold on
Lg_DBC = plot(tp(in_any_polygon,1)*y_range+x_min, tp(in_any_polygon,2)*y_range+y_min,'or','MarkerSize',4.2,'MarkerFaceColor',[1, 0, 0]);
Lg_NBC = plot(tp((~in_any_polygon),1)*y_range+x_min, tp((~in_any_polygon),2)*y_range+y_min,'oc','MarkerSize',4.2,'MarkerFaceColor',[0.3010, 0.7450, 0.9330]); hold on
Lg_Well = plot(ctps_dimless(:,1)*y_range+x_min, ctps_dimless(:,2)*y_range+y_min,'s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6]); hold on
hold off;
axis equal
legend([Lg_Well Lg_ND Lg_DBC Lg_NBC],{'Designated wells','Meshless nodes','Dirichlet BC','Neumann BC'},'Location','northwest','FontSize',20);
xlabel('{\it x} (m)', 'FontSize', 28); ylabel('{\it y} (m)', 'FontSize', 28);
ax = gca;
xtickformat('%.0f') % Force whole numbers
ytickformat('%.0f') % Force whole numbers
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;
ax.YTickLabelRotation = 90;
tick_min = ceil(min(ylim)/5000) * 5000; % Round up to the nearest 5000
tick_max = floor(max(ylim)/5000) * 5000; % Round down to the nearest 5000
yticks(tick_min:5000:tick_max); % Define Y-ticks with 5000 interval
ax.FontSize = max(28, min(14, 0.02 * ax.Position(3)*ax.Position(4))); % Adjust based on figure size
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')

%% Riverbed Elevation
scrsz = get(0,'ScreenSize');
current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
Lg_dajia = plot(tp(in_dajia,1)*y_range+x_min, BC_value(in_dajia,1)*u_range+min(d_est),'ok','MarkerSize',8.2,'MarkerFaceColor',[0, 0, 0]); hold on
Lg_dali = plot(tp(in_dali,1)*y_range+x_min, BC_value(in_dali,1)*u_range+min(d_est),'ok','MarkerSize',8.2,'MarkerFaceColor',[1, 1, 1]);
hold off;
margin = 0.05 * (max(tp(in_dajia,1)*y_range+x_min) - min(tp(in_dali,1)*y_range+x_min)); % 5% margin
xlim([min(tp(in_dali,1)*y_range+x_min) - margin, max(tp(in_dajia,1)*y_range+x_min) + margin]);
legend([Lg_dajia Lg_dali],{'Dajia River','Dali River'},'Location','northwest','FontSize',24);
xlabel('{\it x} (m)', 'FontSize', 28); ylabel('Riverbed Elevation (m)', 'FontSize', 28);
ax = gca;
xtickformat('%.0f') % Force whole numbers
ax.XAxis.Exponent = 0;
ax.FontSize = max(28, min(18, 0.02 * ax.Position(3)*ax.Position(4))); % Adjust based on figure size
tick_min = ceil(min(xlim)/5000) * 5000; % Round up to the nearest 5000
tick_max = floor(max(xlim)/5000) * 5000; % Round down to the nearest 5000
xticks(tick_min:5000:tick_max); % Define X-ticks with 5000 interval
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')

%% Mean of head
scrsz = get(0,'ScreenSize');
resol = 800;
[Xm,Ym] = ndgrid(linspace(min(tp(:,1)*y_range+x_min), max(tp(:,1)*y_range+x_min), resol), linspace(min(tp(:,2)*y_range+y_min), max(tp(:,2)*y_range+y_min), resol*1.2));
ShpLogic = (~reshape(inpolygon(Xm(:),Ym(:),bdy_dimless(:,1)*y_range+x_min,bdy_dimless(:,2)*y_range+y_min),size(Xm)));
Zm = griddata(tp(:,1)*y_range+x_min, tp(:,2)*y_range+y_min, u_mean_ME*u_range+min(d_est), Xm, Ym);
Zm(ShpLogic) = nan;
current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
contourf(Xm, Ym, Zm, 30,'LineColor','none','ShowText','off'); grid off; hold on
plot(bdy_dimless(:,1)*y_range+x_min, bdy_dimless(:,2)*y_range+y_min, 'Color', [0.5 0.5 0.5], 'LineWidth', 1); 
axis equal
colorbar;
fullColormap = othercolor('GnBu7');
colormap(fullColormap);
z_min = min(Zm(:), [], 'omitnan');
z_max = max(Zm(:), [], 'omitnan');
clim([z_min z_max]);
chb = colorbar; 
ylabel(chb, 'Mean of hydraulic head (m)','FontSize',28, 'FontName', 'Times');
tick_interval = 20;
start_tick = ceil(z_min / tick_interval) * tick_interval;  % First tick at or above min
end_tick = floor(z_max / tick_interval) * tick_interval;   % Last tick at or below max
Tspc = start_tick:tick_interval:end_tick;
set(chb,'Ticks',Tspc); 
TL = arrayfun(@(x) sprintf('%.0f',x), Tspc, 'un', 0); % Whole numbers
set(chb,'TickLabels',TL)
xlabel('{\it x} (m)', 'FontSize', 28); ylabel('{\it y} (m)', 'FontSize', 28);
ax = gca;
xtickformat('%.0f') % Force whole numbers
ytickformat('%.0f') % Force whole numbers
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;
ax.YTickLabelRotation = 90;
tick_min = ceil(min(ylim)/5000) * 5000; % Round up to the nearest 5000
tick_max = floor(max(ylim)/5000) * 5000; % Round down to the nearest 5000
yticks(tick_min:5000:tick_max); % Define Y-ticks with 5000 interval
ax.FontSize = max(28, min(14, 0.02 * ax.Position(3)*ax.Position(4))); % Adjust based on figure size
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')


%% Standard deviation of head
scrsz = get(0,'ScreenSize');
resol = 800;
[Xm,Ym] = ndgrid(linspace(min(tp(:,1)*y_range+x_min), max(tp(:,1)*y_range+x_min), resol), linspace(min(tp(:,2)*y_range+y_min), max(tp(:,2)*y_range+y_min), resol*1.2));
ShpLogic = (~reshape(inpolygon(Xm(:),Ym(:),bdy_dimless(:,1)*y_range+x_min,bdy_dimless(:,2)*y_range+y_min),size(Xm)));
Zm = griddata(tp(:,1)*y_range+x_min, tp(:,2)*y_range+y_min, sqrt(u_var_ME*(u_range^2)), Xm, Ym);
Zm(ShpLogic) = nan;
current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
contourf(Xm, Ym, Zm, 30,'LineColor','none','ShowText','off'); grid off; hold on
plot(bdy_dimless(:,1)*y_range+x_min, bdy_dimless(:,2)*y_range+y_min, 'Color', [0.5 0.5 0.5], 'LineWidth', 1); 
axis equal
colorbar;
fullColormap = othercolor('RdPu9');
colormap(fullColormap);
z_min = min(Zm(:), [], 'omitnan');
z_max = max(Zm(:), [], 'omitnan');
clim([0 z_max]);
chb = colorbar; 
ylabel(chb, 'Standard deviation of hydraulic head (m)','FontSize',28, 'FontName', 'Times');
tick_interval = 2;
start_tick = ceil(0 / tick_interval) * tick_interval;  % First tick at or above min
end_tick = floor(z_max / tick_interval) * tick_interval;   % Last tick at or below max
Tspc = start_tick:tick_interval:end_tick;
set(chb,'Ticks',Tspc); 
TL = arrayfun(@(x) sprintf('%.1f',x), Tspc, 'un', 0);
set(chb,'TickLabels',TL)
xlabel('{\it x} (m)', 'FontSize', 28); ylabel('{\it y} (m)', 'FontSize', 28);
ax = gca;
xtickformat('%.0f') % Force whole numbers
ytickformat('%.0f') % Force whole numbers
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;
ax.YTickLabelRotation = 90;
tick_min = ceil(min(ylim)/5000) * 5000; % Round up to the nearest 5000
tick_max = floor(max(ylim)/5000) * 5000; % Round down to the nearest 5000
yticks(tick_min:5000:tick_max); % Define Y-ticks with 5000 interval
ax.FontSize = max(28, min(14, 0.02 * ax.Position(3)*ax.Position(4))); % Adjust based on figure size
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')

%% CDF (2+1)-d domain and BCs
H_min = 0;
H_max = 240;
visualize_complex_extruded_shape(bdy, BC_type, H_min, H_max)

%% Nodal solution of Fh
scrsz = get(0,'ScreenSize');
current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
sz = 5;
scatter3(tpH(:,1)*y_range+x_min, tpH(:,2)*y_range+y_min, tpH(:,3)*u_range+min(d_est), sz, Fh, 'filled','MarkerEdgeColor','k',...
    'MarkerFaceAlpha',1.0,'MarkerEdgeAlpha',.05);
view_az = 27;    % New azimuth angle
view_el = 60;    % New elevation angle
view(view_az, view_el);
grid on;
axis equal;
daspect auto
ax = gca;
pos = ax.Position;
pos(1) = pos(1) - 0.00;
pos(3) = pos(3) - 0.15; % reduce width by 5% (adjust as needed)
ax.Position = pos;
xtickformat('%.0f') % Force whole numbers
ytickformat('%.0f') % Force whole numbers
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;
% ax.YTickLabelRotation = 60;
tick_min = ceil(min(ylim)/5000) * 5000; % Round up to the nearest 5000
tick_max = floor(max(ylim)/5000) * 5000; % Round down to the nearest 5000
yticks(tick_min:5000:tick_max); % Define Y-ticks with 5000 interval
tick_min = ceil(min(zlim)/100) * 100; % Round up to the nearest 5000
tick_max = floor(max(zlim)/100) * 100; % Round down to the nearest 5000
zticks(tick_min:100:tick_max); % Define Y-ticks with 5000 interval
xlabel('{\it x} (m)', 'FontSize', 24);
ylabel('{\it y} (m)', 'FontSize', 24);
zlabel('{\it H} (m)', 'FontSize', 24);
ax.FontSize = max(20, min(14, 0.02 * ax.Position(3)*ax.Position(4))); % Adjust based on figure size
% colorbar('eastoutside');
fullColormap = othercolor('Spectral7');
nColors = size(fullColormap, 1); % Total number of colors in the colormap
upperHalf = fullColormap(ceil(nColors / 2):end, :); % Extract the upper half
colormap(flipud(upperHalf));
clim([0 1]);
chb = colorbar('eastoutside'); ylabel(chb, 'Cumulative distribution function, {\it F}_{\ith}','FontSize',24, 'FontName', 'Times');
chb.Position = [0.85, 0.1, 0.02, 0.8];
tlm = get(chb,'Limits'); Tspc = linspace(0,1,6);
set(chb,'Ticks',Tspc); TL = arrayfun(@(x) sprintf('%.1f',x),Tspc,'un',0);
set(chb,'TickLabels',TL)
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')

%% Exceedance map
% P[h(x) > H = 0.8] = 1 − Fh(H = 0.8; x)
H = 0.40;
select = abs(tpH(:,3) - H) < eps;
h_cut = H*u_range+min(d_est);
tp_select = tpH(select,:);
scrsz = get(0,'ScreenSize');
current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
[Xm,Ym] = ndgrid(linspace(min(tp_select(:,1)*y_range+x_min), max(tp_select(:,1)*y_range+x_min), resol), linspace(min(tpH(select,2)*y_range+y_min), max(tpH(select,2)*y_range+y_min), resol*1.2));
ShpLogic = (~reshape(inpolygon(Xm(:),Ym(:),bdy_dimless(:,1)*y_range+x_min,bdy_dimless(:,2)*y_range+y_min),size(Xm)));
Zm = griddata(tp_select(:,1)*y_range+x_min, tp_select(:,2)*y_range+y_min, 1-Fh(select,1), Xm, Ym);
Zm(ShpLogic) = nan;
contourf(Xm, Ym, Zm, 20,'LineColor','none','ShowText','off'); grid off; hold on;
plot(bdy_dimless(:,1)*y_range+x_min, bdy_dimless(:,2)*y_range+y_min, 'Color', [0.5 0.5 0.5], 'LineWidth', 1); 
axis equal
ax = gca; ax.FontSize = 28;
chb = colorbar; 
fullColormap = othercolor('Spectral7');
nColors = size(fullColormap, 1);       % Total number of colors in the colormap
upperHalf = fullColormap(ceil(nColors / 2):end, :); % Extract the upper half
colormap(flipud(upperHalf));
clim([0 1]);
ylabel(chb, 'P [{\it h} > 100 m] ','FontSize',28, 'FontName', 'Times');
tick_interval = 0.2;
start_tick = ceil(0 / tick_interval) * tick_interval;  % First tick at or above min
end_tick = floor(1 / tick_interval) * tick_interval;   % Last tick at or below max
Tspc = start_tick:tick_interval:end_tick;
set(chb,'Ticks',Tspc); 
TL = arrayfun(@(x) sprintf('%.1f',x), Tspc, 'un', 0);
set(chb,'TickLabels',TL)
xlabel('{\it x} (m)', 'FontSize', 28); ylabel('{\it y} (m)', 'FontSize', 28);
ax = gca;
xtickformat('%.0f') % Force whole numbers
ytickformat('%.0f') % Force whole numbers
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;
ax.YTickLabelRotation = 90;
tick_min = ceil(min(ylim)/5000) * 5000; % Round up to the nearest 5000
tick_max = floor(max(ylim)/5000) * 5000; % Round down to the nearest 5000
yticks(tick_min:5000:tick_max); % Define Y-ticks with 5000 interval
ax.FontSize = max(28, min(14, 0.02 * ax.Position(3)*ax.Position(4))); % Adjust based on figure size
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')
hold off

%% CDF at wells
for i = 1:size(ctps_dimless,1)
    current = figure('OuterPosition',[0 0 (scrsz(4)+100)/1.5 scrsz(4)/1.5]);
    N_well = i;
    x1 = ctps_dimless(N_well,1);
    x2 = ctps_dimless(N_well,2);
    pick = tpH(:,1) == x1 & tpH(:,2) == x2;
    CDF(:,i) = [Fh(pick)];
    uq = plot(tpH(pick,3)*u_range+min(d_est), CDF(:,i),"-k",'LineWidth',2);
    axis([min(u_mean_ME*u_range+min(d_est)) max(u_mean_ME*u_range+min(d_est)) 0 1]);
    text(180,0.1,Wells.myid(i), 'FontName', 'Times', 'FontSize', 28)
    tick_min = ceil(min(xlim)/50) * 50; % Round up to the nearest 5000
    tick_max = floor(max(xlim)/50) * 50; % Round down to the nearest 5000
    xticks(tick_min:50:tick_max); % Define Y-ticks with 5000 interval
    xlabel('{\it H}', 'FontSize', 28);
    % ylabel("{\it F_h}({\it H}; {\it x} = " + num2str(x1*y_range+x_min) + ", {\it y} = " + num2str(x2*y_range+y_min) + ")", 'FontSize', 28);
    ylabel("{\it F_h}({\it H}; {\it x}, {\it y})", 'FontSize', 28);
    ax = gca; ax.FontSize = 28;
    set(findall(gcf,'type','text'), 'FontName', 'Times');
    set(gca, 'FontName', 'Times')
end
%% CDF Together
colors = [
    0 0.4470 0.7410;    % blue
    0.8500 0.3250 0.0980;  % orange
    0.9290 0.6940 0.1250;  % yellow
    0.4940 0.1840 0.5560;  % purple
    0.4660 0.6740 0.1880;  % green
    0.3010 0.7450 0.9330;  % light blue
    0.6350 0.0780 0.1840;  % burgundy
    0.75 0.75 0;          % olive
];
lineStyles = {'-', '-', '-', '-', '-', ':', ':', ':'};
figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
hold on;
for i = 1:size(ctps_dimless,1)
    N_well = i;
    x1 = ctps_dimless(N_well,1);
    x2 = ctps_dimless(N_well,2);
    pick = tpH(:,1) == x1 & tpH(:,2) == x2;
    CDF(:,i) = [Fh(pick)];
    plot(tpH(pick,3)*u_range+min(d_est), CDF(:,i), lineStyles{i}, 'LineWidth', 2, 'Color', colors(i,:));
end
axis([min(u_mean_ME*u_range+min(d_est)) max(u_mean_ME*u_range+min(d_est)) 0 1]);
tick_min = ceil(min(xlim)/50) * 50;
tick_max = floor(max(xlim)/50) * 50;
xticks(tick_min:50:tick_max);
xlabel('{\it H}', 'FontSize', 28);
ylabel("{\it F_h}({\it H}; {\it x}, {\it y})", 'FontSize', 28);
legend(Wells.myid, 'Location', 'southeast', 'FontSize', 24);
ax = gca; ax.FontSize = 28;
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times');
grid off; box on;
hold off;
%% Standard deviation of head + Wells
scrsz = get(0,'ScreenSize');
resol = 800;
[Xm,Ym] = ndgrid(linspace(min(tp(:,1)*y_range+x_min), max(tp(:,1)*y_range+x_min), resol), linspace(min(tp(:,2)*y_range+y_min), max(tp(:,2)*y_range+y_min), resol*1.2));
ShpLogic = (~reshape(inpolygon(Xm(:),Ym(:),bdy_dimless(:,1)*y_range+x_min,bdy_dimless(:,2)*y_range+y_min),size(Xm)));
Zm = griddata(tp(:,1)*y_range+x_min, tp(:,2)*y_range+y_min, sqrt(u_var_ME*(u_range^2)), Xm, Ym);
Zm(ShpLogic) = nan;
current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
contourf(Xm, Ym, Zm, 30,'LineColor','none','ShowText','off'); grid off; hold on
plot(bdy_dimless(:,1)*y_range+x_min, bdy_dimless(:,2)*y_range+y_min, 'Color', [0.5 0.5 0.5], 'LineWidth', 1); 
Lg_Well = plot(ctps_dimless(:,1)*y_range+x_min, ctps_dimless(:,2)*y_range+y_min,'s','MarkerSize',16,'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6]); hold on
axis equal
colorbar;
fullColormap = othercolor('RdPu9');
colormap(fullColormap);
z_min = min(Zm(:), [], 'omitnan');
z_max = max(Zm(:), [], 'omitnan');
clim([0 z_max]);
chb = colorbar; 
ylabel(chb, 'Standard deviation of hydraulic head (m)','FontSize',28, 'FontName', 'Times');
tick_interval = 2;
start_tick = ceil(0 / tick_interval) * tick_interval;  % First tick at or above min
end_tick = floor(z_max / tick_interval) * tick_interval;   % Last tick at or below max
Tspc = start_tick:tick_interval:end_tick;
set(chb,'Ticks',Tspc); 
TL = arrayfun(@(x) sprintf('%.1f',x), Tspc, 'un', 0);
set(chb,'TickLabels',TL)
xlabel('{\it x} (m)', 'FontSize', 28); ylabel('{\it y} (m)', 'FontSize', 28);
ax = gca;
xtickformat('%.0f') % Force whole numbers
ytickformat('%.0f') % Force whole numbers
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;
ax.YTickLabelRotation = 90;
tick_min = ceil(min(ylim)/5000) * 5000; % Round up to the nearest 5000
tick_max = floor(max(ylim)/5000) * 5000; % Round down to the nearest 5000
yticks(tick_min:5000:tick_max); % Define Y-ticks with 5000 interval
ax.FontSize = max(28, min(14, 0.02 * ax.Position(3)*ax.Position(4))); % Adjust based on figure size
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')
%% Depth PDF
H_cord = (HH(1,:)')*u_range+min(d_est);
scrsz = get(0,'ScreenSize');
current = figure('OuterPosition',[0 0 (1.5*scrsz(4)) scrsz(4)]); hold on;
for i = 1:8
    top_limit_H(i) = Wells.Elevation(i);
end
% Convert CDF to PDF by taking the derivative (using diff)
pdf1 = zeros(size(H_cord));
pdf2 = zeros(size(H_cord));
pdf3 = zeros(size(H_cord));
pdf4 = zeros(size(H_cord));
pdf5 = zeros(size(H_cord));
pdf6 = zeros(size(H_cord));
pdf7 = zeros(size(H_cord));
pdf8 = zeros(size(H_cord));
for i = 2:length(H_cord)
    pdf1(i) = (CDF(i,1) - CDF(i-1,1)) / (H_cord(i) - H_cord(i-1));
    pdf2(i) = (CDF(i,2) - CDF(i-1,2)) / (H_cord(i) - H_cord(i-1));
    pdf3(i) = (CDF(i,3) - CDF(i-1,3)) / (H_cord(i) - H_cord(i-1));
    pdf4(i) = (CDF(i,4) - CDF(i-1,4)) / (H_cord(i) - H_cord(i-1));
    pdf5(i) = (CDF(i,5) - CDF(i-1,5)) / (H_cord(i) - H_cord(i-1));
    pdf6(i) = (CDF(i,6) - CDF(i-1,6)) / (H_cord(i) - H_cord(i-1));
    pdf7(i) = (CDF(i,7) - CDF(i-1,7)) / (H_cord(i) - H_cord(i-1));
    pdf8(i) = (CDF(i,8) - CDF(i-1,8)) / (H_cord(i) - H_cord(i-1));
end
pdf1 = pdf1 / max(pdf1(pdf1 > 0));  % Ignore zeros when finding max
pdf2 = pdf2 / max(pdf2(pdf2 > 0));
pdf3 = pdf3 / max(pdf3(pdf3 > 0));
pdf4 = pdf4 / max(pdf4(pdf4 > 0));
pdf5 = pdf5 / max(pdf5(pdf5 > 0));
pdf6 = pdf6 / max(pdf6(pdf6 > 0));
pdf7 = pdf7 / max(pdf7(pdf7 > 0));
pdf8 = pdf8 / max(pdf8(pdf8 > 0));
pdf_names = {'EW01', 'EW02', 'EW03', 'EW04', 'EW05', 'EW06', 'EW07', 'EW08'};
total_width = 12;
num_pdfs = 8;
spacing = total_width / (num_pdfs + 1);
offsets = spacing:spacing:total_width-spacing;
offsets = offsets+0.3*spacing;
colors = [
    0 0.4470 0.7410;  % blue
    0.8500 0.3250 0.0980;  % orange
    0.9290 0.6940 0.1250;  % yellow
    0.4940 0.1840 0.5560;  % purple
    0.4660 0.6740 0.1880;  % green
    0.3010 0.7450 0.9330;  % light blue
    0.6350 0.0780 0.1840;  % burgundy
    0.75 0.75 0;  % olive
];
all_pdfs = {pdf1, pdf2, pdf3, pdf4, pdf5, pdf6, pdf7, pdf8};
padding = max(top_limit_H) * 0.05;  % 5% of the max elevation
for i = 1:num_pdfs
    offset = offsets(i);
    pdf = all_pdfs{i};
    top_limit = top_limit_H(i);
    % Plot the PDF as a continuous line
    valid_indices = find(pdf > 0);  % Only plot where PDF has values
    if ~isempty(valid_indices)
        valid_H = H_cord(valid_indices); % Plot only up to the top limit
        valid_pdf = pdf(valid_indices);
        below_limit = valid_H <= top_limit; % Filter points below the top limit
        plot_H = valid_H(below_limit);
        plot_pdf = valid_pdf(below_limit);
        plot(offset-plot_pdf, plot_H, 'Color', colors(i,:), 'LineWidth', 1.5); % Draw continuous line
        x_fill = [offset - plot_pdf; repmat(offset, size(plot_pdf))];
        y_fill = [plot_H; flipud(plot_H)];
        patch(x_fill, y_fill, colors(i,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
    plot([offset, offset], [0, max(plot_H)], '--', 'Color', colors(i,:), 'LineWidth', 0.8);  
    text(offset-0.15*spacing, max(plot_H)+padding/1.5, pdf_names{i}, 'HorizontalAlignment', 'center', 'Color', colors(i,:), 'FontWeight', 'bold', 'FontSize', 20, 'FontName', 'Times');
end
ylabel('Elevation (m)', 'FontSize', 28);
set(gca, 'XTick', []);  % Remove x-axis ticks as they don't represent numerical values
grid on; box on;
ax = gca; ax.FontSize = 28;
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times')
ylim([0, max(top_limit_H) + 2*padding]);
hold off;

%% PDF Gaussian Fit
scrsz = get(0,'ScreenSize');
current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]); hold on;
colors = [
    0 0.4470 0.7410      % EW01 - blue
];
labels_base = {'EW01'};
markers = {'s'};
targetCols = [1];

[numPoints, ~] = size(CDF);
H = linspace(0, 1, numPoints)';
H_pdf = (H(1:end-1) + H(2:end)) / 2; 

legend_handles = gobjects(1, 2);
legend_labels = cell(1, 2);
hasPlottedGaussian = false;

for i = 1:length(targetCols)
    col = targetCols(i);
    F = CDF(:, col);

    f = diff(F) ./ diff(H);
    f = smooth(f, 0.1, 'loess');  % smoother
    init_params = [mean(H), std(H)];
    opts = optimset('Display','off');
    gauss_cdf = @(params, x) 0.5 * (1 + erf((x - params(1)) ./ (sqrt(2) * params(2))));
    fit_params = fminsearch(@(p) sum((gauss_cdf(p, H) - F).^2), init_params, opts);
    mu_fit = fit_params(1);
    sigma_fit = fit_params(2);
    F_fit = gauss_cdf(fit_params, H);
    f_fit = normpdf(H_pdf, mu_fit, sigma_fit);
    ks_dist = max(abs(F - F_fit));

    if strcmp(markers{i}, 'none')
        h = plot(H_pdf, f, '-', 'Color', colors(i,:), 'LineWidth', 2);
    else
        h = plot(H_pdf, f, '-', 'Color', colors(i,:), 'LineWidth', 2, ...
                 'Marker', markers{i}, 'MarkerSize', 5, ...
                 'MarkerEdgeColor', colors(i,:), ...
                 'MarkerFaceColor', colors(i,:));
    end
    legend_handles(i+1) = h;
    legend_labels{i+1} = sprintf('%s (KS=%.3f)', ...
                                 labels_base{i}, ks_dist);
    if ~hasPlottedGaussian
        h_fit = plot(H_pdf, f_fit, '--', 'Color', [0.5 0.5 0.5], ...
                     'LineWidth', 2, 'DisplayName', 'Gaussian Fit');
        legend_handles(1) = h_fit;
        legend_labels{1} = sprintf('%s (\\mu=%.2f, \\sigma=%.2f)', ...
                                 'Gaussian Fit', mu_fit, sigma_fit);
        hasPlottedGaussian = true;
    else
        plot(H_pdf, f_fit, '--', 'Color', [0.5 0.5 0.5], ...
             'LineWidth', 2, 'HandleVisibility', 'off');
    end
end
legend(legend_handles, legend_labels, 'Location', 'northeast', 'FontSize', 24);
xlabel('{\it H}', 'FontSize', 28);
ylabel("{\it f_h}({\it H}; {\it x}, {\it y})", 'FontSize', 28);
ylim([-1, 7]);
ax = gca;
ax.FontSize = 28;
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times');
grid off; box on;
