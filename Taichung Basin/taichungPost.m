% Target F_h value and tolerance
F_h_target = 0.15;
tolerance = 0.01;

% Filter points where F_h is approximately 0.85
target_points = abs(Fh - F_h_target) <= tolerance;

% Get coordinates of these points
x_target = tpH(target_points,1)*y_range+x_min;
y_target = tpH(target_points,2)*y_range+y_min;
z_target = tpH(target_points,3)*u_range+min(d_est);

% Create the shrunk boundary for filtering interior points
boundary_shrink_factor = 0.97; % More restrictive to avoid edge artifacts
center_x = mean(bdy_dimless(:,1)*y_range+x_min);
center_y = mean(bdy_dimless(:,2)*y_range+y_min);

boundary_adj = zeros(size(bdy_dimless));
for i = 1:size(bdy_dimless,1)
    boundary_adj(i,1) = center_x + boundary_shrink_factor * (bdy_dimless(i,1)*y_range+x_min - center_x);
    boundary_adj(i,2) = center_y + boundary_shrink_factor * (bdy_dimless(i,2)*y_range+y_min - center_y);
end

% Find points inside the shrunk boundary
is_inside = inpolygon(x_target, y_target, boundary_adj(:,1), boundary_adj(:,2));

% Keep only interior points
x_filtered = x_target(is_inside);
y_filtered = y_target(is_inside);
z_filtered = z_target(is_inside);

% Identify and smooth points near the boundary
dist_to_boundary = zeros(size(x_filtered));
for i = 1:length(x_filtered)
    [min_dist, ~] = min(sqrt((x_filtered(i) - boundary_adj(:,1)).^2 + ...
                           (y_filtered(i) - boundary_adj(:,2)).^2));
    dist_to_boundary(i) = min_dist;
end

% Define edge points
edge_threshold = 500; % Distance in meters
is_edge = dist_to_boundary < edge_threshold;

% Apply smoothing to edge points
if any(is_edge)
    edge_neighbors = 10; % Number of neighbors to consider
    for i = find(is_edge)'
        % Find nearest neighbors
        dists = sqrt((x_filtered(i) - x_filtered).^2 + (y_filtered(i) - y_filtered).^2);
        [~, idx] = sort(dists);
        neighbors = idx(2:min(edge_neighbors+1, length(idx))); % Skip first (self)
        % Replace with average
        z_filtered(i) = mean(z_filtered(neighbors));
    end
end

% Setup figure
scrsz = get(0,'ScreenSize');
current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);

% Create a grid strictly within the boundary
[X_grid, Y_grid] = meshgrid(linspace(min(x_filtered), max(x_filtered), 200), ...
                           linspace(min(y_filtered), max(y_filtered), 200));

% Create a mask to exclude grid points outside the boundary
grid_inside = inpolygon(X_grid(:), Y_grid(:), boundary_adj(:,1), boundary_adj(:,2));
grid_inside = reshape(grid_inside, size(X_grid));

% Interpolate only within the boundary
Z_grid = griddata(x_filtered, y_filtered, z_filtered, X_grid, Y_grid, 'linear');

% Apply the mask - set values outside boundary to NaN
Z_grid(~grid_inside) = NaN;

% Create surface with masked grid
surface = surf(X_grid, Y_grid, Z_grid);
set(surface, 'EdgeColor', 'none');

% Apply colormap for elevation
fullColormap = othercolor('GnBu7');
colormap(fullColormap);
c = colorbar;
ylabel(c, 'Elevation (m)', 'FontSize', 20, 'FontName', 'Times');

% Set color axis limits
caxis([min(z_filtered) max(z_filtered)]);

% View settings
view_az = 27;
view_el = 60;
view(view_az, view_el);
grid on;
axis equal;
daspect auto;

% Adjust axis position to leave space for colorbar
ax = gca;
pos = ax.Position;
pos(3) = pos(3) * 0.82; % Shrink width to make room for colorbar
ax.Position = pos;

% Add colorbar outside plot
chb = colorbar('eastoutside');
ylabel(chb, 'Hydraulic head (m)', ...
       'FontSize', 24, 'FontName', 'Times');
chb.Position = [0.85, 0.1, 0.02, 0.8]; % Adjust colorbar position manually

% % Set color limits and ticks
% clim([0 1]); % Assuming normalized Fh values
% Tspc = linspace(0,1,6);
% set(chb,'Ticks',Tspc);
% TL = arrayfun(@(x) sprintf('%.1f',x),Tspc,'un',0);
% set(chb,'TickLabels',TL);

% % Add annotation for the surface
% text_x = max(x_filtered);  % 右側
% text_y = max(y_filtered);  % 上側
% text_z = max(z_filtered) + 50;  % 提高一點以免被曲面擋住
% text(text_x, text_y, text_z, 'Isosurface {\it F}_{\ith} = 0.85', ...
%     'FontSize', 20, 'FontName', 'Times', 'HorizontalAlignment', 'right');

% Set font and formatting
xtickformat('%.0f');
ytickformat('%.0f');
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;

% Tick settings
tick_min = ceil(min(ylim)/5000) * 5000;
tick_max = floor(max(ylim)/5000) * 5000;
yticks(tick_min:5000:tick_max);
tick_min = ceil(min(zlim)/100) * 100;
tick_max = floor(max(zlim)/100) * 100;
zticks(tick_min:100:tick_max);

xlabel('{\it x} (m)', 'FontSize', 24);
ylabel('{\it y} (m)', 'FontSize', 24);
zlabel('{\it H} (m)', 'FontSize', 24);
ax.FontSize = max(20, min(14, 0.02 * ax.Position(3)*ax.Position(4)));

% Add boundary overlay
hold on;
plot3(bdy_dimless(:,1)*y_range+x_min, bdy_dimless(:,2)*y_range+y_min, ...
     ones(size(bdy_dimless,1),1)*min(zlim), 'Color', [0.5 0.5 0.5], ...
     'LineWidth', 1.5, 'LineStyle', ':');

% Apply font settings
set(findall(gcf,'type','text'), 'FontName', 'Times');
set(gca, 'FontName', 'Times');

% Enhance 3D perception
lighting gouraud;
material dull;
camlight('headlight');