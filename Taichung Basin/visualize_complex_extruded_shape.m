function visualize_complex_extruded_shape(xy_nodes, BC_type, z_min, z_max)
    % Set figure size to match screen height
    scrsz = get(0, 'ScreenSize');
    fig = figure('OuterPosition', [0 0 (scrsz(4)+100) scrsz(4)]);
    hold on;
    
    % Set font to Times
    set(fig, 'DefaultTextFontName', 'Times');
    set(fig, 'DefaultAxesFontName', 'Times');
    set(fig, 'DefaultAxesFontSize', 28);
    set(fig, 'DefaultTextFontSize', 28);
    
    % Get number of points in the 2D shape
    n_points = size(xy_nodes, 1);
    
    % Create bottom and top vertices
    bottom_vertices = [xy_nodes, ones(n_points, 1) * z_min];
    top_vertices = [xy_nodes, ones(n_points, 1) * z_max];
    
    % Define boundary condition colors using the new hex colors from the image
    % Convert from hex to RGB (divide by 255)
    bc_colors = {
        [171/255, 224/255, 240/255],   % Light Blue #ABE0F0 (ABE0F0)
        [197/255, 219/255, 107/255],   % Green #C5DB6B (C5DB6B)      
        [255/255, 209/255, 116/255],   % Yellow #FFD174 (FFD174)
        [255/255, 195/255, 195/255],   % Pink #FFC3C3 (FFC3C3)
        [0.5, 0.5, 0.5],               % Gray (Top face)
        [0.3, 0.3, 0.3]                % Dark Gray (Bottom face)
    };
    
    % Set transparency
    alpha_value = 0.6;
    
    % Find transitions between Dirichlet and Neumann
    % This will identify the segments
    transitions = find(diff([BC_type; BC_type(1)]) ~= 0);
    
    % If no transitions (all same type), create artificial division
    if isempty(transitions)
        if BC_type(1) == 1  % All Dirichlet
            segment_D = {1:n_points};
            segment_N = {};
        else  % All Neumann
            segment_D = {};
            segment_N = {1:n_points};
        end
    else
        % Initialize segment containers
        segment_D = {};
        segment_N = {};
        
        % Process each segment
        for i = 1:length(transitions)
            start_idx = transitions(i) + 1;
            if i < length(transitions)
                end_idx = transitions(i+1);
            else
                end_idx = transitions(1) + n_points;
            end
            
            % Handle wrap-around for the last segment
            if end_idx > n_points
                indices = [start_idx:n_points, 1:(end_idx-n_points)];
            else
                indices = start_idx:end_idx;
            end
            
            % Add to appropriate segment list based on boundary type
            if BC_type(start_idx) == 1  % Dirichlet
                segment_D{end+1} = indices;
            else  % Neumann
                segment_N{end+1} = indices;
            end
        end
    end
    
    % Draw bottom face - with no edge lines
    patch('Vertices', bottom_vertices, 'Faces', 1:n_points, ...
          'FaceColor', bc_colors{6}, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    
    % Draw top face - with no edge lines
    patch('Vertices', top_vertices, 'Faces', 1:n_points, ...
          'FaceColor', bc_colors{5}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    % Draw Dirichlet boundary segments
    h_legend = [];
    legend_labels = {};
    
    for s = 1:length(segment_D)
        h_temp = draw_segment(segment_D{s}, bottom_vertices, top_vertices, bc_colors{s}, alpha_value, n_points);
        h_legend(end+1) = h_temp;
        % Correct Greek letter notation with proper subscripts (italic Gamma, smaller subscripts)
        legend_labels{end+1} = ['$\it{\Gamma}_{D_{\scriptscriptstyle ' num2str(s) '}}$'];
    end
    
    % Draw Neumann boundary segments
    for s = 1:length(segment_N)
        h_temp = draw_segment(segment_N{s}, bottom_vertices, top_vertices, bc_colors{2+s}, alpha_value, n_points);
        h_legend(end+1) = h_temp;
        % Correct Greek letter notation with proper subscripts (italic Gamma, smaller subscripts)
        legend_labels{end+1} = ['$\it{\Gamma}_{N_{\scriptscriptstyle ' num2str(s) '}}$'];
    end
    

    
    % Set labels with units in parentheses using \rm for roman (non-italic) font
    ax = gca;
    xtickformat('%.0f') % Force whole numbers
    ytickformat('%.0f') % Force whole numbers
    ax.XAxis.Exponent = 0;
    ax.YAxis.Exponent = 0;
    % ax.YTickLabelRotation = 60;
    tick_min = ceil(min(ylim)/5000) * 5000; % Round up to the nearest 5000
    tick_max = floor(max(ylim)/5000) * 5000; % Round down to the nearest 5000
    yticks(tick_min:5000:tick_max); % Define Y-ticks with 5000 interval
    ax.FontSize = 20;
    xlabel('{\it x} (m)', 'FontSize', 24);
    ylabel('{\it y} (m)', 'FontSize', 24);
    zlabel('{\it H} (m)', 'FontSize', 24);
    % xlabel('$\rm{x}$ $\rm{(m)}$', 'FontSize', 28, 'Interpreter', 'latex');
    % ylabel('$\rm{y}$ $\rm{(m)}$', 'FontSize', 28, 'Interpreter', 'latex');
    % zlabel('$\rm{H}$ $\rm{(m)}$', 'FontSize', 28, 'Interpreter', 'latex');

    % Add top and bottom to legend with correct notation
    h_legend(end+1) = patch(NaN, NaN, bc_colors{5}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    legend_labels{end+1} = '$\rm\it{\Gamma}_{top}$';

    h_legend(end+1) = patch(NaN, NaN, bc_colors{6}, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    legend_labels{end+1} = '$\it{\Gamma}_{bottom}$';  

    % Add legend with adjusted font size and latex interpreter for proper Greek symbols
    % legend(h_legend, legend_labels, 'Location', 'northwest', 'FontSize', 28, 'Interpreter', 'latex');
    
    % Set the new viewing angles as specified
    view_az = 17;    % New azimuth angle
    view_el = 80;    % New elevation angle
    view(view_az, view_el);
    
    % Adjust axis settings
    grid on;
    axis equal;
    daspect auto
    
    % Force font to Times for all text elements
    all_text = findall(fig, 'Type', 'Text');
    for i = 1:length(all_text)
        set(all_text(i), 'FontName', 'Times');
    end
    
    all_axes = findall(fig, 'Type', 'Axes');
    for i = 1:length(all_axes)
        set(all_axes(i), 'FontName', 'Times');
    end
    
    hold off;
end

function h = draw_segment(segment, bottom_vertices, top_vertices, color, alpha_value, n_points)
    % Helper function to draw a segment of the boundary
    % Returns a handle for the legend
    
    for i = 1:length(segment)
        % Handle the case where we need to connect back to the beginning
        if i == length(segment)
            if segment(i) == n_points
                j_idx = 1;  % Connect back to first point
            else
                j_idx = segment(i) + 1;  % Connect to next point
            end
        else
            j_idx = segment(i+1);
        end
        
        i_idx = segment(i);
        
        % Create the four corners of this side face
        face_vertices = [
            bottom_vertices(i_idx,:);
            bottom_vertices(j_idx,:);
            top_vertices(j_idx,:);
            top_vertices(i_idx,:)
        ];
        
        % Create side face with NO edge lines
        h = patch('Vertices', face_vertices, 'Faces', [1 2 3 4], ...
              'FaceColor', color, 'FaceAlpha', alpha_value, 'EdgeColor', 'none');
    end
end

% Function to load data from file and visualize
function visualize_from_file(coord_file, bc_file, z_min, z_max)
    % Load coordinates
    coords = load(coord_file);
    xy_nodes = coords(:, 1:2);
    
    % Load boundary condition types
    BC_type = load(bc_file);
    
    % Call visualization function
    visualize_complex_extruded_shape(xy_nodes, BC_type, z_min, z_max);
end

% Function to save the visualization as a high-quality image
function save_visualization(filename, fig_handle, resolution)
    if nargin < 3
        resolution = 600; % Default DPI
    end
    
    if nargin < 2
        fig_handle = gcf; % Get current figure
    end
    
    % Save the figure
    print(fig_handle, filename, '-dpng', ['-r' num2str(resolution)]);
    
    fprintf('Figure saved as %s with resolution %d DPI\n', filename, resolution);
end