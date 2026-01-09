% If you don't have the 'hatch' function, this simple version can be added to your script
% (it's a simplified version of the full hatch function available on MATLAB File Exchange)
function hatch(h, style, angle, spacing, color, background, width)
    % Simple hatching function
    % Get patch vertices and faces
    v = get(h, 'Vertices');
    f = get(h, 'Faces');
    
    % Get the patch boundaries
    xmin = min(v(:,1));
    xmax = max(v(:,1));
    ymin = min(v(:,2));
    ymax = max(v(:,2));
    
    % Calculate number of hatch lines needed
    width_patch = max(xmax-xmin, ymax-ymin);
    num_lines = ceil(width_patch / spacing) * 2;
    
    % Create hatch lines
    theta = angle * pi/180;
    center_x = (xmin + xmax)/2;
    center_y = (ymin + ymax)/2;
    
    % Start at a distance from the center
    r = width_patch;
    
    for i = -num_lines:num_lines
        % Create a line passing through the center
        x1 = center_x - r*cos(theta) + i*spacing*sin(theta);
        y1 = center_y - r*sin(theta) - i*spacing*cos(theta);
        x2 = center_x + r*cos(theta) + i*spacing*sin(theta);
        y2 = center_y + r*sin(theta) - i*spacing*cos(theta);
        
        % Create the line
        line([x1 x2], [y1 y2], 'Color', color, 'LineWidth', width, 'Clipping', 'on');
    end
end