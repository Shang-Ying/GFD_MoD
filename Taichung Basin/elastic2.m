%% 2D Linear Elasticity Solver using Finite Difference Method with Direct Solver
% Script version with corrected boundary conditions based on the reference image

clear all;
close all;
clc;

%% Problem parameters
nx = 50;             % Number of grid points in x direction
ny = 50;             % Number of grid points in y direction
Lx = 1.0;            % Domain length in x direction
Ly = 1.0;            % Domain length in y direction
E = 210e9;           % Young's modulus (Pa) - steel
nu = 0.3;            % Poisson's ratio

% Grid spacing
dx = Lx / (nx - 1);
dy = Ly / (ny - 1);

% Lame parameters
mu = E / (2 * (1 + nu));
lambda = E * nu / ((1 + nu) * (1 - 2 * nu));

% Create meshgrid for plotting
[X, Y] = meshgrid(linspace(0, Lx, nx), linspace(0, Ly, ny));

% Total number of unknowns
n_dof = 2 * nx * ny;

% Build the global stiffness matrix and force vector
K = sparse(n_dof, n_dof);
F = zeros(n_dof, 1);

%% Define index mapping function
% Map (i,j) grid coordinates to linear index in the global system
% u(i,j) -> index_map(i,j)
% v(i,j) -> index_map(i,j) + nx*ny
index_map = @(i, j) (i-1)*nx + j;

%% Assemble the stiffness matrix for internal nodes
for i = 2:ny-1
    for j = 2:nx-1
        % Equation for u at node (i,j)
        row_u = index_map(i,j);
        
        % Coefficients for u
        K(row_u, index_map(i,j-1)) = K(row_u, index_map(i,j-1)) + mu/(dx^2);
        K(row_u, index_map(i-1,j)) = K(row_u, index_map(i-1,j)) + mu/(dy^2);
        K(row_u, index_map(i,j))   = K(row_u, index_map(i,j))   - 2*mu/(dx^2) - 2*mu/(dy^2) - (lambda+mu)/(dx^2) - (lambda+mu)/(dy^2);
        K(row_u, index_map(i+1,j)) = K(row_u, index_map(i+1,j)) + mu/(dy^2);
        K(row_u, index_map(i,j+1)) = K(row_u, index_map(i,j+1)) + mu/(dx^2);
        
        % Mixed derivative terms (coupling with v)
        K(row_u, index_map(i-1,j-1) + nx*ny) = K(row_u, index_map(i-1,j-1) + nx*ny) - (lambda+mu)/(4*dx*dy);
        K(row_u, index_map(i-1,j+1) + nx*ny) = K(row_u, index_map(i-1,j+1) + nx*ny) + (lambda+mu)/(4*dx*dy);
        K(row_u, index_map(i+1,j-1) + nx*ny) = K(row_u, index_map(i+1,j-1) + nx*ny) + (lambda+mu)/(4*dx*dy);
        K(row_u, index_map(i+1,j+1) + nx*ny) = K(row_u, index_map(i+1,j+1) + nx*ny) - (lambda+mu)/(4*dx*dy);
        
        % Equation for v at node (i,j)
        row_v = index_map(i,j) + nx*ny;
        
        % Coefficients for v
        K(row_v, index_map(i,j-1) + nx*ny) = K(row_v, index_map(i,j-1) + nx*ny) + mu/(dx^2);
        K(row_v, index_map(i-1,j) + nx*ny) = K(row_v, index_map(i-1,j) + nx*ny) + mu/(dy^2);
        K(row_v, index_map(i,j)   + nx*ny) = K(row_v, index_map(i,j)   + nx*ny) - 2*mu/(dx^2) - 2*mu/(dy^2) - (lambda+mu)/(dx^2) - (lambda+mu)/(dy^2);
        K(row_v, index_map(i+1,j) + nx*ny) = K(row_v, index_map(i+1,j) + nx*ny) + mu/(dy^2);
        K(row_v, index_map(i,j+1) + nx*ny) = K(row_v, index_map(i,j+1) + nx*ny) + mu/(dx^2);
        
        % Mixed derivative terms (coupling with u)
        K(row_v, index_map(i-1,j-1)) = K(row_v, index_map(i-1,j-1)) - (lambda+mu)/(4*dx*dy);
        K(row_v, index_map(i-1,j+1)) = K(row_v, index_map(i-1,j+1)) + (lambda+mu)/(4*dx*dy);
        K(row_v, index_map(i+1,j-1)) = K(row_v, index_map(i+1,j-1)) + (lambda+mu)/(4*dx*dy);
        K(row_v, index_map(i+1,j+1)) = K(row_v, index_map(i+1,j+1)) - (lambda+mu)/(4*dx*dy);
    end
end

%% Apply Dirichlet boundary conditions
bc_nodes = [];
bc_values = [];

% Left boundary (x=0): u=0, v=0 (fixed)
for i = 1:ny
    bc_nodes = [bc_nodes; index_map(i,1); index_map(i,1) + nx*ny];
    bc_values = [bc_values; 0; 0];
end

% Right boundary (x=Lx): Free boundary (no constraints)
% We don't add any constraints here

% Bottom boundary (y=0): u=0, v=0 (fixed)
for j = 1:nx
    bc_nodes = [bc_nodes; index_map(1,j); index_map(1,j) + nx*ny];
    bc_values = [bc_values; 0; 0];
end

% Top boundary (y=Ly): v=0 (fixed in y), u=prescribed displacement
u_top = 0.1*Lx;  % Prescribed displacement at top edge
for j = 1:nx
    bc_nodes = [bc_nodes; index_map(ny,j); index_map(ny,j) + nx*ny];
    bc_values = [bc_values; u_top; 0];
end

% Apply Dirichlet boundary conditions using the penalty method
penalty = 1e10;  % Penalty coefficient
for k = 1:length(bc_nodes)
    node = bc_nodes(k);
    value = bc_values(k);
    
    % Replace row with penalty approach
    K(node, :) = 0;  % Zero out the row
    K(node, node) = penalty;
    F(node) = penalty * value;
end

%% Add Neumann boundary conditions for right edge (traction-free)
% We don't need to explicitly apply Neumann BCs for the right edge
% as traction-free is the natural boundary condition of the weak form

%% Solve the system using the backslash operator
fprintf('Solving linear system with %d unknowns...\n', n_dof);
tic;
sol = K \ F;
solve_time = toc;
fprintf('System solved in %.3f seconds\n', solve_time);

% Check if solution contains NaN
if any(isnan(sol))
    fprintf('Warning: Solution contains NaN values!\n');
    
    % Try with a more numerically stable approach
    fprintf('Trying with a more stable approach...\n');
    
    % Rebuild system with a different formulation
    K = sparse(n_dof, n_dof);
    F = zeros(n_dof, 1);
    
    % Set all diagonal elements to a small value to start
    for i = 1:n_dof
        K(i,i) = 1e-10;
    end
    
    % Reassemble for interior nodes
    for i = 2:ny-1
        for j = 2:nx-1
            % Skip boundary nodes
            if any(bc_nodes == index_map(i,j)) || any(bc_nodes == index_map(i,j)+nx*ny)
                continue;
            end
            
            % Equation for u
            row_u = index_map(i,j);
            K(row_u, index_map(i,j-1)) = K(row_u, index_map(i,j-1)) + mu/(dx^2);
            K(row_u, index_map(i-1,j)) = K(row_u, index_map(i-1,j)) + mu/(dy^2);
            K(row_u, index_map(i,j))   = K(row_u, index_map(i,j))   - 2*mu/(dx^2) - 2*mu/(dy^2);
            K(row_u, index_map(i+1,j)) = K(row_u, index_map(i+1,j)) + mu/(dy^2);
            K(row_u, index_map(i,j+1)) = K(row_u, index_map(i,j+1)) + mu/(dx^2);
            
            % Additional terms for div(grad(u)) + grad(div(u))
            K(row_u, index_map(i,j-1)) = K(row_u, index_map(i,j-1)) + (lambda+mu)/(dx^2);
            K(row_u, index_map(i,j+1)) = K(row_u, index_map(i,j+1)) + (lambda+mu)/(dx^2);
            K(row_u, index_map(i,j))   = K(row_u, index_map(i,j))   - 2*(lambda+mu)/(dx^2);
            
            % Cross-coupling terms
            K(row_u, index_map(i-1,j-1) + nx*ny) = K(row_u, index_map(i-1,j-1) + nx*ny) - (lambda+mu)/(4*dx*dy);
            K(row_u, index_map(i-1,j+1) + nx*ny) = K(row_u, index_map(i-1,j+1) + nx*ny) + (lambda+mu)/(4*dx*dy);
            K(row_u, index_map(i+1,j-1) + nx*ny) = K(row_u, index_map(i+1,j-1) + nx*ny) + (lambda+mu)/(4*dx*dy);
            K(row_u, index_map(i+1,j+1) + nx*ny) = K(row_u, index_map(i+1,j+1) + nx*ny) - (lambda+mu)/(4*dx*dy);
            
            % Equation for v
            row_v = index_map(i,j) + nx*ny;
            K(row_v, index_map(i,j-1) + nx*ny) = K(row_v, index_map(i,j-1) + nx*ny) + mu/(dx^2);
            K(row_v, index_map(i-1,j) + nx*ny) = K(row_v, index_map(i-1,j) + nx*ny) + mu/(dy^2);
            K(row_v, index_map(i,j)   + nx*ny) = K(row_v, index_map(i,j)   + nx*ny) - 2*mu/(dx^2) - 2*mu/(dy^2);
            K(row_v, index_map(i+1,j) + nx*ny) = K(row_v, index_map(i+1,j) + nx*ny) + mu/(dy^2);
            K(row_v, index_map(i,j+1) + nx*ny) = K(row_v, index_map(i,j+1) + nx*ny) + mu/(dx^2);
            
            % Additional terms for div(grad(v)) + grad(div(v))
            K(row_v, index_map(i-1,j) + nx*ny) = K(row_v, index_map(i-1,j) + nx*ny) + (lambda+mu)/(dy^2);
            K(row_v, index_map(i+1,j) + nx*ny) = K(row_v, index_map(i+1,j) + nx*ny) + (lambda+mu)/(dy^2);
            K(row_v, index_map(i,j)   + nx*ny) = K(row_v, index_map(i,j)   + nx*ny) - 2*(lambda+mu)/(dy^2);
            
            % Cross-coupling terms
            K(row_v, index_map(i-1,j-1)) = K(row_v, index_map(i-1,j-1)) - (lambda+mu)/(4*dx*dy);
            K(row_v, index_map(i-1,j+1)) = K(row_v, index_map(i-1,j+1)) + (lambda+mu)/(4*dx*dy);
            K(row_v, index_map(i+1,j-1)) = K(row_v, index_map(i+1,j-1)) + (lambda+mu)/(4*dx*dy);
            K(row_v, index_map(i+1,j+1)) = K(row_v, index_map(i+1,j+1)) - (lambda+mu)/(4*dx*dy);
        end
    end
    
    % Apply boundary conditions directly
    for k = 1:length(bc_nodes)
        node = bc_nodes(k);
        value = bc_values(k);
        
        % Zero out row and set diagonal to 1
        K(node, :) = 0;
        K(node, node) = 1;
        F(node) = value;
    end
    
    % Solve again
    fprintf('Resolving system with more stable approach...\n');
    sol = K \ F;
    
    if any(isnan(sol))
        fprintf('Still getting NaN values. Trying a different approach...\n');
        
        % Try solving with regularization
        K_reg = K + speye(size(K)) * 1e-6;
        sol = K_reg \ F;
        
        if any(isnan(sol))
            error('Failed to obtain a valid solution. Please check your boundary conditions and equations.');
        end
    end
end

%% Extract the solution
u = reshape(sol(1:nx*ny), [ny, nx]);
v = reshape(sol(nx*ny+1:end), [ny, nx]);

%% Calculate displacement magnitude and strain
displacement = sqrt(u.^2 + v.^2);

% Calculate deformed mesh
X_def = X + u;
Y_def = Y + v;

% Strain components
[du_dx, du_dy] = gradient(u, Lx/(size(u,2)-1), Ly/(size(u,1)-1));
[dv_dx, dv_dy] = gradient(v, Lx/(size(v,2)-1), Ly/(size(v,1)-1));

% Calculate von Mises strain
e_xx = du_dx;
e_yy = dv_dy;
e_xy = 0.5 * (du_dy + dv_dx);
von_mises_strain = sqrt(e_xx.^2 + e_yy.^2 + 2*e_xy.^2 - e_xx.*e_yy);

%% Create plots
figure('Position', [100, 100, 1200, 800]);

% Plot displacement magnitude
subplot(2, 2, 1);
surf(X, Y, displacement);
title('Displacement Magnitude');
xlabel('x'); ylabel('y'); zlabel('|u|');
colorbar;

% Plot deformed mesh
subplot(2, 2, 2);
mesh(X_def, Y_def, zeros(size(X_def)));
hold on;
mesh(X, Y, -0.1*ones(size(X)), 'EdgeColor', [0.7 0.7 0.7]);
title('Deformed Mesh');
xlabel('x'); ylabel('y');
view(2); % 2D view
axis equal;

% Plot displacement in x-direction
subplot(2, 2, 3);
surf(X, Y, u);
title('x-displacement (u)');
xlabel('x'); ylabel('y');
colorbar;

% Plot von Mises strain
subplot(2, 2, 4);
surf(X, Y, von_mises_strain);
title('von Mises Strain');
xlabel('x'); ylabel('y');
colorbar;

% Create a new figure for the quiver plot
figure;
quiver(X, Y, u, v, 2);
title('Displacement Vector Field');
xlabel('x'); ylabel('y');
axis equal;

% Print maximum displacement
fprintf('Maximum displacement: %.6e\n', max(displacement(:)));