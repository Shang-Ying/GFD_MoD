function [u] = fct_Steady_CdfHead_FDM(tpH,alpha_CDF,beta_CDF,g_CDF,Y_mean_CDF,u_mean_CDF,dH,Sten,xy_res)
% % GE: K(4+7) + [alpha(x,y)(H-h(x,y))+beta(x,y)+g(x,y)](3) = 0
% 
% C = zeros((size(tpH,1)),(size(tpH,1)));
% 
% G_loop = Sten.G_loop;
% D_loop = Sten.D_loop;
% Nxm_loop = Sten.Nxm_loop;
% NxM_loop = Sten.NxM_loop;
% Nym_loop = Sten.Nym_loop;
% NyM_loop = Sten.NyM_loop;
% Nzm_loop = Sten.Nzm_loop;
% NzM_loop = Sten.NzM_loop;
% % n_right = Sten.n_right;
% % n_left = Sten.n_left;
% % n_front = Sten.n_front;
% % n_back = Sten.n_back;
% % n_upper = Sten.n_upper;
% % n_lower = Sten.n_lower;
% 
% d = tpH(1,2)-tpH(2,2); % x & y
% % h = waitbar(0, 'Interior Nodes Processing...'); % Create a progress bar
% for ith = 1:length(G_loop)
%     C(G_loop(ith),(G_loop(ith)+xy_res)) = exp(Y_mean_CDF(G_loop(ith)))/(d^2);
%     C(G_loop(ith),(G_loop(ith)-xy_res)) = exp(Y_mean_CDF(G_loop(ith)))/(d^2);
%     C(G_loop(ith),(G_loop(ith)+1)) = exp(Y_mean_CDF(G_loop(ith)))/(d^2);
%     C(G_loop(ith),(G_loop(ith)-1)) = exp(Y_mean_CDF(G_loop(ith)))/(d^2);
%     Cib = (alpha_CDF(G_loop(ith))*(tpH(G_loop(ith),3) - u_mean_CDF(G_loop(ith)))); % alpha(x,y)(H-h(x,y))
%     Cic = (beta_CDF(G_loop(ith)) + g_CDF(G_loop(ith))); % beta(x,y)+g(x,y) nearly 0
%     C(G_loop(ith),(G_loop(ith)+xy_res*xy_res)) = exp(Y_mean_CDF(G_loop(ith)))/(dH^2) + (Cib+Cic)/(2*dH);
%     C(G_loop(ith),(G_loop(ith)-xy_res*xy_res)) = exp(Y_mean_CDF(G_loop(ith)))/(dH^2) - (Cib+Cic)/(2*dH);
%     C(G_loop(ith),G_loop(ith)) = (-4*exp(Y_mean_CDF(G_loop(ith)))/(d^2)) -(2*exp(Y_mean_CDF(G_loop(ith)))/(dH^2)); 
%     f(G_loop(ith),1) = 0; % No Sink/Source
%     % waitbar(ith/length(G_loop), h, ['Interior Nodes Progress: ', num2str(round(ith/length(G_loop)*100)), '%']);
% end
% % close(h); % Close the progress bar
% 
% for ith = 1:length(D_loop)
%     C(D_loop(ith),D_loop(ith)) = 1;
%     f(D_loop(ith),1) = tpH(D_loop(ith),5);
% end
% 
% for ith = 1:length(Nxm_loop)
%         C(Nxm_loop(ith),(Nxm_loop(ith)+xy_res)) = 1/d;
%         C(Nxm_loop(ith),Nxm_loop(ith)) = -1/d;
%         f(Nxm_loop(ith),1) = 0;
% end
% 
% for ith = 1:length(NxM_loop)
%         C(NxM_loop(ith),(NxM_loop(ith)-xy_res)) = 1/d;
%         C(NxM_loop(ith),NxM_loop(ith)) = -1/d;
%         f(NxM_loop(ith),1) = 0;
% end
% 
% for ith = 1:length(Nym_loop)
%         C(Nym_loop(ith),(Nym_loop(ith)+1)) = 1/d;
%         C(Nym_loop(ith),Nym_loop(ith)) = -1/d;
%         f(Nym_loop(ith),1) = 0;
% end
% 
% for ith = 1:length(NyM_loop)
%         C(NyM_loop(ith),(NyM_loop(ith)-1)) = 1/d;
%         C(NyM_loop(ith),NyM_loop(ith)) = -1/d;
%         f(NyM_loop(ith),1) = 0;
% end
% 
% for ith = 1:length(Nzm_loop)
%         C(Nzm_loop(ith),(Nzm_loop(ith)+xy_res*xy_res)) = 1/dH;
%         C(Nzm_loop(ith),Nzm_loop(ith)) = -1/dH;
%         f(Nzm_loop(ith),1) = 0;
% end
% 
% for ith = 1:length(NzM_loop)
%         C(NzM_loop(ith),(NzM_loop(ith)-xy_res*xy_res)) = 1/dH;
%         C(NzM_loop(ith),NzM_loop(ith)) = -1/dH;
%         f(NzM_loop(ith),1) = 0;
% end
% 
% u=C\f;

% Grid parameters
n = size(tpH, 1);
d = tpH(1, 2) - tpH(2, 2); % Spatial discretization (x, y)

% Preallocate storage for sparse matrix and RHS vector
max_entries = 15 * n; % Estimate maximum number of non-zero entries
row = zeros(max_entries, 1);
col = zeros(max_entries, 1);
val = zeros(max_entries, 1);
f = zeros(n, 1);

% Initialize index for sparse matrix entry tracking
entry_idx = 1;

% Extract stencil groups
G_loop = Sten.G_loop;
D_loop = Sten.D_loop;
Nxm_loop = Sten.Nxm_loop;
NxM_loop = Sten.NxM_loop;
Nym_loop = Sten.Nym_loop;
NyM_loop = Sten.NyM_loop;
Nzm_loop = Sten.Nzm_loop;
NzM_loop = Sten.NzM_loop;

% Interior nodes (G_loop)
for ith = 1:length(G_loop)
    idx = G_loop(ith);
    exp_Y = exp(Y_mean_CDF(idx));
    Cib = alpha_CDF(idx) * (tpH(idx, 3) - u_mean_CDF(idx));
    Cic = beta_CDF(idx) + g_CDF(idx);
    
    % Add contributions to sparse matrix
    row(entry_idx:entry_idx+6) = idx;
    col(entry_idx:entry_idx+6) = [idx+xy_res, idx-xy_res, idx+1, idx-1, idx+xy_res*xy_res, idx-xy_res*xy_res, idx];
    val(entry_idx:entry_idx+6) = [exp_Y/d^2, exp_Y/d^2, exp_Y/d^2, exp_Y/d^2, ...
                                  exp_Y/dH^2 + (Cib+Cic)/(2*dH), exp_Y/dH^2 - (Cib+Cic)/(2*dH), ...
                                  -4*exp_Y/d^2 - 2*exp_Y/dH^2];
    entry_idx = entry_idx + 7;
    
    % Right-hand side
    f(idx) = 0; % No Sink/Source
end

% Dirichlet boundary nodes (D_loop)
for ith = 1:length(D_loop)
    idx = D_loop(ith);
    row(entry_idx) = idx;
    col(entry_idx) = idx;
    val(entry_idx) = 1;
    entry_idx = entry_idx + 1;
    
    % Right-hand side
    f(idx) = tpH(idx, 5);
end

% Boundary conditions (Neumann loops)
neumann_loops = {Nxm_loop, NxM_loop, Nym_loop, NyM_loop, Nzm_loop, NzM_loop};
offsets = [+xy_res, -xy_res, +1, -1, +xy_res*xy_res, -xy_res*xy_res];
coeffs = [1/d, -1/d, 1/d, -1/d, 1/dH, -1/dH];

for i = 1:length(neumann_loops)
    loop = neumann_loops{i};
    offset = offsets(i);
    coeff = coeffs(i);
    for ith = 1:length(loop)
        idx = loop(ith);
        row(entry_idx:entry_idx+1) = [idx, idx];
        col(entry_idx:entry_idx+1) = [idx+offset, idx];
        val(entry_idx:entry_idx+1) = [coeff, -coeff];
        entry_idx = entry_idx + 2;
        
        % Right-hand side
        f(idx) = 0;
    end
end

% Assemble sparse matrix
C = sparse(row(1:entry_idx-1), col(1:entry_idx-1), val(1:entry_idx-1), n, n);

% Solve system
u = C \ f;

end