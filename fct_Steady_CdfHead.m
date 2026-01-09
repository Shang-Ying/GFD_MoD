% function [u] = fct_Steady_CdfHead(Wcoef_3D,idx,ns,tpH,alpha,beta,g,Y_mean,u_mean)
% % GE: K(4+7) + [alpha(x,y)(H-h(x,y))+beta(x,y)+g(x,y)](3) = 0
% 
% C = sparse((size(tpH,1)),(size(tpH,1)));
% 
% G_In = (tpH(:,4) == 0);
% D_Bd = (tpH(:,4) == 1);
% N_Bd = (tpH(:,4) == 2);
% 
% G_loop = find(G_In);
% D_loop = find(D_Bd);
% N_loop = find(N_Bd);
% 
% for ith = 1:length(G_loop)
%     % Cia = (Wcoef_3D(4,1:ns,G_loop(ith)) + Wcoef_3D(7,1:ns,G_loop(ith))) * exp(Y_mean(G_loop(ith))); % K(4+7)
%     Cia = (Wcoef_3D(4,1:ns,G_loop(ith)) + Wcoef_3D(7,1:ns,G_loop(ith)) + Wcoef_3D(9,1:ns,G_loop(ith))) * exp(Y_mean(G_loop(ith))); % K(4+7+9)
%     Cib = alpha(G_loop(ith))*(tpH(G_loop(ith),3) - u_mean(G_loop(ith))); % alpha(x,y)(H-h(x,y))
%     Cic = beta(G_loop(ith)) + g(G_loop(ith)); % beta(x,y)+g(x,y) nearly 0
%     Cid = (Cib+Cic) * Wcoef_3D(3,1:ns,G_loop(ith)); % (Cib + Cic)(3)
%     GE = Cia + Cid;
%     C(G_loop(ith),(idx(G_loop(ith),1:ns))) = GE;
%     f(G_loop(ith),1) = 0; % No Sink/Source
% end
% 
% for ith = 1:length(D_loop)
%     C(D_loop(ith),D_loop(ith)) = 1;
%     f(D_loop(ith),1) = tpH(D_loop(ith),5);
% end
% 
% for ith = 1:length(N_loop)
%     C(N_loop(ith),(idx(N_loop(ith),1:ns))) = tpH(N_loop(ith),6).*Wcoef_3D(1,1:ns,N_loop(ith)) + tpH(N_loop(ith),7).*Wcoef_3D(2,1:ns,N_loop(ith)) + tpH(N_loop(ith),8).*Wcoef_3D(3,1:ns,N_loop(ith));
%     f(N_loop(ith),1) = 0;
% end
% 
% u=C\f;
% 
% function [u] = fct_Steady_CdfHead(Wcoef_3D, idx, ns, tpH, alpha, beta, g, Y_mean, u_mean)
% % GE: K(4+7+9) +  = 0
% 
% % Get the total number of nodes
% nNodes = size(tpH, 1);
% 
% % Preallocate vectors for sparse matrix construction
% max_entries = nnz(tpH(:,4) == 0) * ns + nnz(tpH(:,4) == 1) + nnz(tpH(:,4) == 2) * ns;
% rows = zeros(max_entries, 1);
% cols = zeros(max_entries, 1);
% values = zeros(max_entries, 1);
% 
% % Preallocate right-hand side vector
% f = zeros(nNodes, 1);
% 
% % Initialize counter
% counter = 1;
% 
% % Identify nodes based on their type
% G_In = (tpH(:, 4) == 0); % Interior nodes
% D_Bd = (tpH(:, 4) == 1); % Dirichlet boundary nodes
% N_Bd = (tpH(:, 4) == 2); % Neumann boundary nodes
% 
% % Loop over interior nodes (G)
% G_loop = find(G_In);
% for ith = 1:length(G_loop)
%     node = G_loop(ith);
% 
%     % Compute coefficients
%     Cia = (Wcoef_3D(4, 1:ns, node) + Wcoef_3D(7, 1:ns, node) + Wcoef_3D(9, 1:ns, node)) * exp(Y_mean(node)); % K(4+7+9)
%     Cib = alpha(node) * (tpH(node, 3) - u_mean(node)); % alpha(x,y)(H-h(x,y))
%     Cic = beta(node) + g(node); % beta(x,y) + g(x,y)
%     Cid = (Cib + Cic) * Wcoef_3D(3, 1:ns, node); % (Cib + Cic)(3)
%     GE = Cia + Cid;
% 
%     % Populate sparse matrix triplets
%     rows(counter:counter + ns - 1) = node;
%     cols(counter:counter + ns - 1) = idx(node, 1:ns);
%     values(counter:counter + ns - 1) = GE;
%     counter = counter + ns;
% 
%     % Set the RHS to zero (no sink/source)
%     % f(node, 1) = 0;
% end
% 
% % Loop over Dirichlet boundary nodes (D)
% D_loop = find(D_Bd);
% for ith = 1:length(D_loop)
%     node = D_loop(ith);
% 
%     % Add Dirichlet condition to the sparse matrix
%     rows(counter) = node;
%     cols(counter) = node;
%     values(counter) = 1;
%     counter = counter + 1;
% 
%     % Set RHS to the prescribed value
%     f(node, 1) = tpH(node, 5);
% end
% 
% % Loop over Neumann boundary nodes (N)
% N_loop = find(N_Bd);
% for ith = 1:length(N_loop)
%     node = N_loop(ith);
% 
%     % Compute coefficients for Neumann condition
%     GE_N = tpH(node, 6) .* Wcoef_3D(1, 1:ns, node) + ...
%            tpH(node, 7) .* Wcoef_3D(2, 1:ns, node) + ...
%            tpH(node, 8) .* Wcoef_3D(3, 1:ns, node);
% 
%     % Populate sparse matrix triplets
%     rows(counter:counter + ns - 1) = node;
%     cols(counter:counter + ns - 1) = idx(node, 1:ns);
%     values(counter:counter + ns - 1) = GE_N;
%     counter = counter + ns;
% 
%     % Set RHS to zero
%     % f(node, 1) = 0;
% end
% 
% % Assemble the sparse matrix
% C = sparse(rows(1:counter - 1), cols(1:counter - 1), values(1:counter - 1), nNodes, nNodes);
% 
% % Solve the linear system
% u = C \ f;
% end

function [u] = fct_Steady_CdfHead(Wcoef_3D, idx, ns, tpH, alpha, beta, g, Y_mean_CDF, u_mean_CDF)
% Solve a steady-state CDF head problem with parallel computing.

% Get the total number of nodes
nNodes = size(tpH, 1);

% Preallocate cell arrays for parallelization
numInteriorNodes = nnz(tpH(:, 4) == 0); % Number of interior nodes
numDirichletNodes = nnz(tpH(:, 4) == 1); % Number of Dirichlet nodes
numNeumannNodes = nnz(tpH(:, 4) == 2); % Number of Neumann nodes

rowsInterior = cell(numInteriorNodes, 1);
colsInterior = cell(numInteriorNodes, 1);
valuesInterior = cell(numInteriorNodes, 1);

rowsDirichlet = cell(numDirichletNodes, 1);
colsDirichlet = cell(numDirichletNodes, 1);
valuesDirichlet = cell(numDirichletNodes, 1);
rhsDirichlet = cell(numDirichletNodes, 1);

rowsNeumann = cell(numNeumannNodes, 1);
colsNeumann = cell(numNeumannNodes, 1);
valuesNeumann = cell(numNeumannNodes, 1);

% Preallocate right-hand side vector
f = zeros(nNodes, 1);

% Identify nodes based on their type
G_loop = find(tpH(:, 4) == 0); % Interior nodes
D_loop = find(tpH(:, 4) == 1); % Dirichlet boundary nodes
N_loop = find(tpH(:, 4) == 2); % Neumann boundary nodes

% Parallel loop for interior nodes
parfor ith = 1:length(G_loop)
    node = G_loop(ith);

    % Compute coefficients
    Cia = (Wcoef_3D(4, 1:ns, node) + Wcoef_3D(7, 1:ns, node) + Wcoef_3D(9, 1:ns, node)) * exp(Y_mean_CDF(node));
    Cib = alpha(node) * (tpH(node, 3) - u_mean_CDF(node));
    Cic = beta(node) + g(node);
    Cid = (Cib + Cic) * Wcoef_3D(3, 1:ns, node);
    GE = Cia + Cid;

    % Store results
    rowsInterior{ith} = repmat(node, ns, 1);
    colsInterior{ith} = idx(node, 1:ns)';
    valuesInterior{ith} = GE';
end

% Parallel loop for Dirichlet boundary nodes
parfor ith = 1:length(D_loop)
    node = D_loop(ith);

    % Store results
    rowsDirichlet{ith} = node;
    colsDirichlet{ith} = node;
    valuesDirichlet{ith} = 1;
    rhsDirichlet{ith} = tpH(node, 5);
end

% Parallel loop for Neumann boundary nodes
parfor ith = 1:length(N_loop)
    node = N_loop(ith);

    % Compute coefficients for Neumann condition
    GE_N = tpH(node, 6) .* Wcoef_3D(1, 1:ns, node) + ...
           tpH(node, 7) .* Wcoef_3D(2, 1:ns, node) + ...
           tpH(node, 8) .* Wcoef_3D(3, 1:ns, node);

    % Store results
    rowsNeumann{ith} = repmat(node, ns, 1);
    colsNeumann{ith} = idx(node, 1:ns)';
    valuesNeumann{ith} = GE_N';
end

% Combine results from all parallel loops
rows = [vertcat(rowsInterior{:}); vertcat(rowsDirichlet{:}); vertcat(rowsNeumann{:})];
cols = [vertcat(colsInterior{:}); vertcat(colsDirichlet{:}); vertcat(colsNeumann{:})];
values = [vertcat(valuesInterior{:}); vertcat(valuesDirichlet{:}); vertcat(valuesNeumann{:})];

% Combine right-hand side contributions
f(vertcat(D_loop)) = vertcat(rhsDirichlet{:});

% Assemble the sparse matrix
C = sparse(rows, cols, values, nNodes, nNodes);

% Solve the linear system
u = C \ f;
% u = mldivide(C,f);
end
