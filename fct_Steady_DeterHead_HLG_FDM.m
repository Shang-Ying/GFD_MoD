function [u] = fct_Steady_DeterHead_HLG_FDM(tp_Y,tp_tvV,G_loop,D_loop,Nxm_loop,NxM_loop,Nym_loop,NyM_loop,n_right,n_left,n_upper,n_lower)
% C = sparse((size(tp_tvV,1)),(size(tp_tvV,1)));
% 
% d = tp_tvV(1,2)-tp_tvV(2,2);
% K = exp(tp_Y(:,3));
% 
% for ith = 1:length(G_loop)
%     Kavg_right = geomean([K(G_loop(ith)) K(n_right(G_loop(ith)))]);
%     Kavg_left = geomean([K(G_loop(ith)) K(n_left(G_loop(ith)))]);
%     Kavg_upper = geomean([K(G_loop(ith)) K(n_upper(G_loop(ith)))]);
%     Kavg_lower = geomean([K(G_loop(ith)) K(n_lower(G_loop(ith)))]);
%     C(G_loop(ith),n_right(G_loop(ith))) = Kavg_right/(d^2);
%     C(G_loop(ith),n_left(G_loop(ith))) = Kavg_left/(d^2);
%     C(G_loop(ith),n_upper(G_loop(ith))) = Kavg_upper/(d^2);
%     C(G_loop(ith),n_lower(G_loop(ith))) = Kavg_lower/(d^2);
%     C(G_loop(ith),G_loop(ith)) = -(Kavg_right + Kavg_left + Kavg_upper + Kavg_lower)/(d^2);
%     f(G_loop(ith),1) = 0;
% end
% 
% for ith = 1:length(D_loop)
%         C(D_loop(ith),D_loop(ith)) = 1;
%         f(D_loop(ith),1) = tp_tvV(D_loop(ith),4);
% end
% 
% for ith = 1:length(Nxm_loop)
%         C(Nxm_loop(ith),n_right(Nxm_loop(ith))) = 1;
%         C(Nxm_loop(ith),Nxm_loop(ith)) = -1;
%         f(Nxm_loop(ith),1) = 0;
% end
% 
% for ith = 1:length(NxM_loop)
%         C(NxM_loop(ith),n_left(NxM_loop(ith))) = 1;
%         C(NxM_loop(ith),NxM_loop(ith)) = -1;
%         f(NxM_loop(ith),1) = 0;
% end
% 
% for ith = 1:length(Nym_loop)
%         C(Nym_loop(ith),n_upper(Nym_loop(ith))) = 1;
%         C(Nym_loop(ith),Nym_loop(ith)) = -1;
%         f(Nym_loop(ith),1) = 0;
% end
% 
% for ith = 1:length(NyM_loop)
%         C(NyM_loop(ith),n_lower(NyM_loop(ith))) = 1;
%         C(NyM_loop(ith),NyM_loop(ith)) = -1;
%         f(NyM_loop(ith),1) = 0;
% end
% 
% u=C\f;

% Grid spacing
    d = tp_tvV(1, 2) - tp_tvV(2, 2);

    % Hydraulic conductivity
    K = exp(tp_Y(:, 3));

    % Estimate number of non-zero elements for sparse matrix
    nnz_est = length(G_loop) * 5 + length(D_loop) + 2 * (length(Nxm_loop) + length(NxM_loop) + length(Nym_loop) + length(NyM_loop));
    
    % Preallocate arrays for sparse matrix construction
    rows = zeros(nnz_est, 1);
    cols = zeros(nnz_est, 1);
    vals = zeros(nnz_est, 1);
    f = zeros(size(tp_tvV, 1), 1);

    nnz_idx = 0;

    % Interior nodes (G_loop)
    for ith = 1:length(G_loop)
        node = G_loop(ith);

        % Compute geometric means of conductivities
        Kavg_right = geomean([K(node), K(n_right(node))]);
        Kavg_left = geomean([K(node), K(n_left(node))]);
        Kavg_upper = geomean([K(node), K(n_upper(node))]);
        Kavg_lower = geomean([K(node), K(n_lower(node))]);

        % Fill sparse matrix entries
        nnz_idx = nnz_idx + 5;
        rows(nnz_idx-4:nnz_idx) = [node, node, node, node, node];
        cols(nnz_idx-4:nnz_idx) = [n_right(node), n_left(node), n_upper(node), n_lower(node), node];
        vals(nnz_idx-4:nnz_idx) = [Kavg_right, Kavg_left, Kavg_upper, Kavg_lower, -(Kavg_right + Kavg_left + Kavg_upper + Kavg_lower)] / d^2;

        % Source term
        f(node) = 0;
    end

    % Dirichlet boundary nodes (D_loop)
    for ith = 1:length(D_loop)
        node = D_loop(ith);

        % Fill sparse matrix entries
        nnz_idx = nnz_idx + 1;
        rows(nnz_idx) = node;
        cols(nnz_idx) = node;
        vals(nnz_idx) = 1;

        % Boundary value
        f(node) = tp_tvV(node, 4);
    end

    % Neumann boundary nodes (Nxm_loop, NxM_loop, Nym_loop, NyM_loop)
    for ith = 1:length(Nxm_loop)
        node = Nxm_loop(ith);
        nnz_idx = nnz_idx + 2;
        rows(nnz_idx-1:nnz_idx) = [node, node];
        cols(nnz_idx-1:nnz_idx) = [n_right(node), node];
        vals(nnz_idx-1:nnz_idx) = [1, -1];
        f(node) = 0;
    end

    for ith = 1:length(NxM_loop)
        node = NxM_loop(ith);
        nnz_idx = nnz_idx + 2;
        rows(nnz_idx-1:nnz_idx) = [node, node];
        cols(nnz_idx-1:nnz_idx) = [n_left(node), node];
        vals(nnz_idx-1:nnz_idx) = [1, -1];
        f(node) = 0;
    end

    for ith = 1:length(Nym_loop)
        node = Nym_loop(ith);
        nnz_idx = nnz_idx + 2;
        rows(nnz_idx-1:nnz_idx) = [node, node];
        cols(nnz_idx-1:nnz_idx) = [n_upper(node), node];
        vals(nnz_idx-1:nnz_idx) = [1, -1];
        f(node) = 0;
    end

    for ith = 1:length(NyM_loop)
        node = NyM_loop(ith);
        nnz_idx = nnz_idx + 2;
        rows(nnz_idx-1:nnz_idx) = [node, node];
        cols(nnz_idx-1:nnz_idx) = [n_lower(node), node];
        vals(nnz_idx-1:nnz_idx) = [1, -1];
        f(node) = 0;
    end

    % Assemble sparse matrix
    rows = rows(1:nnz_idx);
    cols = cols(1:nnz_idx);
    vals = vals(1:nnz_idx);
    C = sparse(rows, cols, vals, size(tp_tvV, 1), size(tp_tvV, 1));

    % Solve linear system
    u = C \ f;