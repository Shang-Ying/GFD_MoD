% function [CYY,CYu,Cuu] = fct_Steady_2ndMEHead_HLG(Wcoef,idx,ns,tp_cov,tp_tvV,u_mean,G_loop,D_loop,N_loop)
% 
% CYY = tp_cov(:,4:end);
% 
% % CYu = zeros(size(tp_cov,1),size(tp_cov,1));
% % Cuu = zeros(size(tp_cov,1),size(tp_cov,1));
% Wcoef_1 = squeeze(Wcoef(1,:,:))';
% Wcoef_2 = squeeze(Wcoef(2,:,:))';
% Wcoef_3 = squeeze(Wcoef(3,:,:))';
% % Wcoef_4 = squeeze(Wcoef(4,:,:))';
% Wcoef_5 = squeeze(Wcoef(5,:,:))';
% 
% parfor dx = 1:size(idx,1)
%     % C_1 = zeros(size(tp_cov,1),size(tp_cov,1));
%     f_1 = zeros(size(tp_cov,1),1);
%     nzmax = size(tp_cov,1)*(ns+1);
%     rows = zeros(1, nzmax);
%     cols = zeros(1, nzmax);
%     vals = zeros(1, nzmax);
%     counter = 0;
%     for ith = 1:length(G_loop)
%         idxGl = idx(G_loop(ith),:);
%         W_1 = Wcoef_1(G_loop(ith),:);
%         W_2 = Wcoef_2(G_loop(ith),:);
%         W_3 = Wcoef_3(G_loop(ith),:);
%         % W_4 = Wcoef_4(G_loop(ith),:);
%         W_5 = Wcoef_5(G_loop(ith),:);
%         % Yidx = tp_cov(idxGl,3);
%         CYYdchi = CYY(dx,idxGl);
%         umidx = u_mean(idxGl,1);
%         % dYdx_1 = W_1*Yidx;
%         % dYdy_1 = W_2*Yidx;
%         dCYYdx = W_1*CYYdchi';
%         dCYYdy = W_2*CYYdchi';
%         dudx_1 = W_1*umidx;
%         dudy_1 = W_2*umidx;
%         J1_1 = -1*dudx_1;
%         J2_1 = -1*dudy_1;
%         % Cia = W_3 + W_5;
%         % Cib = dYdx_1*W_1;
%         % Cic = dYdy_1*W_2;
%         % C_1(G_loop(ith),(idx(G_loop(ith),:))) = Cia + Cib + Cic;
%         % C_1(G_loop(ith),idxGl) = W_3 + W_5;
%         len = length(idxGl);
%         rows(counter+1:counter+len) = G_loop(ith);
%         cols(counter+1:counter+len) = idxGl;
%         vals(counter+1:counter+len) = W_3 + W_5;
%         counter = counter + len;
%         f_1(G_loop(ith),1) = (J1_1*dCYYdx)+(J2_1*dCYYdy);
%     end
% 
%     for ith = 1:length(D_loop)
%         counter = counter + 1;
%         rows(counter) = D_loop(ith);
%         cols(counter) = D_loop(ith);
%         vals(counter) = 1;
%         % C_1(D_loop(ith),D_loop(ith)) = 1;
%         f_1(D_loop(ith),1) = 0;
%     end
%     for ith = 1:length(N_loop)
%         idxNl = idx(N_loop(ith), :);
%         len = length(idxNl);
%         rows(counter+1:counter+len) = N_loop(ith);
%         cols(counter+1:counter+len) = idxNl;
%         vals(counter+1:counter+len) = (abs(tp_tvV(N_loop(ith), 5)) .* Wcoef_1(N_loop(ith), :)) + ...
%                                       (abs(tp_tvV(N_loop(ith), 6)) .* Wcoef_2(N_loop(ith), :));
%         counter = counter + len;
%         % C_1(N_loop(ith),(idx(N_loop(ith),:))) = (abs(tp_tvV(N_loop(ith),5)).*Wcoef_1(N_loop(ith),:)) + ...
%         %     (abs(tp_tvV(N_loop(ith),6)).*Wcoef_2(N_loop(ith),:));
%         f_1(N_loop(ith),1) = 0;
%     end
% 
%     % C_1 = sparse(C_1);
%     C_1 = sparse(rows(1:counter), cols(1:counter), vals(1:counter), size(tp_cov, 1), size(tp_cov, 1));
%     f_1 = sparse(f_1);
%     CYu(dx,:)=(C_1\f_1)';
% end
% 
% parfor dchi = 1:size(idx,1)
%     % C_2 = zeros(size(tp_cov,1),size(tp_cov,1));
%     f_2 = zeros(size(tp_cov,1),1);
%     nzmax = size(tp_cov,1)*(ns+1);
%     rows = zeros(1, nzmax);
%     cols = zeros(1, nzmax);
%     vals = zeros(1, nzmax);
%     counter = 0;
%     for ith = 1:length(G_loop)
%         idxGl = idx(G_loop(ith),:);
%         W_1 = Wcoef_1(G_loop(ith),:);
%         W_2 = Wcoef_2(G_loop(ith),:);
%         W_3 = Wcoef_3(G_loop(ith),:);
%         % W_4 = Wcoef_4(G_loop(ith),:);
%         W_5 = Wcoef_5(G_loop(ith),:);
%         CYudchi = CYu(idxGl,dchi);
%         umidx = u_mean(idxGl,1);
%         % dYdx_2 = (Wcoef(1,1:ns,G_loop(ith))*tp_cov((idx(G_loop(ith),1:ns)),3));
%         % dYdy_2 = (Wcoef(2,1:ns,G_loop(ith))*tp_cov((idx(G_loop(ith),1:ns)),3));
%         dCYudx = W_1*CYudchi;
%         dCYudy = W_2*CYudchi;
%         dudx_2 = (W_1*umidx);
%         dudy_2 = (W_2*umidx);
%         J1_2 = -1*dudx_2;
%         J2_2 = -1*dudy_2;
%         % Cia = Wcoef(3,:,G_loop(ith)) + Wcoef(5,:,G_loop(ith));
%         % Cib = dYdx_2*Wcoef(1,:,G_loop(ith));
%         % Cic = dYdy_2*Wcoef(2,:,G_loop(ith));
%         % C_2(G_loop(ith),(idx(G_loop(ith),:))) = Cia + Cib + Cic;
%         % C_2(G_loop(ith),idxGl) = W_3 + W_5;
%         len = length(idxGl);
%         rows(counter+1:counter+len) = G_loop(ith);
%         cols(counter+1:counter+len) = idxGl;
%         vals(counter+1:counter+len) = W_3 + W_5;
%         counter = counter + len;
%         f_2(G_loop(ith),1) = (J1_2*dCYudx)+(J2_2*dCYudy);
%     end
%     for ith = 1:length(D_loop)
%         counter = counter + 1;
%         rows(counter) = D_loop(ith);
%         cols(counter) = D_loop(ith);
%         vals(counter) = 1;
%         % C_2(D_loop(ith),D_loop(ith)) = 1;
%         f_2(D_loop(ith),1) = 0;
%     end
%     for ith = 1:length(N_loop)
%         idxNl = idx(N_loop(ith), :);
%         len = length(idxNl);
%         rows(counter+1:counter+len) = N_loop(ith);
%         cols(counter+1:counter+len) = idxNl;
%         vals(counter+1:counter+len) = (abs(tp_tvV(N_loop(ith), 5)) .* Wcoef_1(N_loop(ith), :)) + ...
%                                       (abs(tp_tvV(N_loop(ith), 6)) .* Wcoef_2(N_loop(ith), :));
%         counter = counter + len;
%         % C_2(N_loop(ith),(idx(N_loop(ith),:))) = abs(tp_tvV(N_loop(ith),5)).*Wcoef_1(N_loop(ith),:) + ...
%         %     abs(tp_tvV(N_loop(ith),6)).*Wcoef_2(N_loop(ith),:);
%         f_2(N_loop(ith),1) = 0;
%     end
%     % C_2 = sparse(C_2);
%     C_2 = sparse(rows(1:counter), cols(1:counter), vals(1:counter), size(tp_cov, 1), size(tp_cov, 1));
%     f_2 = sparse(f_2);
%     Cuu(:,dchi)=(C_2\f_2);
% end
function [CYY, CYu, Cuu] = fct_Steady_2ndMEHead_HLG(Wcoef, idx, ns, tp_cov, tp_tvV, u_mean, G_loop, D_loop, N_loop)

% Initialize Outputs
num_pts = size(tp_cov, 1);
% num_eqs = size(idx, 1);
CYY = tp_cov(:, 4:end);
CYu = zeros(num_pts, num_pts);
Cuu = zeros(num_pts, num_pts);

% Extract Coefficients for Parallel Use
Wcoef_1 = squeeze(Wcoef(1, :, :))';
Wcoef_2 = squeeze(Wcoef(2, :, :))';
Wcoef_3 = squeeze(Wcoef(3, :, :))';
Wcoef_5 = squeeze(Wcoef(5, :, :))';

% Process CYu
parfor dx = 1:num_pts
    CYu(dx, :) = computeCYu(dx, G_loop, D_loop, N_loop, tp_cov, tp_tvV, ...
        Wcoef_1, Wcoef_2, Wcoef_3, Wcoef_5, idx, CYY, u_mean, ns);
end

% Process Cuu
parfor dchi = 1:num_pts
    Cuu(:, dchi) = computeCuu(dchi, G_loop, D_loop, N_loop, tp_cov, tp_tvV, ...
        Wcoef_1, Wcoef_2, Wcoef_3, Wcoef_5, idx, CYu, u_mean, ns);
end
end

function f_out = computeCYu(dx, G_loop, D_loop, N_loop, tp_cov, tp_tvV, ...
    Wcoef_1, Wcoef_2, Wcoef_3, Wcoef_5, idx, CYY, u_mean, ns)

num_pts = size(tp_cov, 1);
nzmax = num_pts * (ns + 1);
rows = zeros(1, nzmax);
cols = zeros(1, nzmax);
vals = zeros(1, nzmax);
f_1 = zeros(num_pts, 1);
counter = 0;

% Governing Equations Loop
for ith = 1:length(G_loop)
    current_G = G_loop(ith);
    idxGl = idx(current_G, :);
    W_1 = Wcoef_1(current_G, :);
    W_2 = Wcoef_2(current_G, :);
    W_3 = Wcoef_3(current_G, :);
    W_5 = Wcoef_5(current_G, :);
    CYYdchi = CYY(dx, idxGl);
    umidx = u_mean(idxGl, 1);

    % Compute Gradients
    dCYYdx = W_1 * CYYdchi';
    dCYYdy = W_2 * CYYdchi';
    dudx_1 = W_1 * umidx;
    dudy_1 = W_2 * umidx;

    % Update f_1
    f_1(current_G) = -dudx_1 * dCYYdx - dudy_1 * dCYYdy;

    % Sparse Matrix Assembly
    len = length(idxGl);
    rows(counter+1:counter+len) = current_G;
    cols(counter+1:counter+len) = idxGl;
    vals(counter+1:counter+len) = W_3 + W_5;
    counter = counter + len;
end

% Dirichlet Boundary
for ith = 1:length(D_loop)
    counter = counter + 1;
    rows(counter) = D_loop(ith);
    cols(counter) = D_loop(ith);
    vals(counter) = 1;
    f_1(D_loop(ith)) = 0;
end

% Neumann Boundary
for ith = 1:length(N_loop)
    idxNl = idx(N_loop(ith), :);
    len = length(idxNl);
    rows(counter+1:counter+len) = N_loop(ith);
    cols(counter+1:counter+len) = idxNl;
    vals(counter+1:counter+len) = (tp_tvV(N_loop(ith), 5)) .* Wcoef_1(N_loop(ith), :) + ...
        (tp_tvV(N_loop(ith), 6)) .* Wcoef_2(N_loop(ith), :);
    counter = counter + len;
    f_1(N_loop(ith)) = 0;
end

% Sparse Matrix Solve
C_1 = sparse(rows(1:counter), cols(1:counter), vals(1:counter), num_pts, num_pts);
f_out = sparse(C_1 \ f_1);
end

function f_out = computeCuu(dchi, G_loop, D_loop, N_loop, tp_cov, tp_tvV, ...
    Wcoef_1, Wcoef_2, Wcoef_3, Wcoef_5, idx, CYu, u_mean, ns)

num_pts = size(tp_cov, 1);
nzmax = num_pts * (ns + 1);
rows = zeros(1, nzmax);
cols = zeros(1, nzmax);
vals = zeros(1, nzmax);
f_2 = zeros(num_pts, 1);
counter = 0;

% Governing Equations Loop
for ith = 1:length(G_loop)
    current_G = G_loop(ith);
    idxGl = idx(G_loop(ith), :);
    W_1 = Wcoef_1(current_G, :);
    W_2 = Wcoef_2(current_G, :);
    W_3 = Wcoef_3(current_G, :);
    W_5 = Wcoef_5(current_G, :);
    CYudchi = CYu(idxGl, dchi);
    umidx = u_mean(idxGl, 1);

    % Compute Gradients
    dCYudx = W_1 * CYudchi;
    dCYudy = W_2 * CYudchi;
    dudx_2 = W_1 * umidx;
    dudy_2 = W_2 * umidx;

    % Update f_2
    f_2(current_G) = -dudx_2 * dCYudx - dudy_2 * dCYudy;

    % Sparse Matrix Assembly
    len = length(idxGl);
    rows(counter+1:counter+len) = current_G;
    cols(counter+1:counter+len) = idxGl;
    vals(counter+1:counter+len) = W_3 + W_5;
    counter = counter + len;
end

% Dirichlet Boundary
for ith = 1:length(D_loop)
    counter = counter + 1;
    rows(counter) = D_loop(ith);
    cols(counter) = D_loop(ith);
    vals(counter) = 1;
    f_2(D_loop(ith)) = 0;
end

% Neumann Boundary
for ith = 1:length(N_loop)
    idxNl = idx(N_loop(ith), :);
    len = length(idxNl);
    rows(counter+1:counter+len) = N_loop(ith);
    cols(counter+1:counter+len) = idxNl;
    vals(counter+1:counter+len) = (tp_tvV(N_loop(ith), 5)) .* Wcoef_1(N_loop(ith), :) + ...
        (tp_tvV(N_loop(ith), 6)) .* Wcoef_2(N_loop(ith), :);
    counter = counter + len;
    f_2(N_loop(ith)) = 0;
end

% Sparse Matrix Solve
C_2 = sparse(rows(1:counter), cols(1:counter), vals(1:counter), num_pts, num_pts);
f_out = sparse(C_2 \ f_2);
end