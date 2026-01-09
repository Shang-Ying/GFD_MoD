% function [CYY,CYu,Cuu] = fct_Steady_2ndMEHead_HLG_FDM(tpY_cov,tp_tv,u_mean,G_loop,D_loop,Nxm_loop,NxM_loop,Nym_loop,NyM_loop,n_right,n_left,n_upper,n_lower)
% 
% CYY = tpY_cov(:,4:end);
% CYu = zeros(size(tpY_cov,1),size(tpY_cov,1));
% Cuu = zeros(size(tpY_cov,1),size(tpY_cov,1));
% 
% d = tp_tv(1,2)-tp_tv(2,2);
% 
% %% x
% parfor dx = 1:size(tpY_cov,1)
%     C_1 = sparse(size(tpY_cov,1),size(tpY_cov,1));
%     f_1 = zeros(size(tpY_cov,1),1);
%     for ith = 1:length(G_loop)
%         C_1(G_loop(ith),n_right(G_loop(ith))) = 1/(d^2);
%         C_1(G_loop(ith),n_left(G_loop(ith))) = 1/(d^2);
%         C_1(G_loop(ith),n_upper(G_loop(ith))) = 1/(d^2);
%         C_1(G_loop(ith),n_lower(G_loop(ith))) = 1/(d^2);
%         C_1(G_loop(ith),G_loop(ith)) = -4/(d^2);
%         dYdx_1 = (tpY_cov(n_right(G_loop(ith)),3)-tpY_cov(n_left(G_loop(ith)),3))/(2*d);
%         dYdy_1 = (tpY_cov(n_upper(G_loop(ith)),3)-tpY_cov(n_lower(G_loop(ith)),3))/(2*d);
%         dCYYdx = (CYY(dx,n_right(G_loop(ith)))'-CYY(dx,n_left(G_loop(ith)))')/(2*d);
%         dCYYdy = (CYY(dx,n_upper(G_loop(ith)))'-CYY(dx,n_lower(G_loop(ith)))')/(2*d);
%         dudx_1 = (u_mean(n_right(G_loop(ith)),1)-u_mean(n_left(G_loop(ith)),1))/(2*d);
%         dudy_1 = (u_mean(n_upper(G_loop(ith)),1)-u_mean(n_lower(G_loop(ith)),1))/(2*d);
%         J1_1 = -1*dudx_1;
%         J2_1 = -1*dudy_1;
%         f_1(G_loop(ith),1) = (J1_1*dCYYdx)+(J2_1*dCYYdy);
%     end
%     for ith = 1:length(D_loop)
%             C_1(D_loop(ith),D_loop(ith)) = 1;
%             f_1(D_loop(ith),1) = 0;
%     end    
%     for ith = 1:length(Nxm_loop)
%             C_1(Nxm_loop(ith),n_right(Nxm_loop(ith))) = 1;
%             C_1(Nxm_loop(ith),Nxm_loop(ith)) = -1;
%             f_1(Nxm_loop(ith),1) = 0;
%     end    
%     for ith = 1:length(NxM_loop)
%             C_1(NxM_loop(ith),n_left(NxM_loop(ith))) = 1;
%             C_1(NxM_loop(ith),NxM_loop(ith)) = -1;
%             f_1(NxM_loop(ith),1) = 0;
%     end    
%     for ith = 1:length(Nym_loop)
%             C_1(Nym_loop(ith),n_upper(Nym_loop(ith))) = 1;
%             C_1(Nym_loop(ith),Nym_loop(ith)) = -1;
%             f_1(Nym_loop(ith),1) = 0;
%     end  
%     for ith = 1:length(NyM_loop)
%             C_1(NyM_loop(ith),n_lower(NyM_loop(ith))) = 1;
%             C_1(NyM_loop(ith),NyM_loop(ith)) = -1;
%             f_1(NyM_loop(ith),1) = 0;
%     end
%     CYu(dx,:)=(C_1\f_1)';
% end
% 
% 
% %% chi
% parfor dchi = 1:size(tpY_cov,1)
%     C_2 = sparse(size(tpY_cov,1),size(tpY_cov,1));
%     f_2 = zeros(size(tpY_cov,1),1);
%     for ith = 1:length(G_loop)
%         C_2(G_loop(ith),n_right(G_loop(ith))) = 1/(d^2);
%         C_2(G_loop(ith),n_left(G_loop(ith))) = 1/(d^2);
%         C_2(G_loop(ith),n_upper(G_loop(ith))) = 1/(d^2);
%         C_2(G_loop(ith),n_lower(G_loop(ith))) = 1/(d^2);
%         C_2(G_loop(ith),G_loop(ith)) = -4/(d^2);
%         dYdx_2 = (tpY_cov(n_right(G_loop(ith)),3)-tpY_cov(n_left(G_loop(ith)),3))/(2*d);
%         dYdy_2 = (tpY_cov(n_upper(G_loop(ith)),3)-tpY_cov(n_lower(G_loop(ith)),3))/(2*d);
%         dCYudx = (CYu(n_right(G_loop(ith)),dchi)-CYu(n_left(G_loop(ith)),dchi))/(2*d);
%         dCYudy = (CYu(n_upper(G_loop(ith)),dchi)-CYu(n_lower(G_loop(ith)),dchi))/(2*d);
%         dudx_2 = (u_mean(n_right(G_loop(ith)),1)-u_mean(n_left(G_loop(ith)),1))/(2*d);
%         dudy_2 = (u_mean(n_upper(G_loop(ith)),1)-u_mean(n_lower(G_loop(ith)),1))/(2*d);
%         J1_2 = -1*dudx_2;
%         J2_2 = -1*dudy_2;
%         f_2(G_loop(ith),1) = (J1_2*dCYudx)+(J2_2*dCYudy);
%     end
%     for ith = 1:length(D_loop)
%             C_2(D_loop(ith),D_loop(ith)) = 1;
%             f_2(D_loop(ith),1) = 0;
%     end    
%     for ith = 1:length(Nxm_loop)
%             C_2(Nxm_loop(ith),n_right(Nxm_loop(ith))) = 1;
%             C_2(Nxm_loop(ith),Nxm_loop(ith)) = -1;
%             f_2(Nxm_loop(ith),1) = 0;
%     end    
%     for ith = 1:length(NxM_loop)
%             C_2(NxM_loop(ith),n_left(NxM_loop(ith))) = 1;
%             C_2(NxM_loop(ith),NxM_loop(ith)) = -1;
%             f_2(NxM_loop(ith),1) = 0;
%     end    
%     for ith = 1:length(Nym_loop)
%             C_2(Nym_loop(ith),n_upper(Nym_loop(ith))) = 1;
%             C_2(Nym_loop(ith),Nym_loop(ith)) = -1;
%             f_2(Nym_loop(ith),1) = 0;
%     end  
%     for ith = 1:length(NyM_loop)
%             C_2(NyM_loop(ith),n_lower(NyM_loop(ith))) = 1;
%             C_2(NyM_loop(ith),NyM_loop(ith)) = -1;
%             f_2(NyM_loop(ith),1) = 0;
%     end
%     Cuu(:,dchi)=(C_2\f_2);
% end
function [CYY, CYu, Cuu] = fct_Steady_2ndMEHead_HLG_FDM(tpY_cov, tp_tv, u_mean, G_loop, D_loop, Nxm_loop, NxM_loop, Nym_loop, NyM_loop, n_right, n_left, n_upper, n_lower)

CYY = tpY_cov(:, 4:end);
CYu = zeros(size(tpY_cov, 1), size(tpY_cov, 1));
Cuu = zeros(size(tpY_cov, 1), size(tpY_cov, 1));

d = tp_tv(1, 2) - tp_tv(2, 2);

%% x
parfor dx = 1:size(tpY_cov, 1)
    rows = [];
    cols = [];
    vals = [];
    f_1 = zeros(size(tpY_cov, 1), 1);
    
    % Assemble sparse matrix entries
    for ith = 1:length(G_loop)
        g_idx = G_loop(ith);
        rows = [rows; g_idx; g_idx; g_idx; g_idx; g_idx];
        cols = [cols; n_right(g_idx); n_left(g_idx); n_upper(g_idx); n_lower(g_idx); g_idx];
        vals = [vals; 1/(d^2); 1/(d^2); 1/(d^2); 1/(d^2); -4/(d^2)];
        
        dCYYdx = (CYY(dx, n_right(g_idx)) - CYY(dx, n_left(g_idx))) / (2 * d);
        dCYYdy = (CYY(dx, n_upper(g_idx)) - CYY(dx, n_lower(g_idx))) / (2 * d);
        dudx_1 = (u_mean(n_right(g_idx)) - u_mean(n_left(g_idx))) / (2 * d);
        dudy_1 = (u_mean(n_upper(g_idx)) - u_mean(n_lower(g_idx))) / (2 * d);
        
        J1_1 = -1 * dudx_1;
        J2_1 = -1 * dudy_1;
        
        f_1(g_idx) = (J1_1 * dCYYdx) + (J2_1 * dCYYdy);
    end
    
    for ith = 1:length(D_loop)
        d_idx = D_loop(ith);
        rows = [rows; d_idx];
        cols = [cols; d_idx];
        vals = [vals; 1];
        f_1(d_idx) = 0;
    end
    
    for ith = 1:length(Nxm_loop)
        nxm_idx = Nxm_loop(ith);
        rows = [rows; nxm_idx; nxm_idx];
        cols = [cols; n_right(nxm_idx); nxm_idx];
        vals = [vals; 1; -1];
        f_1(nxm_idx) = 0;
    end
    
    for ith = 1:length(NxM_loop)
        nxM_idx = NxM_loop(ith);
        rows = [rows; nxM_idx; nxM_idx];
        cols = [cols; n_left(nxM_idx); nxM_idx];
        vals = [vals; 1; -1];
        f_1(nxM_idx) = 0;
    end
    
    for ith = 1:length(Nym_loop)
        nym_idx = Nym_loop(ith);
        rows = [rows; nym_idx; nym_idx];
        cols = [cols; n_upper(nym_idx); nym_idx];
        vals = [vals; 1; -1];
        f_1(nym_idx) = 0;
    end
    
    for ith = 1:length(NyM_loop)
        nyM_idx = NyM_loop(ith);
        rows = [rows; nyM_idx; nyM_idx];
        cols = [cols; n_lower(nyM_idx); nyM_idx];
        vals = [vals; 1; -1];
        f_1(nyM_idx) = 0;
    end
    
    % Solve sparse system
    C_1 = sparse(rows, cols, vals, size(tpY_cov, 1), size(tpY_cov, 1));
    CYu(dx, :) = (C_1 \ f_1)';
end

%% chi
parfor dchi = 1:size(tpY_cov, 1)
    rows = [];
    cols = [];
    vals = [];
    f_2 = zeros(size(tpY_cov, 1), 1);
    
    % Assemble sparse matrix entries
    for ith = 1:length(G_loop)
        g_idx = G_loop(ith);
        rows = [rows; g_idx; g_idx; g_idx; g_idx; g_idx];
        cols = [cols; n_right(g_idx); n_left(g_idx); n_upper(g_idx); n_lower(g_idx); g_idx];
        vals = [vals; 1/(d^2); 1/(d^2); 1/(d^2); 1/(d^2); -4/(d^2)];
        
        dCYudx = (CYu(n_right(g_idx), dchi) - CYu(n_left(g_idx), dchi)) / (2 * d);
        dCYudy = (CYu(n_upper(g_idx), dchi) - CYu(n_lower(g_idx), dchi)) / (2 * d);
        dudx_2 = (u_mean(n_right(g_idx)) - u_mean(n_left(g_idx))) / (2 * d);
        dudy_2 = (u_mean(n_upper(g_idx)) - u_mean(n_lower(g_idx))) / (2 * d);
        
        J1_2 = -1 * dudx_2;
        J2_2 = -1 * dudy_2;
        
        f_2(g_idx) = (J1_2 * dCYudx) + (J2_2 * dCYudy);
    end
    
    for ith = 1:length(D_loop)
        d_idx = D_loop(ith);
        rows = [rows; d_idx];
        cols = [cols; d_idx];
        vals = [vals; 1];
        f_2(d_idx) = 0;
    end
    
    for ith = 1:length(Nxm_loop)
        nxm_idx = Nxm_loop(ith);
        rows = [rows; nxm_idx; nxm_idx];
        cols = [cols; n_right(nxm_idx); nxm_idx];
        vals = [vals; 1; -1];
        f_2(nxm_idx) = 0;
    end
    
    for ith = 1:length(NxM_loop)
        nxM_idx = NxM_loop(ith);
        rows = [rows; nxM_idx; nxM_idx];
        cols = [cols; n_left(nxM_idx); nxM_idx];
        vals = [vals; 1; -1];
        f_2(nxM_idx) = 0;
    end
    
    for ith = 1:length(Nym_loop)
        nym_idx = Nym_loop(ith);
        rows = [rows; nym_idx; nym_idx];
        cols = [cols; n_upper(nym_idx); nym_idx];
        vals = [vals; 1; -1];
        f_2(nym_idx) = 0;
    end
    
    for ith = 1:length(NyM_loop)
        nyM_idx = NyM_loop(ith);
        rows = [rows; nyM_idx; nyM_idx];
        cols = [cols; n_lower(nyM_idx); nyM_idx];
        vals = [vals; 1; -1];
        f_2(nyM_idx) = 0;
    end
    
    % Solve sparse system
    C_2 = sparse(rows, cols, vals, size(tpY_cov, 1), size(tpY_cov, 1));
    Cuu(:, dchi) = (C_2 \ f_2);
end