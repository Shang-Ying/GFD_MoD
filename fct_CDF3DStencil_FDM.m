function [G_loop,D_loop,Nxm_loop,NxM_loop,Nym_loop,NyM_loop,Nzm_loop,NzM_loop,n_right,n_left,n_front,n_back,n_upper,n_lower] = fct_CDF3DStencil_FDM(tpH_sort,xy_res,H_res)

G_In = (tpH_sort(:,4) == 0);
D_Bd = (tpH_sort(:,4) == 1);
N_Bd_xm = (tpH_sort(:,4) == 2 & tpH_sort(:,1) == 0); % x = 0 (Dirichlet 1)
N_Bd_xM = (tpH_sort(:,4) == 2 & tpH_sort(:,1) == 1); % x = 1 (Dirichlet 0)
N_Bd_ym = (tpH_sort(:,4) == 2 & tpH_sort(:,2) == 0); % y = 0 (Neumann 0)
N_Bd_yM = (tpH_sort(:,4) == 2 & tpH_sort(:,2) == 1); % y = 1 (Neumann 0)
N_Bd_zm = (tpH_sort(:,4) == 2 & tpH_sort(:,3) == 0.1); % z = 0.1 (Dirichlet 0)
N_Bd_zM = (tpH_sort(:,4) == 2 & tpH_sort(:,3) == 1.1); % z = 1.1 (Dirichlet 1)

G_loop = find(G_In);
D_loop = find(D_Bd);
Nxm_loop = find(N_Bd_xm);
NxM_loop = find(N_Bd_xM);
Nym_loop = find(N_Bd_ym);
NyM_loop = find(N_Bd_yM);
Nzm_loop = find(N_Bd_zm);
NzM_loop = find(N_Bd_zM);

% G_In = (tpH(:,4) == 0);
% D_Bd = (tpH(:,4) == 1);
% N_Bd_xm = (tpH(:,4) == 2 & tpH(:,1) == 0); % x = 0 (Dirichlet 1)
% N_Bd_xM = (tpH(:,4) == 2 & tpH(:,1) == 1); % x = 1 (Dirichlet 0)
% N_Bd_ym = (tpH(:,4) == 2 & tpH(:,2) == 0); % y = 0 (Neumann 0)
% N_Bd_yM = (tpH(:,4) == 2 & tpH(:,2) == 1); % y = 1 (Neumann 0)
% N_Bd_zm = (tpH(:,4) == 2 & tpH(:,3) == 0.1); % z = 0.1 (Dirichlet 0)
% N_Bd_zM = (tpH(:,4) == 2 & tpH(:,3) == 1.1); % z = 1.1 (Dirichlet 1)
% 
% G_loop = find(G_In);
% D_loop = find(D_Bd);
% Nxm_loop = find(N_Bd_xm);
% NxM_loop = find(N_Bd_xM);
% Nym_loop = find(N_Bd_ym);
% NyM_loop = find(N_Bd_yM);
% Nzm_loop = find(N_Bd_zm);
% NzM_loop = find(N_Bd_zM);
% 
% d = tpH(1,2)-tpH(2,2); % x & y

for ith = 1:length(G_loop)
    n_right(G_loop(ith)) = G_loop(ith)+xy_res;
    n_left(G_loop(ith)) = G_loop(ith)-xy_res;
    n_front(G_loop(ith)) = G_loop(ith)+1;
    n_back(G_loop(ith)) = G_loop(ith)+1;
    n_upper(G_loop(ith)) = G_loop(ith)+xy_res*xy_res;
    n_lower(G_loop(ith)) = G_loop(ith)-xy_res*xy_res;
    % n_right(G_loop(ith)) = find(ismembertol(tpH(:,1:3), [tpH(G_loop(ith),1)+d tpH(G_loop(ith),2) tpH(G_loop(ith),3)],'ByRows',true));
    % n_left(G_loop(ith)) = find(ismembertol(tpH(:,1:3), [tpH(G_loop(ith),1)-d tpH(G_loop(ith),2) tpH(G_loop(ith),3)],'ByRows',true));
    % n_front(G_loop(ith)) = find(ismembertol(tpH(:,1:3), [tpH(G_loop(ith),1) tpH(G_loop(ith),2)+d tpH(G_loop(ith),3)],'ByRows',true));
    % n_back(G_loop(ith)) = find(ismembertol(tpH(:,1:3), [tpH(G_loop(ith),1) tpH(G_loop(ith),2)-d tpH(G_loop(ith),3)],'ByRows',true));
    % n_upper(G_loop(ith)) = find(ismembertol(tpH(:,1:3), [tpH(G_loop(ith),1) tpH(G_loop(ith),2) tpH(G_loop(ith),3)+dH],'ByRows',true));
    % n_lower(G_loop(ith)) = find(ismembertol(tpH(:,1:3), [tpH(G_loop(ith),1) tpH(G_loop(ith),2) tpH(G_loop(ith),3)-dH],'ByRows',true));
end

for ith = 1:length(Nxm_loop)
    n_right(Nxm_loop(ith)) = Nxm_loop(ith)+xy_res;
end

for ith = 1:length(NxM_loop)
    n_left(NxM_loop(ith)) = NxM_loop(ith)-xy_res;
end

for ith = 1:length(Nym_loop)
    n_front(Nym_loop(ith)) = Nym_loop(ith)+1;
end

for ith = 1:length(NyM_loop)
    n_back(NyM_loop(ith)) = NyM_loop(ith)-1;
end

for ith = 1:length(Nzm_loop)
    n_upper(Nzm_loop(ith)) = Nzm_loop(ith)+xy_res*xy_res;
end

for ith = 1:length(NzM_loop)
    n_lower(Nzm_loop(ith)) = Nzm_loop(ith)-xy_res*xy_res;
end