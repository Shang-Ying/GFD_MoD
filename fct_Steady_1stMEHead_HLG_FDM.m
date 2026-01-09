function [u] = fct_Steady_1stMEHead_HLG_FDM(tpY_cov,tp_tvV,G_loop,D_loop,Nxm_loop,NxM_loop,Nym_loop,NyM_loop,n_right,n_left,n_upper,n_lower)
C = sparse((size(tpY_cov,1)),(size(tpY_cov,1)));

d = tp_tvV(1,2)-tp_tvV(2,2);

for ith = 1:length(G_loop)
    C(G_loop(ith),n_right(G_loop(ith))) = 1;
    C(G_loop(ith),n_left(G_loop(ith))) = 1;
    C(G_loop(ith),n_upper(G_loop(ith))) = 1;
    C(G_loop(ith),n_lower(G_loop(ith))) = 1;
    C(G_loop(ith),G_loop(ith)) = -4;
    f(G_loop(ith),1) = 0;
end

for ith = 1:length(D_loop)
        C(D_loop(ith),D_loop(ith)) = 1;
        f(D_loop(ith),1) = tp_tvV(D_loop(ith),4);
end

for ith = 1:length(Nxm_loop)
        C(Nxm_loop(ith),n_right(Nxm_loop(ith))) = 1;
        C(Nxm_loop(ith),Nxm_loop(ith)) = -1;
        f(Nxm_loop(ith),1) = 0;
end

for ith = 1:length(NxM_loop)
        C(NxM_loop(ith),n_left(NxM_loop(ith))) = 1;
        C(NxM_loop(ith),NxM_loop(ith)) = -1;
        f(NxM_loop(ith),1) = 0;
end

for ith = 1:length(Nym_loop)
        C(Nym_loop(ith),n_upper(Nym_loop(ith))) = 1;
        C(Nym_loop(ith),Nym_loop(ith)) = -1;
        f(Nym_loop(ith),1) = 0;
end

for ith = 1:length(NyM_loop)
        C(NyM_loop(ith),n_lower(NyM_loop(ith))) = 1;
        C(NyM_loop(ith),NyM_loop(ith)) = -1;
        f(NyM_loop(ith),1) = 0;
end

u=C\f;