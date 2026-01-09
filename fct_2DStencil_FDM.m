function [G_loop,D_loop,Nxm_loop,NxM_loop,Nym_loop,NyM_loop,n_right,n_left,n_upper,n_lower] = fct_2DStencil_FDM(tp_tvV)

G_In = (tp_tvV(:,3) == 0);
D_Bd = (tp_tvV(:,3) == 1);
N_Bd_xm = (tp_tvV(:,3) == 2 & tp_tvV(:,1) == 0);
N_Bd_xM = (tp_tvV(:,3) == 2 & tp_tvV(:,1) == 1);
N_Bd_ym = (tp_tvV(:,3) == 2 & tp_tvV(:,2) == 0);
N_Bd_yM = (tp_tvV(:,3) == 2 & tp_tvV(:,2) == 1);

G_loop = find(G_In);
D_loop = find(D_Bd);
Nxm_loop = find(N_Bd_xm);
NxM_loop = find(N_Bd_xM);
Nym_loop = find(N_Bd_ym);
NyM_loop = find(N_Bd_yM);

d = tp_tvV(1,2)-tp_tvV(2,2);

n_right = zeros(size(tp_tvV,2),1);
n_left = zeros(size(tp_tvV,2),1);
n_upper = zeros(size(tp_tvV,2),1);
n_lower = zeros(size(tp_tvV,2),1);

for ith = 1:length(G_loop)
    n_right(G_loop(ith)) = find(ismembertol(tp_tvV(:,1:2), [tp_tvV(G_loop(ith),1)+d tp_tvV(G_loop(ith),2)],'ByRows',true));
    n_left(G_loop(ith)) = find(ismembertol(tp_tvV(:,1:2), [tp_tvV(G_loop(ith),1)-d tp_tvV(G_loop(ith),2)],'ByRows',true));
    n_upper(G_loop(ith)) = find(ismembertol(tp_tvV(:,1:2), [tp_tvV(G_loop(ith),1) tp_tvV(G_loop(ith),2)+d],'ByRows',true));
    n_lower(G_loop(ith)) = find(ismembertol(tp_tvV(:,1:2), [tp_tvV(G_loop(ith),1) tp_tvV(G_loop(ith),2)-d],'ByRows',true));
end

for ith = 1:length(Nxm_loop)
    n_right(Nxm_loop(ith)) = find(ismembertol(tp_tvV(:,1:2), [tp_tvV(Nxm_loop(ith),1)+d tp_tvV(Nxm_loop(ith),2)],'ByRows',true));
end

for ith = 1:length(NxM_loop)
    n_left(NxM_loop(ith)) = find(ismembertol(tp_tvV(:,1:2), [tp_tvV(NxM_loop(ith),1)-d tp_tvV(NxM_loop(ith),2)],'ByRows',true));
end

for ith = 1:length(Nym_loop)
    n_upper(Nym_loop(ith)) = find(ismembertol(tp_tvV(:,1:2), [tp_tvV(Nym_loop(ith),1) tp_tvV(Nym_loop(ith),2)+d],'ByRows',true));
end

for ith = 1:length(NyM_loop)
    n_lower(NyM_loop(ith)) = find(ismembertol(tp_tvV(:,1:2), [tp_tvV(NyM_loop(ith),1) tp_tvV(NyM_loop(ith),2)-d],'ByRows',true));
end