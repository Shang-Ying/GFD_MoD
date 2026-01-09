function [G_loop,D_loop,N_loop] = fct_2DStencil_GFDM(tp_tvV)

G_In = (tp_tvV(:,3) == 0);
D_Bd = (tp_tvV(:,3) == 1);
N_Bd = (tp_tvV(:,3) == 2);

G_loop = find(G_In);
D_loop = find(D_Bd);
N_loop = find(N_Bd);
