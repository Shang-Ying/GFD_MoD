function [u] = fct_Steady_1stMEHead_HLG(Wcoef,idx,ns,tp_cov,tp_tvV,G_loop,D_loop,N_loop)
C = sparse((size(tp_cov,1)),(size(tp_cov,1)));

for ith = 1:length(G_loop)
    % dYdx = (Wcoef(1,1:ns,G_loop(ith))*tp_cov((idx(G_loop(ith),1:ns)),3));
    % dYdy = (Wcoef(2,1:ns,G_loop(ith))*tp_cov((idx(G_loop(ith),1:ns)),3));
    % Cia = Wcoef(3,:,G_loop(ith)) + Wcoef(5,:,G_loop(ith));
    % Cib = dYdx*Wcoef(1,:,G_loop(ith));
    % Cic = dYdy*Wcoef(2,:,G_loop(ith));
    % C(G_loop(ith),(idx(G_loop(ith),:))) = Cia + Cib + Cic;
    C(G_loop(ith),(idx(G_loop(ith),:))) = Wcoef(3,:,G_loop(ith)) + Wcoef(5,:,G_loop(ith)); % Stationary
    f(G_loop(ith),1) = 0; % Sink/Source
end

for ith = 1:length(D_loop)
        C(D_loop(ith),D_loop(ith)) = 1;
        f(D_loop(ith),1) = tp_tvV(D_loop(ith),4);
end

for ith = 1:length(N_loop)
        C(N_loop(ith),(idx(N_loop(ith),:))) = (tp_tvV(N_loop(ith),5)).*Wcoef(1,:,N_loop(ith)) + (tp_tvV(N_loop(ith),6)).*Wcoef(2,:,N_loop(ith));
        f(N_loop(ith),1) = tp_tvV(N_loop(ith),4);
end

u=C\f;