function [alpha,beta] = fct_Steady_Closure_FDM(tp_Y,tp_tv,u_mean,u_var,g,G_loop,D_loop,Nxm_loop,NxM_loop,Nym_loop,NyM_loop,n_right,n_left,n_upper,n_lower)

G_In = (tp_tv(:,3) == 0);
Bd_xm = (tp_tv(:,3) ~= 0 & tp_tv(:,1) == 0 & tp_tv(:,2) ~= 0 & tp_tv(:,2) ~= 1);
Bd_xM = (tp_tv(:,3) ~= 0 & tp_tv(:,1) == 1 & tp_tv(:,2) ~= 0 & tp_tv(:,2) ~= 1);
Bd_ym = (tp_tv(:,3) ~= 0 & tp_tv(:,2) == 0 & tp_tv(:,1) ~= 0 & tp_tv(:,1) ~= 1);
Bd_yM = (tp_tv(:,3) ~= 0 & tp_tv(:,2) == 1 & tp_tv(:,1) ~= 0 & tp_tv(:,1) ~= 1);

G_loop = find(G_In);
Nxm_loop = find(Bd_xm);
NxM_loop = find(Bd_xM);
Nym_loop = find(Bd_ym);
NyM_loop = find(Bd_yM);

d = tp_tv(1,2)-tp_tv(2,2);

K_mean = exp(tp_Y(:,3));

% V(x,y)
V = zeros(size(tp_Y,1),1);
% beta(x,y)
beta = zeros(size(tp_Y,1),1);
% alpha
alpha = zeros(size(tp_Y,1),1);
for ith = 1:length(G_loop)
    n_right = find(ismembertol(tp_tv(:,1:2), [tp_tv(G_loop(ith),1)+d tp_tv(G_loop(ith),2)],'ByRows',true));
    n_left = find(ismembertol(tp_tv(:,1:2), [tp_tv(G_loop(ith),1)-d tp_tv(G_loop(ith),2)],'ByRows',true));
    n_upper = find(ismembertol(tp_tv(:,1:2), [tp_tv(G_loop(ith),1) tp_tv(G_loop(ith),2)+d],'ByRows',true));
    n_lower = find(ismembertol(tp_tv(:,1:2), [tp_tv(G_loop(ith),1) tp_tv(G_loop(ith),2)-d],'ByRows',true));
    V(G_loop(ith),1) = -0.5.*K_mean(G_loop(ith)).*...
        ((u_var(n_right,1)+u_var(n_left,1)+u_var(n_upper,1)+u_var(n_lower,1)-4*u_var(G_loop(ith),1))/(d^2));
    beta(G_loop(ith),1) = K_mean(G_loop(ith)).*...
        ((u_mean(n_right,1)+u_mean(n_left,1)+u_mean(n_upper,1)+u_mean(n_lower,1)-4*u_mean(G_loop(ith),1))/(d^2))+...
        g(G_loop(ith),1); % no sink/source
    alpha(G_loop(ith),1) = ((K_mean(G_loop(ith)).*...
        ((((u_mean(n_right,1)-u_mean(n_left,1))/(2*d))^2) + (((u_mean(n_upper,1)-u_mean(n_lower,1))/(2*d))^2)) -...
        V(G_loop(ith),1))) ./ u_var(G_loop(ith),1); 
end

for ith = 1:length(Nxm_loop)
    n_right = find(ismembertol(tp_tv(:,1:2), [tp_tv(Nxm_loop(ith),1)+d tp_tv(Nxm_loop(ith),2)],'ByRows',true));
    n_upper = find(ismembertol(tp_tv(:,1:2), [tp_tv(Nxm_loop(ith),1) tp_tv(Nxm_loop(ith),2)+d],'ByRows',true));
    n_lower = find(ismembertol(tp_tv(:,1:2), [tp_tv(Nxm_loop(ith),1) tp_tv(Nxm_loop(ith),2)-d],'ByRows',true));
    V(Nxm_loop(ith),1) = -0.5.*K_mean(Nxm_loop(ith)).*...
        ((2*u_var(n_right,1)+u_var(n_upper,1)+u_var(n_lower,1)-4*u_var(Nxm_loop(ith),1))/(d^2));
    beta(Nxm_loop(ith),1) = K_mean(Nxm_loop(ith)).*...
        ((2*u_mean(n_right,1)+u_mean(n_upper,1)+u_mean(n_lower,1)-4*u_mean(Nxm_loop(ith),1))/(d^2))+...
        g(Nxm_loop(ith),1); % no sink/source
    alpha(Nxm_loop(ith),1) = ((K_mean(Nxm_loop(ith)).*...
        ((((u_mean(n_right,1)-u_mean(Nxm_loop(ith),1))/(d))^2) + (((u_mean(n_upper,1)-u_mean(n_lower,1))/(2*d))^2)) -...
        V(Nxm_loop(ith),1))) ./ u_var(Nxm_loop(ith),1); 
end

for ith = 1:length(NxM_loop)
    n_left = find(ismembertol(tp_tv(:,1:2), [tp_tv(NxM_loop(ith),1)-d tp_tv(NxM_loop(ith),2)],'ByRows',true));
    n_upper = find(ismembertol(tp_tv(:,1:2), [tp_tv(NxM_loop(ith),1) tp_tv(NxM_loop(ith),2)+d],'ByRows',true));
    n_lower = find(ismembertol(tp_tv(:,1:2), [tp_tv(NxM_loop(ith),1) tp_tv(NxM_loop(ith),2)-d],'ByRows',true));
    V(NxM_loop(ith),1) = -0.5.*K_mean(NxM_loop(ith)).*...
        ((2*u_var(n_left,1)+u_var(n_upper,1)+u_var(n_lower,1)-4*u_var(NxM_loop(ith),1))/(d^2));
    beta(NxM_loop(ith),1) = K_mean(NxM_loop(ith)).*...
        ((2*u_mean(n_left,1)+u_mean(n_upper,1)+u_mean(n_lower,1)-4*u_mean(NxM_loop(ith),1))/(d^2))+...
        g(NxM_loop(ith),1); % no sink/source
    alpha(NxM_loop(ith),1) = ((K_mean(NxM_loop(ith)).*...
        ((((u_mean(NxM_loop(ith),1)-u_mean(n_right,1))/(d))^2) + (((u_mean(n_upper,1)-u_mean(n_lower,1))/(2*d))^2)) -...
        V(NxM_loop(ith),1))) ./ u_var(NxM_loop(ith),1);
end

for ith = 1:length(Nym_loop)
    n_upper = find(ismembertol(tp_tv(:,1:2), [tp_tv(Nym_loop(ith),1) tp_tv(Nym_loop(ith),2)+d],'ByRows',true));
    n_left = find(ismembertol(tp_tv(:,1:2), [tp_tv(Nym_loop(ith),1)-d tp_tv(Nym_loop(ith),2)],'ByRows',true));
    n_right = find(ismembertol(tp_tv(:,1:2), [tp_tv(Nym_loop(ith),1)+d tp_tv(Nym_loop(ith),2)],'ByRows',true));
    V(Nym_loop(ith),1) = -0.5.*K_mean(Nym_loop(ith)).*...
        ((u_var(n_right,1)+u_var(n_left,1)+2*u_var(n_upper,1)-4*u_var(Nym_loop(ith),1))/(d^2));
    beta(Nym_loop(ith),1) = K_mean(Nym_loop(ith)).*...
        ((u_mean(n_right,1)+u_mean(n_left,1)+2*u_mean(n_upper,1)-4*u_mean(Nym_loop(ith),1))/(d^2))+...
        g(Nym_loop(ith),1); % no sink/source
    alpha(Nym_loop(ith),1) = ((K_mean(Nym_loop(ith)).*...
        ((((u_mean(n_right,1)-u_mean(n_left,1))/(2*d))^2) + (((u_mean(n_upper,1)-u_mean(Nym_loop(ith),1))/(d))^2)) -...
        V(Nym_loop(ith),1))) ./ u_var(Nym_loop(ith),1); 
end

for ith = 1:length(NyM_loop)
    n_lower = find(ismembertol(tp_tv(:,1:2), [tp_tv(NyM_loop(ith),1) tp_tv(NyM_loop(ith),2)-d],'ByRows',true));
    n_left = find(ismembertol(tp_tv(:,1:2), [tp_tv(NyM_loop(ith),1)-d tp_tv(NyM_loop(ith),2)],'ByRows',true));
    n_right = find(ismembertol(tp_tv(:,1:2), [tp_tv(NyM_loop(ith),1)+d tp_tv(NyM_loop(ith),2)],'ByRows',true));
    V(Nym_loop(ith),1) = -0.5.*K_mean(Nym_loop(ith)).*...
        ((u_var(n_right,1)+u_var(n_left,1)+2*u_var(n_lower,1)-4*u_var(Nym_loop(ith),1))/(d^2));
    beta(Nym_loop(ith),1) = K_mean(Nym_loop(ith)).*...
        ((u_mean(n_right,1)+u_mean(n_left,1)+2*u_mean(n_lower,1)-4*u_mean(Nym_loop(ith),1))/(d^2))+...
        g(Nym_loop(ith),1); % no sink/source
    alpha(Nym_loop(ith),1) = ((K_mean(Nym_loop(ith)).*...
        ((((u_mean(n_right,1)-u_mean(n_left,1))/(2*d))^2) + (((u_mean(Nym_loop(ith),1)-u_mean(n_lower,1))/(d))^2)) -...
        V(Nym_loop(ith),1))) ./ u_var(Nym_loop(ith),1); 
end


