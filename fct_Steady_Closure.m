function [alpha,beta] = fct_Steady_Closure(Wcoef,idx,ns,Y_mean,u_mean,u_var,g)
% V(x,y)
K_mean = exp(Y_mean);
V = zeros(length(Y_mean),1);
parfor ith =  1:length(Y_mean)
    V(ith,1) = -0.5.*K_mean(ith).*...
        ((Wcoef(3,1:ns,ith)*u_var((idx(ith,1:ns)),1)) + (Wcoef(5,1:ns,ith)*u_var((idx(ith,1:ns)),1)));
end

% beta(x,y)
beta = zeros(length(Y_mean),1);
parfor ith =  1:length(Y_mean)
    beta(ith,1) = K_mean(ith).*...
        ((Wcoef(3,1:ns,ith)*u_mean((idx(ith,1:ns)),1)) + (Wcoef(5,1:ns,ith)*u_mean((idx(ith,1:ns)),1)))+...
        g(ith,1); % no sink/source
end

% alpha
alpha = zeros(length(Y_mean),1);
for ith =  1:length(Y_mean)
    alpha(ith,1) = (K_mean(ith).*...
        (((Wcoef(1,1:ns,ith)*u_mean((idx(ith,1:ns)),1))^2 + (Wcoef(2,1:ns,ith)*u_mean((idx(ith,1:ns)),1))^2)) -...
        V(ith,1)) ./ u_var(ith,1); 
end


