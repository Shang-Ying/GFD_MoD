function   r = radius_3d_ctps_RH(radius, mdis) 
% % Return exclusion radius at location (x,y,z)
% % scale = mean([box(2)-box(1), box(4)-box(3), box(6)-box(5)]);
% % r = radius + 0.05*(mdis^1.5);
% 
% % r = radius + (((1./(1+exp(-(mdis-radius)/radius)))-0.5).*2).*radius;
% 
% % % s-shaped sigmoid function
% % x = -5:0.01:5;
% % y = 1./(1+exp(-x));
% % plot(x,y)
% 
% % example
% % radius = 5;
% % mdis = 0:0.01:30;
r = radius + (((1./(1+exp(-(mdis-radius)/radius)))-0.5).*2).*radius;
% % figure
% % plot(mdis,r)

% r = radius;