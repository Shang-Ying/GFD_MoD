function xy = node_drop_2d_ctps_tol(box, bdy, ninit, dotmax, radius, ctps, btol, ptol, vargin)
% NODE_DROP_2D generates quasi-uniform nodes in a 2-D bounded box with 
% spacing between nodes specified by an exclusion radius function
%
% Copyright (C) 2019 Kiera van der Sande
%
% --- Input parameters ---
%   box         Size of box be filled by nodes; [xmin,xmax,ymin,ymax]
%   ninit       Number of PDP entries: [xinit]
%   dotmax      Upper bound for number of nodes to place
%   radius      The function radius(xy) provides exclusion radius to be used
%               at location (x,y).
%   vargin      Additional parameters for exclusion radius function
%
% --- Output parameter ---
%   xy          Array xy(:,2) with the generated node locations

if nargin < 6
    vargin = [];
end

dotnr   = 0;                            % Counter for the placed nodes
rng(0);                                 % Initialize random number generator
pdp     = [linspace(box(1),box(2),ninit)',box(3)*ones(ninit,1)]; % Array to hold PDPs
% r = radius(xy(dotnr,:),ctps); 

% Exclusion radius for bottom nodes
pdp(:,2) = pdp(:,2)+ 0.001*rand(ninit,1);  % Add random perturbation
xy      = zeros(dotmax,2);              % Array to store produced node locations
nodeindices = zeros(length(pdp),1);     % Array of pointers to the produced node locations 
excessheight = 0.1;                     % Percentage of the height to go over
dx      = pdp(2,1)-pdp(1,1);            % Grid size
[ym,i]  = min(pdp(:,2));                % Locate PDP with lowest y-coordinate

xy_new = pdp(i,:);
mdis = 0;

while ym <= (1+excessheight)*box(4) && dotnr < dotmax 
    % --- Add new node to generated nodes
    dotnr = dotnr + 1;                  
    xy(dotnr,:) = xy_new;             
    nodeindices(i) = dotnr;             
    
    % r = radius_2d_ctps(radius, mdis);  
    r = radius(xy(dotnr,:),ctps,bdy);
    
    % --- Find PDPs inside the new circle
    ileft  = i - floor(r/dx);
    ileft = max(1,ileft);
    iright = i + floor(r/dx);
    iright = min(ninit,iright);

    % --- Update heights of PDPs within radius 
    pdp(ileft:iright,2) = max([pdp(ileft:iright,2), sqrt(r^2-(pdp(ileft:iright,1)-pdp(i,1)).^2)+pdp(i,2)],[],2);
    
    % --- Identify next node location as a local minimum of the PDPs
    if pdp(ileft,2)<pdp(iright,2)
        i = ileft;
    else
        i = iright;
    end
    
    searchr = min(2*ceil(r/dx),floor(ninit/2)-1);     
    
    while 1
        % Wrap around if a boundary is reached
        if i-searchr < 1 
            xsearch = [ninit+i-searchr:ninit,1:i+searchr];
        elseif i+searchr > ninit
            xsearch = [i-searchr:ninit,1:i+searchr-ninit];
        else
            xsearch = i-searchr:i+searchr;
        end
        [ym,ix] = min(pdp(xsearch,2));
        i = xsearch(ix);
        
        % Stop once a local min has been found within the search radius
        if ix > searchr/2 && ix < length(xsearch)-searchr/2
            break
        end
    end
    % calculate distance between the new node and the ctps group
    xy_new = [pdp(i,:)];
end                                     

xy = xy(1:dotnr,:);                     % Remove unused entries in array xy
xy = xy(xy(:,2)<=box(4),:);             % Remove any nodes placed above the box

% check all ctps are in the generated xyz
% if yoou want some fixed nodes 
if ~isempty(ctps)
    xy=setdiff(xy,ctps,'rows'); 
end

pfix=unique(ctps,'row','stable'); % REMOVE DUPLICATE ONES
% p=[pfix; p];                                      
% ib  = find(min(pdist2(bdy,xy))< btol);
% xy(ib,:) = [];
dists = zeros(size(xy,1), 1);
for i = 1:size(xy,1)
    point = xy(i,:);
    minDist = inf;
    
    % Loop through each segment of the boundary
    for j = 1:size(bdy,1)
        % Define the line segment (wrap around for the last point)
        p1 = bdy(j,:);
        p2 = bdy(mod(j,size(bdy,1))+1,:);
        
        % Vector calculations for point-to-segment distance
        v = p1;
        w = p2;
        l2 = sum((w-v).^2);
        
        if l2 == 0
            segDist = sqrt(sum((point-v).^2));
        else
            t = max(0, min(1, sum((point-v).*(w-v))/l2));
            projection = v + t * (w-v);
            segDist = sqrt(sum((point-projection).^2));
        end
        
        minDist = min(minDist, segDist);
    end
    
    dists(i) = minDist;
end

% Remove points that are too close to the boundary
ib = find(dists < btol);
xy(ib,:) = [];
ib  = find(min(pdist2(ctps,xy))< ptol);
xy(ib,:) = [];
xy=[pfix; xy];

