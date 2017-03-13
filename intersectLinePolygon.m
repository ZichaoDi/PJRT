function intersects = intersectLinePolygon(line, poly)
tol = 1e-14;
% create the array of edges
N = size(poly, 1);
edges = [poly(1:N, :) poly([2:N 1], :)];
% compute intersections with supporting lines of polygon edges
supportLines = [edges(:, 1:2) edges(:, 3:4)-edges(:, 1:2)];
% intersects = intersectLines(line, supportLines, tol);
N1 = size(line, 1);
N2 = size(supportLines, 1);
N = max(N1, N2);

dx = bsxfun(@minus, supportLines(:,1), line(:,1));
dy = bsxfun(@minus, supportLines(:,2), line(:,2));

% indices of parallel lines
denom = line(:,3) .* supportLines(:,4) - supportLines(:,3) .* line(:,4);
par = abs(denom) < tol;

% indices of colinear lines
col = abs(dx .* line(:,4) - dy .* line(:,3)) < tol & par ;

% initialize result array
x0 = zeros(N, 1);
y0 = zeros(N, 1);

% initialize result for parallel lines
x0(col) = Inf;
y0(col) = Inf;
x0(par & ~col) = NaN;
y0(par & ~col) = NaN;

% in case all line couples are parallel, return
if all(par)
    point = [x0 y0];
    return;
end


%% Extract coordinates of itnersecting lines

% indices of intersecting lines
inds = ~par;

% extract base coordinates of first lines
if N1 > 1
    line = line(inds,:);
end
x1 =  line(:,1);
y1 =  line(:,2);
dx1 = line(:,3);
dy1 = line(:,4);

% extract base coordinates of second lines
if N2 > 1
    supportLines = supportLines(inds,:);
end
x2 =  supportLines(:,1);
y2 =  supportLines(:,2);
dx2 = supportLines(:,3);
dy2 = supportLines(:,4);

% re-compute coordinate differences of origin points
dx = bsxfun(@minus, supportLines(:,1), line(:,1));
dy = bsxfun(@minus, supportLines(:,2), line(:,2));


%% Compute intersection points

denom = denom(inds);
x0(inds) = (x2 .* dy2 .* dx1 - dy .* dx1 .* dx2 - x1 .* dy1 .* dx2) ./ denom ;
y0(inds) = (dx .* dy1 .* dy2 + y1 .* dx1 .* dy2 - y2 .* dx2 .* dy1) ./ denom ;

% concatenate result
intersects = [x0 y0];
% find edges that are not parallel to the input line
inds = find(isfinite(intersects(:, 1)));

% compute position of intersection points on corresponding lines
% pos = linePosition(intersects(inds, :), supportLines(inds, :), 'diag');

    np = length(inds);
    vx = supportLines(inds, 3);
    vy = supportLines(inds, 4);
    
    % difference of coordinates between point and line origins
    dx = intersects(inds, 1) - supportLines(inds, 1);
    dy = intersects(inds, 2) - supportLines(inds, 2);
% squared norm of direction vector, with a check of validity 
delta = vx .* vx + vy .* vy;
invalidLine = delta < eps;
delta(invalidLine) = 1; 

% compute position of points projected on the line, by using normalised dot
% product (NP-by-NE array) 
pos = bsxfun(@rdivide, bsxfun(@times, dx, vx) + bsxfun(@times, dy, vy), delta);

% ensure degenerated edges are correclty processed (consider the first
% vertex is the closest)
pos(:, invalidLine) = 0;
% and keep only intersection points located on edges
b = pos > -tol & pos < 1+tol;
inds = inds(b);
intersects = intersects(inds, :);

% remove multiple vertices (can occur for intersections located at polygon
% vertices)
[intersects, I, J] = unique(intersects, 'rows'); %#ok<NASGU>

