% Number of points
N = 100;

% Generate random points on the sphere with uniform distribution
theta = 2 * pi * rand(N, 1);
phi = acos(2 * rand(N, 1) - 1);
x = sin(phi) .* cos(theta);
y = sin(phi) .* sin(theta);
z = cos(phi);

% Points matrix
points = [x, y, z];

% Create a nearest neighbor searcher using a k-d tree
ns = createns(points, 'Distance', 'euclidean');

% Example: Find the 10 closest points to the first point in the list
[idx, dist] = knnsearch(ns, points(1, :), 'K', 11);

% Display the indices and distances of the nearest neighbors
disp('Indices of nearest neighbors:');
disp(idx(2:end));  % Ignore the first index as it is the point itself

disp('Distances to nearest neighbors:');
disp(dist(2:end));
