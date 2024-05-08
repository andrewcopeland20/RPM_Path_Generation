%% Initialization
clc; clear; close all;

% Global Variables
N_SIM = 150;            % number of simulations
LENGTH_R = 500000;       % number of random points
DIM_BOX_SEARCH = 0.06;  % search box dimension
MAX_ANGLE_DEG = 15;     % maximum deviation angle in degrees
MAX_COS_ANGLE = cosd(MAX_ANGLE_DEG);
DS = 0.05;              % constant distance to travel along trajectory  

% Generate random unit vectors
r = randn(3, LENGTH_R);
r = r ./ vecnorm(r);

% Extract x, y, z coordinates from random vector
x = r(1,:);
y = r(2,:);
z = r(3,:);

% Choose random initial point from vector of random points
% randomInt = randi(LENGTH_R,1);

%figure 1 - Vizualisation of distribution of points on a sphere
figure(1)
scatter3(x(1:500:end), y(1:500:end), z(1:500:end), 'bo', 'filled') %
grid off
axis equal
shading interp
view(3)
hold on;


% Select initial location defined by the randomly chosen point and plot
% first_point = [x(randomInt); y(randomInt); z(randomInt)];
first_point = r(:, randi(LENGTH_R));
plot3(first_point(1), first_point(2), first_point(3), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
plot3(first_point(1), first_point(2), first_point(3), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
%hold off;

%%
% V2: Distance calculation using logical indexing
% Preallocate the distance array with default value
distance = 5 * ones(1, LENGTH_R);

% Logical condition to check proximity within DIM_BOX_SEARCH
in_box = abs(x - first_point(1)) < DIM_BOX_SEARCH | ...
         abs(y - first_point(2)) < DIM_BOX_SEARCH | ...
         abs(z - first_point(3)) < DIM_BOX_SEARCH;

% Calculate distance only for points within the box region
distance(in_box) = sqrt((x(in_box) - first_point(1)).^2 + ...
                        (y(in_box) - first_point(2)).^2 + ...
                        (z(in_box) - first_point(3)).^2);

% Finding and vizualizing the 10 closest points
non_zero_distances = distance > 0; % Logical indexing for non-zero distances
[closest_distances, closest_indices] = mink(distance(non_zero_distances), 10);

% Extract coordinates of the closest points
closest_points_x = x(closest_indices);
closest_points_y = y(closest_indices);
closest_points_z = z(closest_indices);

% Store indices of distances corresponding to first 10 points to match
% original code
result = closest_indices; 

%% Second Point

% Generate random index from 1 to number of elements in array of closest
% point, and define the second point from this:

random_index = randi(length(closest_indices));
second_point = [closest_points_x(random_index); closest_points_y(random_index); closest_points_z(random_index)];

%% First Vector
total_vector = zeros(N_SIM, 6); % Pre-allocate total_vector

first_vector = [second_point(1)-first_point(1);second_point(2)-first_point(2);second_point(3)-first_point(3)];
%quiver3(location(1),location(2),location(3),firstVector(1),firstVector(2),firstVector(3))
total_vector = [first_vector(1),first_vector(2),first_vector(3),first_point(1),first_point(2),first_point(3)]; %running storage of vector and coordinates
%scatter3(secondPoint(1),secondPoint(2),secondPoint(3),'g','filled');
%shg
%% Continuing Path
%scatter3(x,y,z)
%scatter3(location(1),location(2),location(3),'r','x');
%set(gca,'XLim',[-1 1],'YLim',[-1 1],'ZLim',[-1 1]) %%%%%%%%%%%%%%%%%%%%%%%%%%%

next_point = second_point; %starts off code
previous_point = first_point;
last_vector = first_vector;

%% Simulation Loop

for sim = 1:N_SIM
    % Compute distances to all other points using vecnorm
    vectors_to_next = r - previous_point;
    distances = vecnorm(vectors_to_next, 2, 1);
    
    % Find points within search box and under max angle
    in_box = distances < DIM_BOX_SEARCH;
    distances(~in_box) = inf;  % Exclude points outside the box

    % Find indices of the points within this distance
    interim_distances = distances(distances~=inf);
    [~, interim_valid_indices] = mink(interim_distances, 50);  % Choose 50 closest points from this subset

    % Extract indices from closest 50 points based on distance criteria and
    % convert to actual points
    interim_closest_indices = find(interim_distances);
    interim_closest_valid_indices = interim_closest_indices(interim_valid_indices);
    interim_valid_points = r(:, interim_closest_valid_indices);

    % plot first step of valid points based on distance criteria
    scatter3(interim_valid_points(1,:), interim_valid_points(2,:), interim_valid_points(3,:), 'yo', 'filled');

    cosine_angles = dot(vectors_to_next, repmat(last_vector, 1, LENGTH_R), 1) ./ (distances .* norm(last_vector));
    %cosine_angles = max(min(cosine_angles, 1), -1);  % Clamp cosine values to valid range
    %angles = acosd(cosine_angles);

    valid_points = (cosine_angles <= MAX_COS_ANGLE) & (distances < DIM_BOX_SEARCH);

    % Break if no valid points are found
    if all(~valid_points)
        disp('No valid points found under angle and distance constraints.');
        break;
    end

    % Find indices of the closest valid points
    valid_distances = distances(valid_points);
    [~, valid_indices] = mink(valid_distances, 10);  % Get 10 closest valid points

    % Extract actual indices from valid points indexing
    closest_indices = find(valid_points);
    closest_valid_indices = closest_indices(valid_indices);

    % Extract actual valid points and plot
    interim_new_points = r(:, closest_valid_indices)
    scatter3(interim_new_points(1,:), interim_new_points(2,:), interim_new_points(3,:), 'mo', 'filled')

    % Select a new point randomly from 10 closest valid points
    new_point_index = closest_valid_indices(randi(length(closest_valid_indices)));
    new_point = r(:, new_point_index);

    % Set up new vector to interim chosen point and normalize
    vector_to_new = new_point - previous_point;
    vector_to_new_norm = vector_to_new / norm(vector_to_new);

    % Find next point along trajectory and project back onto the sphere
    next_point = previous_point + vector_to_new_norm*DS;
    next_point = next_point ./ norm(next_point);

    % Set up vector to next
    vector_to_next = next_point - previous_point;

    % Set up next point and last vector based on this chosen next point
    % which is about to become the last point
    last_vector = vector_to_next;

    % Visualization
    plot3([previous_point(1), next_point(1)], [previous_point(2), next_point(2)], [previous_point(3), next_point(3)], 'g-');  % Green line for trajectory
    scatter3(next_point(1), next_point(2), next_point(3), 'go', 'filled');
    drawnow;  % Update the plot

    if ~ishandle(gcf)  % Stop if the figure is closed
        break;
    end

    previous_point = next_point;
    pause(0.1);  % Visual effect

end
