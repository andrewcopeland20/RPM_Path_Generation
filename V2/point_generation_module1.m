%% Initialization
clc; clear all; close all;

% Global Variables
N_SIM = 130;            % number of points in the trajectory generation
LENGTH_R = 400;         % number of points in the distribution
DIM_BOX_SEARCH = 0.05;  % distance defining box region
MAX_ANGLE_DEG = 30;     % maximum allowed trajectory deviation angle in degrees

% Generate random vector of points
r = randn(3, LENGTH_R); % Use a large n
r = r ./ vecnorm(r);
% r = r ./ sqrt(sum(r.^2, 1));
%r = bsxfun(@rdivide,r,sqrt(sum(r.^2,1))); from original code, replaced with improved version

% Extract x, y, z coordinates from random vector
x = r(1,:);
y = r(2,:);
z = r(3,:);

% Choose random initial point from vector of random points
% randomInt = randi(LENGTH_R,1);

%figure 1 - Vizualisation of distribution of points on a sphere
scatter3(x, y, z, 'bo', 'filled') %
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
hold off;

%histogram(y)
%shg
%% Distance Calculation
% for a = 1:LENGTH_R
%     if abs(x(a)-location(1)) < DIM_BOX_SEARCH || abs(y(a)-location(2)) < DIM_BOX_SEARCH || abs(z(a)-location(3)) < DIM_BOX_SEARCH
%         distance(a) = ((x(a)-location(1))^2 + (y(a)-location(2))^2 + (z(a)-location(3))^2)^(1/2);
%     else
%         distance(a) = 5;
%     end
% end


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


% for a = i:LENGTH_R
%     distance(a) = ((x(a)-location(1))^2 + (y(a)-location(2))^2 + (z(a)-location(3))^2)^(1/2);
% end

%histogram(distance,'Normalization','probability')

% Finding and vizualizing the 10 closest points
non_zero_distances = distance > 0; % Logical indexing for non-zero distances
[closest_distances, closest_indices] = mink(distance(non_zero_distances), 10);

% Extract coordinates of the closest points
closest_points_x = x(closest_indices);
closest_points_y = y(closest_indices);
closest_points_z = z(closest_indices);

% Visualize these points
scatter3(closest_points_x, closest_points_y, closest_points_z, 'filled');
title('10 Closest Points');
xlabel('X');
ylabel('Y');
zlabel('Z');

% Store indices of distances corresponding to first 10 points to match
% original code
result = closest_indices(2:end); % NOTE: Figure out why it starts from 2 not 1!!!

% closePoints = mink(distance,40);
% for b = 2:length(closePoints)
%     result(b) = find(distance==closePoints(b));
%     %scatter3(x(result(b)),y(result(b)),z(result(b)),'b','filled');
% end
% 

%shg
%% Second Point
% indexPoint = randperm(length(result),1); %selects random point from 1 to length of the result vector
% secondPointElem = result(indexPoint);    %returns the corresponding element index from previous line
% if secondPointElem == 0
%     indexPoint = randperm(length(result),1); %selects random point from 1 to length of the result vector
%     secondPointElem = result(indexPoint);    %returns the corresponding element index from previous line
% else
%     1
% end
% secondPoint = [x(secondPointElem);y(secondPointElem);z(secondPointElem)];   %gives point value indicated by element index

% Generate random index from 1 to number of elements in array of closest
% point, and define the second point from this:

random_index = randi(length(closest_points_x));
second_point = [closest_points_x(random_index), closest_points_y(random_index), closest_points_z(random_index)];

%scatter3(second_point(1),second_point(2),second_point(3),'g','filled');    %plots random next point in green
%shg
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
previous_point = second_point;
last_vector = first_vector;
%clear the plot and replot just the points and second point with first vector, research for new close points
%this should eventually turn into a loop
%shg
%% LOOOOOOP
% for sim = 2:N_SIM
%     for a = 1:LENGTH_R
%         if abs(x(a)-next_point(1)) < DIM_BOX_SEARCH || abs(y(a)-next_point(2)) < DIM_BOX_SEARCH || abs(z(a)-next_point(3)) < DIM_BOX_SEARCH
%             distance(a) = ((x(a)-next_point(1))^2 + (y(a)-next_point(2))^2 + (z(a)-next_point(3))^2)^(1/2);
%         else
%             distance(a) = 5;
%         end
%     end
%     
%     closePoints = mink(distance,N_SIM/2);
%     for b = 1:length(closePoints)
%         secondResult(b) = find(distance== closePoints(b));
%         %scatter3(x(secondResult(b)),y(secondResult(b)),z(secondResult(b)),'b','filled');
%     end
%     closepoints = mink(distance,10)
%     for q = 1:length(closepoints)
%             closeindeces(q) = find(distance==closepoints(q));
%             closedistance(q,sim-1) = distance(closeindeces(q));
%     end
%     
%     for c = 1:length(closePoints)
%         point(c,:) = [x(secondResult(c)),y(secondResult(c)),z(secondResult(c))];
%         vectors(c,:) = [point(c,1)-next_point(1),point(c,2)-next_point(2),point(c,3)-next_point(3)];
%         vectorRadian(c,:) = atan2(norm(cross(vectors(c,:),last_vector)), dot(vectors(c,:),last_vector));
%         vectorAngles(c,1) = rad2deg(vectorRadian(c,1)); % convert to degrees
%         
%         if vectorAngles(c,1) < 30 %invalidates all vectors that deviate by more than 30 degrees
%             validPoints(c,:) = point(c,:);
%             %scatter3(point(c,1),point(c,2),point(c,3),'b')
%             
%             if validPoints(c,1) == previous_point(1) || validPoints(c,2) == previous_point(2) || validPoints(c,3) == previous_point(3)
%                 validPoints(c,:) = 0;
%             end
%         else
%             %failingPoints(c,:) = point(c,:);
%             validPoints(c,:) = 0;
%         end
%         %validIndexes(c,1) = vectorAngles(c,1) < 30;
%         %scatter3(x(secondResult(b)),y(secondResult(b)),z(secondResult(b)),'b','filled');
%     end
% 
%     cleanValidPoints = validPoints(any(validPoints,2),:); %this is where I would implement the weighting average
%     
%     %usedVector = (nextPoint(1),nextPoint(2),nextPoint(3))
%     
%     total_vector(sim,4:6) = [previous_point(1),previous_point(2),previous_point(3)];
%     
%     indexPoint2 = randperm(size(cleanValidPoints,1),1);
%     next_point(:) = cleanValidPoints(indexPoint2,:);
%     %scatter3(nextPoint(1),nextPoint(2),nextPoint(3),'g','filled');
%     total_vector(sim,1:3) = [next_point(1)-previous_point(1),next_point(2)-previous_point(2),next_point(3)-previous_point(3)]; %updates the total vector
%     last_vector = [next_point(1)-previous_point(1),next_point(2)-previous_point(2),next_point(3)-previous_point(3)];
%     %quiver3(previousPoint(1),previousPoint(2),previousPoint(3),totalVector(sim,1),totalVector(sim,2),totalVector(sim,3))
%     
%     previous_point = next_point;
%     
%     F = getframe; % used for capturing current frame which is a snapshot of the current axes
% end

% VERSION 2 OF THE LOOP:
for sim = 1:N_SIM
    % Compute vectors to all possible next points
    vectors_to_next = r - previous_point;

    % Compute distances
    distances = vecnorm(vectors_to_next, 2, 1);

    % Apply box search condition
    in_box = distances < DIM_BOX_SEARCH;
    distances(~in_box) = 100;  % Set distances outside the box search to large integer

    % Calculate angles using vector operations
    cosine_angles = dot(vectors_to_next, last_vector) ./ (vecnorm(vectors_to_next, 2, 1) .* norm(last_vector));
    angles = acosd(cosine_angles);  % Convert cosine values to angles in degrees

    % Find indices of points where angle and distance conditions are met
    valid_points = (angles < MAX_ANGLE_DEG) & (distances < DIM_BOX_SEARCH);

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

    % Select a new point randomly from 10 closest valid points
    new_point_index = closest_valid_indices(randi(length(closest_valid_indices)));
    next_point = r(:, new_point_index);

    % Compute the vector to the new point
    vector_to_next = next_point - previous_point;

    % Update last_vector to the new movement direction
    last_vector = vector_to_next / norm(vector_to_next);

    % Store the trajectory information
    total_vector(sim, :) = [vector_to_next' previous_point'];

    % Update for next iteration
    previous_point = next_point;
end

% % Optionally plot the trajectory
% plot3(total_vector(:, 4), total_vector(:, 5), total_vector(:, 6), 'rx-');
% xlabel('X'); ylabel('Y'); zlabel('Z');
% title('Trajectory Plot');

%imshow(F.cdata) % cdata from getframe is a matrix containing the image data
%p0 = [1 2 3];
%p1 = [4 5 6];
%vectarrow(p0,p1)
%find vectors between the new point and surrounding points, select a point based on what has < a certain angle difference between the vectors

%% histogram of point distance distribution
% closedistance = closedistance(any(closedistance,2),:);
% histogram(closedistance,'Normalization','probability')
% ylabel('probability')
% xlabel('distance')
% 
% i = 1;
% m = 0;
% for i = 1:length(closedistance)
%     m(i) = mean(closedistance(:,i))
% end
% 
% normplot(closedistance)


%% To do
%see if we can search by an area rather than calculating distance to
%points-- do a comparison of point distances first, i.e. check to see if x
%is within a certain distance
%this will allow us to increase the point density without exponentially
%increasing compute time
%should we order the points on the sphere first?