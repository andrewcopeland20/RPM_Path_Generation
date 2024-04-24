%% Initialization
clc; clear all; close all;

% Global Variables
N_SIM = 130;
LENGTH_R = 400;
DIM_BOX_SEARCH = 0.05;

%
r = randn(3,LENGTH_R); % Use a large n
r = r./sqrt(sum(r.^2,1));
%r = bsxfun(@rdivide,r,sqrt(sum(r.^2,1)));
x = r(1,:);
y = r(2,:);
z = r(3,:);
randomInt = randi(LENGTH_R,1);

%figure
scatter3(x,y,z,'blue','filled') %
grid off
axis equal
shading interp

hold on;
location = [x(randomInt); y(randomInt); z(randomInt)];
%scatter3(location(1),location(2),location(3),'r','x');
%histogram(y)
%shg
%% Distance Calculation
for a = 1:LENGTH_R
    if abs(x(a)-location(1)) < DIM_BOX_SEARCH || abs(y(a)-location(2)) < DIM_BOX_SEARCH || abs(z(a)-location(3)) < DIM_BOX_SEARCH
        distance(a) = ((x(a)-location(1))^2 + (y(a)-location(2))^2 + (z(a)-location(3))^2)^(1/2);
    else
        distance(a) = 5;
    end
end

% for a = i:LENGTH_R
%     distance(a) = ((x(a)-location(1))^2 + (y(a)-location(2))^2 + (z(a)-location(3))^2)^(1/2);
% end

%histogram(distance,'Normalization','probability')

closePoints = mink(distance,40);
for b = 2:length(closePoints)
    result(b) = find(distance==closePoints(b));
    %scatter3(x(result(b)),y(result(b)),z(result(b)),'b','filled');
end
%shg
%% Second Point
indexPoint = randperm(length(result),1); %selects random point from 1 to length of the result vector
secondPointElem = result(indexPoint);    %returns the corresponding element index from previous line
if secondPointElem == 0
    indexPoint = randperm(length(result),1); %selects random point from 1 to length of the result vector
    secondPointElem = result(indexPoint);    %returns the corresponding element index from previous line
else
    1
end
secondPoint = [x(secondPointElem);y(secondPointElem);z(secondPointElem)];   %gives point value indicated by element index
%scatter3(secondPoint(1),secondPoint(2),secondPoint(3),'g','filled');    %plots random next point in green
%shg
%% First Vector
firstVector = [secondPoint(1)-location(1);secondPoint(2)-location(2);secondPoint(3)-location(3)];
%quiver3(location(1),location(2),location(3),firstVector(1),firstVector(2),firstVector(3))
totalVector = [firstVector(1),firstVector(2),firstVector(3),location(1),location(2),location(3)]; %running storage of vector and coordinates
%scatter3(secondPoint(1),secondPoint(2),secondPoint(3),'g','filled');
%shg
%% Continuing Path
%scatter3(x,y,z)
%scatter3(location(1),location(2),location(3),'r','x');
axis equal
%set(gca,'XLim',[-1 1],'YLim',[-1 1],'ZLim',[-1 1]) %%%%%%%%%%%%%%%%%%%%%%%%%%%
nextPoint = secondPoint; %starts off code
previousPoint = secondPoint;
lastVector = firstVector;
%clear the plot and replot just the points and second point with first vector, research for new close points
%this should eventually turn into a loop
%shg
%% LOOOOOOP
for sim = 2:N_SIM
    %previous code
    %{
for a = 1:LENGTH_R 
distance(a) = ((x(a)-nextPoint(1))^2 + (y(a)-nextPoint(2))^2 + (z(a)-nextPoint(3))^2)^(1/2); 
end  
    %}
    %new proposal
    for a = 1:LENGTH_R
        if abs(x(a)-nextPoint(1)) < DIM_BOX_SEARCH || abs(y(a)-nextPoint(2)) < DIM_BOX_SEARCH || abs(z(a)-nextPoint(3)) < DIM_BOX_SEARCH
            distance(a) = ((x(a)-nextPoint(1))^2 + (y(a)-nextPoint(2))^2 + (z(a)-nextPoint(3))^2)^(1/2);
        else
            distance(a) = 5;
        end
    end
    
    closePoints = mink(distance,N_SIM/2);
    for b = 1:length(closePoints)
        secondResult(b) = find(distance== closePoints(b));
        %scatter3(x(secondResult(b)),y(secondResult(b)),z(secondResult(b)),'b','filled');
    end
    closepoints = mink(distance,10)
    for q = 1:length(closepoints)
            closeindeces(q) = find(distance==closepoints(q));
            closedistance(q,sim-1) = distance(closeindeces(q));
    end
    
    for c = 1:length(closePoints)
        point(c,:) = [x(secondResult(c)),y(secondResult(c)),z(secondResult(c))];
        vectors(c,:) = [point(c,1)-nextPoint(1),point(c,2)-nextPoint(2),point(c,3)-nextPoint(3)];
        vectorRadian(c,:) = atan2(norm(cross(vectors(c,:),lastVector)), dot(vectors(c,:),lastVector));
        vectorAngles(c,1) = rad2deg(vectorRadian(c,1)); % convert to degrees
        
        if vectorAngles(c,1) < 30 %invalidates all vectors that deviate by more than 30 degrees
            validPoints(c,:) = point(c,:);
            %scatter3(point(c,1),point(c,2),point(c,3),'b')
            
            if validPoints(c,1) == previousPoint(1) || validPoints(c,2) == previousPoint(2) || validPoints(c,3) == previousPoint(3)
                validPoints(c,:) = 0;
            end
        else
            %failingPoints(c,:) = point(c,:);
            validPoints(c,:) = 0;
        end
        %validIndexes(c,1) = vectorAngles(c,1) < 30;
        %scatter3(x(secondResult(b)),y(secondResult(b)),z(secondResult(b)),'b','filled');
    end

    cleanValidPoints = validPoints(any(validPoints,2),:); %this is where I would implement the weighting average
    
    %usedVector = (nextPoint(1),nextPoint(2),nextPoint(3))
    
    totalVector(sim,4:6) = [previousPoint(1),previousPoint(2),previousPoint(3)];
    
    indexPoint2 = randperm(size(cleanValidPoints,1),1);
    nextPoint(:) = cleanValidPoints(indexPoint2,:);
    %scatter3(nextPoint(1),nextPoint(2),nextPoint(3),'g','filled');
    totalVector(sim,1:3) = [nextPoint(1)-previousPoint(1),nextPoint(2)-previousPoint(2),nextPoint(3)-previousPoint(3)]; %updates the total vector
    lastVector = [nextPoint(1)-previousPoint(1),nextPoint(2)-previousPoint(2),nextPoint(3)-previousPoint(3)];
    %quiver3(previousPoint(1),previousPoint(2),previousPoint(3),totalVector(sim,1),totalVector(sim,2),totalVector(sim,3))
    
    previousPoint = nextPoint;
    
    F = getframe; % used for capturing current frame which is a snapshot of the current axes
end

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