%% Initialization
clc; clear all; close all;

% Global Variables
N_SIM = 130;
LENGTH_R = 50000;
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
%scatter3(x,y,z) %

hold on;
location = [x(randomInt); y(randomInt); z(randomInt)];
scatter3(location(1),location(2),location(3),'r','x');
%histogram(y)
shg
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
    scatter3(x(result(b)),y(result(b)),z(result(b)),'b','filled');
end
shg
%% Second Point
indexPoint = randperm(length(result),1);    %selects random point from 1 to length of the result vector
secondPointElem = result(indexPoint);    %returns the corresponding element index from previous line
secondPoint = [x(secondPointElem);y(secondPointElem);z(secondPointElem)];   %gives point value indicated by element index
scatter3(secondPoint(1),secondPoint(2),secondPoint(3),'g','filled');    %plots random next point in green
shg
%% First Vector
firstVector = [secondPoint(1)-location(1);secondPoint(2)-location(2);secondPoint(3)-location(3)];
quiver3(location(1),location(2),location(3),firstVector(1),firstVector(2),firstVector(3))
totalVector = [firstVector(1),firstVector(2),firstVector(3),location(1),location(2),location(3)]; %running storage of vector and coordinates
scatter3(secondPoint(1),secondPoint(2),secondPoint(3),'g','filled');
shg