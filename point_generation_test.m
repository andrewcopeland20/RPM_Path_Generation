%% Initialization
clc; clear all; close all;

% Global Variables
N_SIM = 130;
LENGTH_R = 5;
DIM_BOX_SEARCH = 0.05;

%
r = randn(3,LENGTH_R); % Use a large n
%r = randi(5,3,LENGTH_R);

%%
r = bsxfun(@rdivide,r,sqrt(sum(r.^2,1)));
x = r(1,:);
y = r(2,:);
z = r(3,:);
randomInt = randi(LENGTH_R,1);

%figure
%scatter3(x,y,z) %

%hold on;
location = [x(randomInt); y(randomInt); z(randomInt)];
%scatter3(location(1),location(2),location(3),'r','x');

%% Distance Calculation

delta_r = r-location;
dist_r = sqrt(sum(delta_r.^2,1));
B = mink(dist_r,40);

for a = 1:LENGTH_R
    if abs(x(a)-location(1)) < DIM_BOX_SEARCH || abs(y(a)-location(2)) < DIM_BOX_SEARCH || abs(z(a)-location(3)) < DIM_BOX_SEARCH
        distance(a) = ((x(a)-location(1))^2 + (y(a)-location(2))^2 + (z(a)-location(3))^2)^(1/2);
    else
        distance(a) = 5;
    end
end

G = mink(distance,40);
%result = find(distance==G);
%scatter3(x(result),y(result),z(result),'b','filled');


r
distance
dist_r

