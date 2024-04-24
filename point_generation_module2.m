load('simSphere.mat')
solutions =[];
x = r(1,:);
y = r(2,:);
z = r(3,:);
averageGravity = 0.5;
zzz = 1;
vectorMaxLength = 30000;
while averageGravity > 0.0015
    %% Initialization
    tic
    %figure
    %scatter3(x,y,z)
    randomInt = randperm(length(r),1);
    hold on;
    randomPoint = [x(randomInt); y(randomInt); z(randomInt)];
    %scatter3(randomPoint(1),randomPoint(2),randomPoint(3),'r','x');
    location = randomPoint; %starts at a random point on the sphere
    %% Distance Calculation
    boxSearchDims = 0.05;
    for a = 1:length(x)
        if abs(x(a)-location(1)) < boxSearchDims || abs(y(a)-location(2)) < boxSearchDims || abs(z(a)-location(3)) < boxSearchDims
            distance(a) = ((x(a)-location(1))^2 + (y(a)-location(2))^2 + (z(a)-location(3))^2)^(1/2);
        else
            distance(a) = 5;
        end
    end
    closePoints = mink (distance,40);
    for b = 2:length(closePoints)
        result(b) = find(distance== closePoints(b));
        %scatter3(x(result(b)),y(result(b)),z(result(b)),'b','filled');
    end
    %% Second Point
    indexPoint = randperm(length(result),1);    %selects random point from 1 to length of the result vector
    secondPointElem = result(indexPoint);    %returns the corresponding element index from previous line
    secondPoint = [x(secondPointElem);y(secondPointElem);z(secondPointElem)];   %gives point value
    indicated by element index
    %scatter3(secondPoint(1),secondPoint(2),secondPoint(3),'g','fillsed');    %plots random next point in
    green
    %% First Vector
    firstVector = [secondPoint(1)-location(1);secondPoint(2)-location(2);secondPoint(3)-location(3)];
    %quiver3(location(1),location(2),location(3),firstVector(1),firstVector(2),firstVector(3))
    totalVector = [firstVector(1),firstVector(2),firstVector(3),location(1),location(2),location(3)]; %running
    storage of vector and coordinates
    %scatter3(secondPoint(1),secondPoint(2),secondPoint(3),'g','filled');
    %% Continuing Path
    %scatter3(x,y,z)
    %scatter3(randomPoint(1),randomPoint(2),randomPoint(3),'r','x');
    nextPoint = secondPoint; %starts off code
    previousPoint = secondPoint;
    lastVector = firstVector;
    %clear the plot and replot just the points and second point with first vector, research for new close points
    %this should eventually turn into a loop
    for sim = 2:vectorMaxLength
        %previous code
        %{
        for a = 1:length(x) 
        distance(a) = ((x(a)-nextPoint(1))^2 + (y(a)-nextPoint(2))^2 + (z(a)-nextPoint(3))^2)^(1/2); 
        end  
        %}
        %new proposal
        for a = 1:length(x)
            if abs(x(a)-nextPoint(1)) < boxSearchDims || abs(y(a)-nextPoint(2)) < boxSearchDims || abs(z(a)-nextPoint(3)) < boxSearchDims
                distance(a) = ((x(a)-nextPoint(1))^2 + (y(a)-nextPoint(2))^2 + (z(a)-nextPoint(3))^2)^(1/2);
            else
                distance(a) = 5;
            end
        end
        closePoints = mink (distance,75);
        for b = 1:length(closePoints)
            secondResult(b) = find(distance== closePoints(b));
            %scatter3(x(secondResult(b)),y(secondResult(b)),z(secondResult(b)),'b','filled');
        end
        for c = 1:length(closePoints)
            point(c,:) = [x(secondResult(c)),y(secondResult(c)),z(secondResult(c))];
            vectors(c,:) = [point(c,1)-nextPoint(1),point(c,2)-nextPoint(2),point(c,3)-nextPoint(3)];
            vectorRadian(c,:) = atan2(norm(cross(vectors(c,:),lastVector)), dot(vectors(c,:),lastVector));
            vectorAngles(c,1) = rad2deg(vectorRadian(c,1));

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
        cleanValidPoints = validPoints(any(validPoints,2),:); %this is where I would implement the weighting
        average
        %usedVector = (nextPoint(1),nextPoint(2),nextPoint(3)
        %totalVector(sim,4:6) = [nextPoint(1),nextPoint(2),nextPoint(3)];
        totalVector(sim,4:6) = [previousPoint(1),previousPoint(2),previousPoint(3)];
        indexPoint2 = randperm(size(cleanValidPoints,1),1);
        nextPoint(:) = cleanValidPoints(indexPoint2,:);
        %scatter3(nextPoint(1),nextPoint(2),nextPoint(3),'g','filled');
        totalVector(sim,1:3)= [nextPoint(1)-previousPoint(1),nextPoint(2)-previousPoint(2),nextPoint(3)
            previousPoint(3)]; %updates the total vector
        lastVector = [nextPoint(1)-previousPoint(1),nextPoint(2)-previousPoint(2),nextPoint(3)
            previousPoint(3)];
        %quiver3(previousPoint(1),previousPoint(2),previousPoint(3),totalVector(sim,1),totalVector(sim,2),totalVector(sim,3))
        previousPoint = nextPoint;
        %F = getframe
    end
    timer = toc;
    %aimshow(F.cdata)
    totalSolution(:,:,zzz) = totalVector;
    %quiver3(totalVector(:,4),totalVector(:,5),totalVector(:,6),totalVector(:,1),totalVector(:,2),totalVector(:,3))
    for a = 1:length(totalVector)
        u = [totalVector(a,4) totalVector(a,5) totalVector(a,6)];
        v = [0 0 -9.8];
        ThetaInDegrees = atan2d(norm(cross(u,v)),dot(u,v));
        Fref(a,1) = -9.81*cosd(ThetaInDegrees); %z component of gravity
        Fref(a,2) = 9.81*sind(ThetaInDegrees); %perpendicular component of gravity
        if u(2)<0
            Fref(a,2) = Fref(a,2)*-1;
        end

    end


    movingSum = [];

    movingSum(1,1) = sum(Fref(1:vectorMaxLength,1)./vectorMaxLength)/9.81;
    movingSum(1,2) = sum(Fref(1:vectorMaxLength,2)./vectorMaxLength)/9.81;


    averageGravity = sqrt(movingSum(1,1)^2 + movingSum(1,2)^2);

    disp(['Finished Solution ' num2str(zzz) ' in ' num2str(timer) ' seconds. Average gravity is '
        num2str(averageGravity) 'g.']);

    zzz = zzz + 1;

    remainder = mod(zzz,10);


end

%p0 = [1 2 3];
%p1 = [4 5 6];
%vectarrow(p0,p1)

%find vectors between the new point and surrounding points, select a point based on what has < a certain angle difference between the vectors


%% To do

%solution for range of microgravity
%slow motors down near the poles

%how do we convert vector output to motor controls
%project vector onto planes orthogonal to axes, the perpendicular distance
%crossed with magnitude of vector will give the speed